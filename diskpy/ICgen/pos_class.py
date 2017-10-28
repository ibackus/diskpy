# -*- coding: utf-8 -*-
"""
Defines a function to randomly generate particle positions according to 
the desired surface density profile (sigma vs r) and the vertical profile
(rho vs r,z).

Created on Mon Jan 27 18:48:04 2014

@author: ibackus
"""

__version__ = "$Revision: 1 $"
# $Source$

__iversion__ = int(filter(str.isdigit,__version__))

# External packages
import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
from scipy.integrate import cumtrapz
import os

# ICgen packages
from diskpy import global_settings
from diskpy.utils import strip_units
from diskpy.pdmath import interp1dunits
import ICgen_utils

    
class pos:
    """
    position class.  Generates particle positions from rho and sigma
    
    USAGE:
    # method = 'grid' or 'random'    
    pos = pos_class.pos(ICobj, method)
    
    ICobj should be an initial conditions object (ICgen.IC) with rho already
    calculated.

    """
    
    def __init__(self, ICobj, method = None, generate=True, seed=None):
        
        self._seed = seed
        # Set version
        self.__version__ = __iversion__
        # Link to parent initial conditions object
        self._parent = ICobj
        # Check that sigma and rho have been generated
        if not hasattr(ICobj, 'rho'):
            
            raise NameError,'rho could not be found in the IC object'
        
        if not hasattr(ICobj,'sigma'):
            
            raise NameError,'sigma could not be found in the IC object'
            
        if method == None:
        
            self.method = ICobj.settings.pos_gen.method
            
        else:
            
            self.method = method
            # Update settings in ICobj
            ICobj.settings.pos_gen.method = method
            
        self.nParticles = ICobj.settings.pos_gen.nParticles
        print 'Generating {0} particle positions using method: {1}'.format(\
        self.nParticles, self.method)
        
        # Generate positions
        self._generate_r()
        self.xyz = SimArray(np.zeros([self.nParticles, 3], dtype=np.float32), self.r.units)
        self._generate_z()
        self._generate_theta()
        self._cartesian_pos()
        
        # To save on memory, delete theta.  It can be re-calculated later
        # if absolutely needed
        del self.theta
        
    def __getstate__(self):
        """
        This is required to make the object pickle-able
        """
        
        # Define a dictionary containing everything needed.  
        # Ignore self.parent
        state = self.__dict__.copy()
        state.pop('_parent', None)
        
        # Now handle the possibly large arrays (too large to pickle)
        for key,val in state.iteritems():
            
            if isinstance(val, np.ndarray):
                
                state[key] = ICgen_utils.listify(val, 1001)
                
        return state
        
    def __setstate__(self, d):
        """
        This is required to make the object un-pickleable
        """
        for key, val in d.iteritems():
            
            if isinstance(val, ICgen_utils.larray):
                
                d[key] = val.delistify()
                
        self.__dict__ = d
        
    def _generate_r_glass(self):
        r"""
        Generates x-y positions and z-values mapped to -1 to 1 by tiling a
        small peridic cube glass and mapping to physical space.
        
        This is done by generating a uniformly distributed glass cylinder in
        what I'll call CDF space, for lack of a better word.  By using the
        inverse radial and vertical CDFs we can map points on that cylinder to
        points in real space which will be distributed according to the PDF
        (i.e. density) in real space.
        
        Since the ICs are vaguely cylindrical, this can be done without 
        distorting the glass too-much (locally).  The final result is
        significantly more glass-like than using other methods.
        
        The only trick is to approximate the height of the cylinder in CDF
        space in a reasonable way.  One method that works is to just
        take the average <H/R>, mass weighted over the disk (For a unit radius
        cylinder).
        
        I use a slightly more convoluted method which seems to work pretty well,
        but either way should be good.
        
        To map from :math:`\tilde{z}, \tilde{r}` to :math:`z, r`, (i.e.
        CDF to physical space) we can do the following:
        
        .. math::
            
            z = \mathrm{CDF}_z^{-1}[\tilde{z}/\tilde{z}_{max}]
            
            r = \mathrm{CDF}_r^{-1}[(\tilde{r}/\tilde{r}_{max})^2]
            
        where the square of :math:`(\tilde{r}/\tilde{r}_{max})` comes from the
        fact that the particles are uniformly distributed along x-y in our
        CDF space, not along r, i.e. 
        :math:`\mathrm{CDF}_{\tilde{r}} = (\tilde{r}/\tilde{r}_{max})^2`
        """
        cdf_inv_r = self._parent.sigma.cdf_inv
        # Load the small glass cube template
        glassfile = os.path.join(global_settings['misc']['bin-dir'],
                                 global_settings['misc']['icgen-glass-file'])
        sn = pynbody.load(glassfile)
        n_glass_cube = len(sn)
        nParticles = self.nParticles
        rmax = 1.
        zmax = cdf_space_height(self._parent, rmax)
        L_cube = ((2 * np.pi*zmax * rmax**2 * n_glass_cube)/nParticles)**(1./3)
        print 'zmax:', zmax
        print 'L_cube:', L_cube
        # Re-scale the cube
        sn.g['pos'] *= L_cube
        # Round up a little extra to make sure num boxes is not zero
        nbox_r = int(rmax/L_cube+.6)
        # Generate a grid of box centers for a single z-layer
        x = np.arange(-nbox_r, nbox_r+1) * L_cube
        y = x
        z = np.array([0.5*L_cube]) # First z-layer location
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        centers = np.stack([X.flatten(), Y.flatten(), Z.flatten()],axis=1)
        
        # Filter out cubes that definitely don't overlap the cylinder at all
        cube_r = np.sqrt(centers[:,0]**2 + centers[:,1]**2)
        keep = (cube_r - L_cube/np.sqrt(2)) <= rmax
        centers = centers[keep]
        # How many complete layers
        n_layers = int(2*zmax/L_cube)
        
        # Get a sub-section of the input glass blox to use for the top-layer
        glass_z = sn.g['z'].copy()
        glass_z.units = None
        top_z_middle = L_cube * (n_layers + 0.5)
        keep = (glass_z + top_z_middle) <= float(2 * zmax)
        sn_top = sn.g[keep]
        n_glass_top = len(sn_top)
        
        # Total maximum number of particles
        n_per_layer = len(centers) * n_glass_cube
        n_top_layer = len(centers) * n_glass_top
        npMax = n_per_layer*n_layers + n_top_layer
        
        # Tile the glass, layer by layer
        pos = np.zeros([npMax, 3])
        # Handle the complete layers
        for iLayer in range(n_layers):
            i0 = iLayer * n_per_layer
            layer = pos[i0:i0 + n_per_layer]
            for i, center in enumerate(centers):
                # Fill up this layer by tiling
                ilo = i*n_glass_cube
                ihi = (i+1)*n_glass_cube
                layer[ilo:ihi] = sn.g['pos'] + center[None,:]
            # Shift z-center up by cube length for next layer
            centers[:, 2] += L_cube
        # Handle the top layer
        layer = pos[n_layers*n_per_layer:]
        for i, center in enumerate(centers):
            ilo = i*len(sn_top)
            ihi = (i+1)*len(sn_top)
            layer[ilo:ihi] = sn_top['pos'] + center[None,:]
            
        assert(np.all(pos[:,2] <= 2*zmax))
        pos[:, 2] -= zmax
        # Find radius that would include nParticles
        r = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
        rsort = np.sort(r)
        rcut = 0.5*(rsort[nParticles-1]+rsort[nParticles])
        # If all the scaling is done right, rcut should be nearly 1
        print "(rmax, rcut)", rmax, rcut
        if abs(rcut-rmax)/rmax > 1e-2:
            print 'Warning: rcut not good...might try differet nParticles'
        # Assign output
        mask = r < rcut
        pos = pos[mask]
        self.theta = np.arctan2(pos[:,1], pos[:,0])
        # Store z/zmax which will be used by the inverse CDF 
        # later to map to z in physical space
        self.m_z = pos[:,2]/zmax
        r = (r[mask]/rcut)
        # Here we use r^2 since our particles are not uniformly distributed
        # along r
        self.r = cdf_inv_r(r**2).astype(np.float32)
        
    def _generate_r_grid(self):
        """
        Generates positions on a spiral 'grid' in x-y following the method
        of Cartwright, Stamatellos & Whitworth 2009
        
        r positions are generated directly from the inverse radial CDF
        """
        cdf_inv_r = self._parent.sigma.cdf_inv
        # Generate linearly increasing values of m, using 2 more than
        # necessary to avoid boundary issues
        m = np.linspace(0,1,self.nParticles + 2)
        # Calculate r from inverse CDF
        r = cdf_inv_r(m[1:-1]).astype(np.float32)
        # Assign output
        self.r = r
        
    def _generate_r_random(self):
        """
        Randomly generate radial positions according to the radial PDF
        
        PDF ~ 2*pi*R*sigma
        """
        cdf_inv_r = self._parent.sigma.cdf_inv
        np.random.seed(self._seed)
        m = np.random.rand(self.nParticles)
        r = cdf_inv_r(m).astype(np.float32)
        self.r = r
    
    def _generate_r(self):
        """
        Wrapper function to generate radial positions depending on chosen 
        method
        """
        
        print 'Generating r positions'
        method = self.method
        
        if method == 'glass':
            
            self._generate_r_glass()
            
        elif method == 'grid':
            
            self._generate_r_grid()
            
        elif method == 'random':
            
            self._generate_r_random()
            
        else:
            
            raise ValueError, 'Unrecognized method for generating positions {}'\
            .format(method)
            
    def _generate_z(self):
        """
        Generate z positions
        """
        
        print 'Generating z positions'
        if self.method == 'glass':
            # The inverse CDF over z as a function of r
            cdf_inv_z = self._parent.rho.cdf_inv
            # Glassy numbers between -1 and 1 already in self.z
            # Calculate z
            z = cdf_inv_z(np.abs(self.m_z), self.r)
            z = z * np.sign(self.m_z)
            # For memory reasons, delete self.m_z (it's not needed now)
            del self.m_z
            # Assign output
            self.xyz[:,2] = z

        else:
            # The inverse CDF over z as a function of r
            cdf_inv_z = self._parent.rho.cdf_inv
            # Random numbers between 0 and 1
            np.random.seed(self._seed)
            m = np.random.rand(self.nParticles)
            # Calculate z
            z = cdf_inv_z(m, self.r)
            # Randomly select sign of z
            z = z * np.random.choice(np.array([-1,1]), self.nParticles)
            # Assign output
            self.xyz[:,2] = z
        
    def _generate_theta(self):
        """
        Generate angular positions
        """
        
        nParticles = self.nParticles
        
        if self.method == 'glass':
            
            #already done in generate_r
            assert(len(self.theta)==nParticles)
            
        if self.method == 'grid':
            
            r = self.r
            
            dtheta = np.sqrt(2*np.pi*(1 - r[0:-1]/r[1:]))
            dtheta = strip_units(dtheta)
            theta = np.zeros(nParticles)
            
            for n in range(nParticles - 1):
                
                # NOTE: it's import to subtract (not add) dtheta.  The particles
                # will be moving counter-clockwise.  To prevent the particle
                # spirals from kinking, the particle spirals must go out
                # clockwise
                theta[n+1] = theta[n] - dtheta[n]
                
            self.theta = theta
            
        if self.method == 'random':
            
            np.random.seed(self._seed)
            theta = 2*np.pi*np.random.rand(nParticles)
            self.theta = theta
            
    def _cartesian_pos(self):
        """
        Generate x,y
        """
        
        r = self.r
        theta = self.theta
        self.xyz[:,0] = r*np.cos(theta)
        self.xyz[:,1] = r*np.sin(theta)

def cdf_space_height(IC, rmax=1.):
    """
    Calculates a reasonable zmax for a uniformly distributed cylinder in 
    CDF space (used by the 'glass' method for generating positions).
    
    The cylinder is of radius rmax assumed to extend to +/- zmax
    """
    rho = IC.rho.rho_binned
    z = IC.rho.z_bins
    rbins = IC.rho.r_bins
    
    # Integral of rho^2 along z
    rho2int = np.trapz(rho**2, z, axis=0)
    # Integral of Prob(z)^2 along z
    pz2int = rho2int/(0.5*IC.sigma(rbins))**2
    pz2int_spl = interp1dunits(rbins, pz2int, bounds_error=False,
                               fill_value='nearest')
    # integrate along r
    r = IC.sigma.r_bins
    pz2int_r = pz2int_spl(r)
    sigma_r = IC.sigma(r)
    cdf_r = cumtrapz(r*sigma_r, r, initial=0.)
    cdf_r /= cdf_r[-1]
    return (0.5 *  rmax) /np.trapz(pz2int_r * np.sqrt(cdf_r), r)


def mean_height_ratio(IC):
    """
    Calculate (disk height/R) averaged over the disk, weighted by mass
    """
    z = IC.rho.z_bins
    rho = IC.rho.rho_binned
    r = IC.rho.r_bins
    sigma = IC.sigma(r)
    H = np.sqrt(np.trapz(rho*z**2, z, axis=0)/(0.5*sigma))
    h_mean = np.trapz(sigma*H, r)/np.trapz(r*sigma, r)
    return float(h_mean.in_units('1'))