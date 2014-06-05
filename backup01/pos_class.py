# -*- coding: utf-8 -*-
"""
Defines a function to randomly generate particle positions according to 
the desired surface density profile (sigma vs r) and the vertical profile
(rho vs r,z).

Created on Mon Jan 27 18:48:04 2014

@author: ibackus
"""

import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import cPickle as pickle
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
import scipy.optimize as opt

# LOCAL IMPORTING
import calc_sigma
import calc_rho_zr
import isaac

class pos:
    
    def __init__(self, ICobj, method = None):
        
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
        print 'Generating {} particle positions using method: {}'.format(\
        self.nParticles, self.method)
        
        # Generate positions
        self._generate_r()
        self._generate_z()
        self._generate_theta()
        self._cartesian_pos()
        
    def __getstate__(self):
        """
        This is required to make the object pickle-able
        """
        
        # Define a dictionary containing everything needed.  Ignore self.parent
        state = self.__dict__.copy()
        state.pop('_parent', None)
        
        return state
        
    
    def _generate_r(self):
        """
        Generate radial positions
        """
        
        print 'Generating r positions'
        cdf_inv_r = self._parent.sigma.cdf_inv
        
        if self.method == 'grid':
            
            # Generate linearly increasing values of m, using 2 more than
            # necessary to avoid boundary issues
            m = np.linspace(0,1,self.nParticles + 2)
            # Calculate r from inverse CDF
            r = cdf_inv_r(m[1:-1])
            # Assign output
            self.r = r
            
        if self.method == 'random':
            
            m = np.random.rand(self.nParticles)
            r = cdf_inv_r(m)
            self.r = r
            
    def _generate_z(self):
        """
        Generate z positions
        """
        
        print 'Generating z positions'
        # The inverse CDF over z as a function of r
        cdf_inv_z = self._parent.rho.cdf_inv
        # Random numbers between 0 and 1
        m = np.random.rand(self.nParticles)
        # Calculate z
        z = cdf_inv_z(m, self.r)
        # Randomly select sign of z
        z = z * np.random.choice(np.array([-1,1]), self.nParticles)
        # Assign output
        self.z = z
        
    def _generate_theta(self):
        """
        Generate angular positions
        """
        
        nParticles = self.nParticles
        
        if self.method == 'grid':
            
            r = self.r
            
            dtheta = np.sqrt(2*np.pi*(1 - r[0:-1]/r[1:]))
            dtheta = isaac.strip_units(dtheta)
            theta = np.zeros(nParticles)
            
            for n in range(nParticles - 1):
                
                # NOTE: it's import to subtract (not add) dtheta.  The particles
                # will be moving counter-clockwise.  To prevent the particle
                # spirals from kinking, the particle spirals must go out
                # clockwise
                theta[n+1] = theta[n] - dtheta[n]
                
            self.theta = theta
            
        if self.method == 'random':
            
            theta = 2*np.pi*np.random.rand(nParticles)
            self.theta = theta
            
    def _cartesian_pos(self):
        """
        Generate x,y
        """
        
        r = self.r
        z = self.z
        theta = self.theta
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        
        xyz = np.zeros([self.nParticles, 3])
        xyz = isaac.match_units(xyz, r)[0]
        
        xyz[:,0] = x
        xyz[:,1] = y
        xyz[:,2] = isaac.match_units(z, r)[0]
        
        self.x = x
        self.y = y
        self.xyz = xyz

def make(rhoFileName,sigFileName,nParticles,rlim = None,zlim = None,\
savename = None):
    """
    Generate particle positions on a pseudo-grid according to the surface density
    profile contained in sigFileName and the vertical profile contained in
    rhoFileName.  PDFs are interpolated using a cubic spline, and are 
    assumed to be 0 outside the regions defined the the input files
    
    * Arguments *
    rhoFileName - file name for a file containing a pickled dictionary
        for rho calculated at points on a grid.  The dictionary should
        contain:
            dict['rho'] : 2D array, rho at all pairs of points (z,r)
            dict['z']   : a 1D array of z points
            dict['r']   : a 1D array of r points
            
    sigFileName - file name for a file containing a pickled dictionary
        for sigma (surface density) calculated at points r.  Dictionary
        should contain:
            dict['sigma']   : 1D array, rho at all points r
            dict['r']       : 1D array of r points
    
    nParticles - total number of particles to generate
    
    rlim - cylindrical radius bounds for generating particles.  if None,
        then the default rlims are the range of values contained in 
        rhoFileName.
        
    zlim - positive z-bounds for generating particles.  Particles will be 
        always generated symmetrically, within these bounds (the sign of
        z is randomly chosen).  If None, bounds are taken from rhoFileName
    
    savename - if None, no file is saved.  Otherwise, the output is pickled
        and saved to 'savename'
    
    * Output *
    Returns a dictionary containing the randomly generated positions, with
    the following keys:
        
        'x'     : x-positions of the particles
        'y'     : y-positions
        'z'     : z-positions
        'r'     : monotonically increasing R-positions
        'theta' : angular positions
    
    
    !!!Everything should be in units of Msol and au!!!
    pynbody simarrays are acceptable
    """
    # ------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------
    a = pickle.load(open(rhoFileName,'rb'))
    rho = a['rho']
    z = a['z']
    r = a['r']
    if rlim == None:
        rlim = [r.min(),r.max()]
    if zlim == None:
        zlim = [z.min(),z.max()]
    # ------------------------------------------------------------
    # Generate PDFs
    # ------------------------------------------------------------
    #print('Generating rho(z,r) spline interpolation')
    #rho_zr = RectBivariateSpline(z,r,rho)
    cdfinv_zr = calc_rho_zr.cdfinv_zr(rhoFileName)
    
    # Probability as a function of r (2*pi*r*sigma, up to a constant factor)
    # Pr is normalized, but we need the MAX value of it
    print 'Calculating probability vs r spline'    
    Pr = calc_sigma.prob(sigFileName)
    def PrNeg(x):
        return -1.0*Pr(x)
    print 'Finding maximum probability for r = {0} to {1}'.format(rlim[0],rlim[1])
    #print PrNeg(rlim[1])
    xmax = opt.fminbound(PrNeg,float(rlim[0]),float(rlim[1]))
    Prmax = Pr(xmax)
    print 'Max. prob = {0} \nat r = {1}'.format(Prmax,xmax)
    # ------------------------------------------------------------
    # Generate random R values according to PDF from sigma vs r
    # ------------------------------------------------------------
    print('Calculating R') 
    R = np.zeros(nParticles)
    cdfinv_r = calc_sigma.cdfinv_r(sigFileName,Pr)
    # Calculate R.  Note, we drop the first particle position because it is
    # at the origin at the last particle position to keep it from being at the
    # edge of our radial bins
    m = np.linspace(0,1,nParticles + 2)
    R = cdfinv_r(m[1:-1])
    # ------------------------------------------------------------
    # Generate random Z values for the R values found above,
    # according to rho(z,r)
    # ------------------------------------------------------------
    print('Calculating Z')
    Z = np.zeros(nParticles)
    cdfvals = np.random.uniform(size=nParticles)
    rbins = np.digitize(R,r)
    dr = r[1]-r[0]
    for i in range(nParticles):
        rbin = rbins[i]
        zlo = cdfinv_zr[rbin-1](cdfvals[i])
        zhi = cdfinv_zr[rbin](cdfvals[i])
        Z[i] = zlo + ((zhi-zlo)/dr)*(R[i] - r[rbin-1])
    # ------------------------------------------------------------
    # Generate angles
    # ------------------------------------------------------------
    dtheta = np.sqrt(2*np.pi*(1 - R[0:-1]/R[1:]))
    theta = np.zeros(nParticles)
    for n in range(0,nParticles-1):
        theta[n+1] = theta[n] + dtheta[n]
    # ------------------------------------------------------------
    # Format the output
    # ------------------------------------------------------------
    # Assign units
    if pynbody.units.has_units(r):
        rUnits = r.units
    else:
        rUnits = pynbody.units.au
    if pynbody.units.has_units(z):
        zUnits = z.units
    else:
        zUnits = pynbody.units.au
    if pynbody.units.has_units(Z):
        Z.convert_units(zUnits)
    else:
        Z = SimArray(Z,zUnits)
    if pynbody.units.has_units(R):
        R.convert_units(rUnits)
    else:
        R = SimArray(R,rUnits)
    Z = Z*np.sign(np.random.randn(nParticles))
    x = R*np.cos(theta)
    y = R*np.sin(theta)
    # Add random noise to positions
    var = np.zeros(nParticles)
    var[1:] = (x[1:]-x[0:-1])**2 + (y[1:]-y[0:-1])**2
    var[0] = var[1]
    mu = np.zeros(2)
    xydelta = np.zeros([nParticles,2])
    rand_scale = 0.01
    for n in range(nParticles):
        cov = var[n]*np.identity(2)
        xydelta[n,:] = np.random.multivariate_normal(mu,cov)
    xydelta *= rand_scale
    x += xydelta[:,0]
    y += xydelta[:,1]
    R = np.sqrt(x**2 + y**2)
    
    # Output
    outDict = {'x': x,'y': y, 'z': Z, 'r': R, 'theta': theta}
    if savename is not None:
        pickle.dump(outDict,open(savename,'wb'))
    return outDict
            
    
    