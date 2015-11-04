# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:28:37 2015

@author: ibackus
"""
# external modules
import cPickle as pickle
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
import numpy as np
import pynbody as pb
SimArray = pb.array.SimArray

# diskpy modules
from diskpy.pdmath import meshinterp, resolvedbins
from diskpy.disk import rho0_est, h_est
from diskpy.utils import strip_units
from vertical_solver import vertical_solver

# Constants
G = SimArray(1.0, 'G')
kB = SimArray(1.0, 'k')

class rhosolver():
    """
    Defines the rho class that allows the solving of vertical hydrostatic 
    equilibrium over the disk, and generates callable methods for estimating
    density and the normalized inverse CDF over the disk.
    
    Examples
    --------
    
    Initialize rho, solve vertical equilibrium
    
    >>> IC.maker.sigma_gen()
    >>> rho = rhosolver(IC)
    >>> rho.solve(maxiter=100)
    
    Rho is now callable:
    
    >>> rho(z, r)
    >>> rho.cdf_inv(r, m)
    
    Save rho
    
    >>> rho_dict = {'z': rho.z_bins, 'r': rho.r_bins, 'rho': rho.rho_binned}
    >>> pickle.dump(rho_dict, open('rhofile.p','w'))
    
    Load rho
    
    >>> rho.load('rhofile.p')
    """
    def __init__(self, IC):
        
        self._parent = IC
        self.solver_options = {}
        self.rho = None
        
    def __call__(self, z, r):
        """
        Call method for rhosolver objects.  Returns rho estimated at (z,r) if
        equilibrium as already been solved
        """
        if self.rho is not None:
            
            return self.rho(z, r)
    
    def solve(self, **kwargs):
        """
        Solves the hydrostatic vertical equilibrium for the disk to find the
        density.  Also calculates the normalized inverse CDF
        
        Parameters
        ----------
        
        kwargs :
            (optional) key word arguments to pass to the root finder
            [scipy.optimize.newton_krylov]
        """
        
        options = self.solver_options.copy()
        options.update(kwargs)
        z, r, rho = calc_rho(self._parent, **options)
        self.r_bins = r
        self.z_bins = z
        self.rho_binned = rho
        self.rho = rhointerp(z, r, rho)
        self.cdf_inv = cdf_inv(z, r, rho)
        
        return
        
    def load(self, f):
        """
        Load/initialize from a dictionary or a file
        
        Parameters
        ----------
        
        f : dict or str
            Either a dictionary containing everything required to load rho
            or a filename pointing to a pickled dictionary
        """
        
        if isinstance(f, str):
        
            f = pickle.load(open(f,'r'))
    
        if 'z' in f:
            # The dictionary should contain r, z, and rho (binned)
            z = f['z']
            r = f['r']
            rho = f['rho']
            
            # Older versions of rho only had one set of z-bins
            z = _updatez(z, f)
            
            # Initialize
            self.r_bins = r
            self.z_bins = z
            self.rho_binned = rho
            self.rho = rhointerp(z, r, rho)
            self.cdf_inv = cdf_inv(z, r, rho)
            
        else:
            # Not sure what's going on!
            raise ValueError, 'Could not load rho'
            
        return
        
def setup_r_bins(IC, r=None):
        
    if IC.settings.rho_calc.nr is not None:
        
        nr = IC.settings.rho_calc.nr
        rmax = IC.sigma.r_bins[[-1]]
        rbins = np.linspace(0, rmax, nr)
        return rbins
    
    if r is None:
        # Setup the initial r bins
        rmax = IC.sigma.r_bins[[-1]]
        nr = len(IC.sigma.r_bins) * 10
        #r = np.linspace(0, rmax, nr)
        # dflemin3 Nov 4, 2015: made units more explicit
        # via SimArrays
        r_units = IC.sigma.r_bins.units
        r = SimArray(np.linspace(0, rmax, nr),r_units)        
        
    bin_error_tol = IC.settings.rho_calc.r_bin_tol
    minbins = IC.settings.rho_calc.min_r_bins
        
    # Estimate the disk height
    M = IC.settings.physical.M
    m = IC.settings.physical.m
    T = IC.T(r)
    h = h_est(r, M, m, T, gamma=1)
    
    # Estimate midplane density
    sigma = IC.sigma(r)
    rho0 = rho0_est(h, sigma)
    
    # Estimate a reasonable function tolerance for the bins
    # This is done by taking a weighted mean of midplane density: weighted
    # by the number of particles at each radius (the PDF)
    w = abs(IC.sigma.pdf(r))
    w /= w.sum()
    w = strip_units(w)
    # also weight the midplane density
    rho0 = w*rho0
    # Now do the mean
    rho0mean = (rho0*w).sum()
    ftol = bin_error_tol * rho0mean
    
    # Estimate reasonable bins.  This is done by removing as many bins as
    # possible to still allow the midplane density to be well resolved
    rbins = resolvedbins(r, rho0, minbins=minbins, ftol=ftol)
    rbins = rbins.copy()
    
    print '{} radial bins used for density calculation'.format(len(rbins))
    
    return rbins
        

def _updatez(z, rhoDict):
    """
    Older versions of rho only had one set of z-bins
    ie, z.shape = (nz, )
    Newer versions should have z.shape = (nz, nr)
    """
    r = rhoDict['r']
    rho = rhoDict['rho']
    
    if z.shape != rho.shape:

        if np.ndim(z) == 1:
            z = z[:,None]
        if (np.ndim(z) != 2) or (z.shape[-1] != 1):
            
            raise ValueError, 'Could not understand z input'
            
        else:
            
            nr = len(r)
            z = np.dot(z, np.ones([1, nr]))
            
    return z
     
            
def loadrho(IC, f):
    """
    Loads a rhosolver object from f (a file or dictionary) and saves it to
    IC.rho
    
    Parameters
    ----------
    
    f : dict or str
        Either a dictionary containing everything required to load rho
        or a filename pointing to a pickled dictionary
        
    Returns
    -------
    
    rho : rhosolver
        An intialied rho solver object
    """
    
    rho = rhosolver(IC)
    rho.load(f)
    
    IC.rho = rho
    return
        
        
def rhointerp(z, r, rho):
    """
    Generates a callable interpolation of rho on the z, r points
    
    Parameters
    ----------
    
    z : 2D SimArray or array
        z[i,j] is the ith z value at r[j]
    r : 1D SimArray or array
        Radial positions
    rho : SimArray
        density at points z[i,j], r[j]
        
    Returns
    -------
    
    rhospline : function
        An interpolation function for estimating rho as a function of z, r
    """
    
    f = meshinterp(r, z.T, rho.T)
    
    def rhospline(Z, R):
        """
        Returns rho evaluated at Z, R.  Z and R must be 1D and the same length
        
        Parameters
        ----------
        
        Z, R : SimArray, array, or float
            Z, R positions to estimate the density at.  Must be same length
            
        Returns
        -------
        
        rho : SimArray
            Density evaluated at the Z,R positions
        """
        return f(R, Z)
        
    return rhospline
    
    
def calc_rho(IC, r=None, **kwargs):
    """
    Calculates the density for the initial conditions object IC by numerically
    solving hydrostatic equilibrium (see vertical_solver)
    
    Parameters
    ----------
    
    IC : IC object
        Initial conditions object
    r : SimArray
        (optional) intial bin radii: not all will be used
    **kwargs : keyword arguments to pass to the root finder used
        (scipy.optimize.newton_krylov)
        
    Returns
    -------
    
    R : 1D SimArray
        Radial bins the density is calculated at
    z : 2D SimArray
        z points the density is calculated at.  z[i,j] is the ith z position
        at R[j]
    rho : 2D SimArray
        Density as a function of (z,R).  rho[i,j] is calculated at 
        (z[i,j], R[j])
    """
    R = setup_r_bins(IC, r)
    nr = len(R)
    nz = IC.settings.rho_calc.nz
    # Initialize z and rho.  these get transposed at the end!
    z = np.zeros([nr, nz])
    rho = np.zeros([nr, nz])
    
    nPrint = 10
    iPrint = nr/nPrint
    for i in range(nr):
        
        if (i%iPrint) == 0:
            
            print 'Calculating rho:\t{0:.1f} %'.format(100.*i/nr)
        
        if R[[i]] > 0:
            
            rtf = vertical_solver(IC, R[[i]])
            rtf.fit(**kwargs)
            rho[i] = rtf.results['rho']
            z[i] = rtf.results['z']
    
    z = SimArray(z, rtf.results['z'].units)
    z[0] = z[1]
    rho = SimArray(rho, rtf.results['rho'].units)
    
    # Transpose rho and z to return rho as a function of z, R
    z = z.T
    rho = rho.T
    
    return z, R, rho
    
def cdf_inv(z, r, rho):
    """
    Calculates the inverse CDF as a function of r over the whole disk
    
    Parameters
    ----------
    
    z : SimArray or array
        2D z-positions that rho is calculated at.  z[i,j] is the ith z bin
        at r[j]
    r : SimArray or array
        Radial bins (1D array) the z and rho are calculated at
    rho : SimArray or array
        2D array of density values.  rho[i,j] is rho at z[i,j], r[j]
        
    Returns
    -------
    
    f : function
        Inverse CDF.  f(m, r) = z returns the z value for radius and 0<m<1.
        r and m are 1-D arrays of the same length, or numbers.
    """
    
    nz, nr = z.shape
    f = np.zeros(z.shape)
    zout = 0.*z
    for i in range(nr):
        
        f[:,i], zout[:,i] = cdf_inv_z(z[:,i], rho[:,i])
        
    cdf_inv_spl = meshinterp(r, f.T, zout.T)
    
    def fspl(m, R):
        """
        Normalized inverse CDF at R.  Calcuates z as a function of m
        
        Parameters
        ----------
        
        m : 1D array or float
            Number(s) between 0 and 1
        R : 1D SimArray, array, or float
            Radius at which to calculate the CDF inverse along rho
        
        Returns
        -------
        
        z : SimArray
            z positions 
        """
        return cdf_inv_spl(R, m)
        
    return fspl

def cdf_inv_z(z,rho):
    """
    Calculates the inverse of the cumulative distribution function for
    probability as a function of z for a given r (ie gives you z as a function
    of the CDF)
    
    Parameters
    ----------
    
    z : array or SimArray
        Positions to calculate over.  1D array
            
    rho: array or SimArray 
        Density as a function of z.  Treated as an un-normalized PDF. 1D array
    
    Returns
    -------
        
    f : array
        Normalized CDF
    z : array or SimArray
        z as a functin of the normalized CDF.  Monotonically increasing
        
    Notes
    -----
    
    To ensure z, f are montonically increasing, some values are dropped.
    The remaining values of f are padded with 2, and the values of z are
    padded with z.max()
    """
        
    # Calculate the CDF from prob
    nz = len(z)
    f = np.zeros(nz)
    f[1:] = cumtrapz(rho,z)
    if f.max() <= 0.0:
        # The density (rho) is zero here for all z or neg or something.
        # Make all particles go to z = 0.0
        f = np.linspace(0, 1, nz)
        z = 0.*z
        return f, z
    f /= f.max()
    # Calculate the inverse CDF.
    # Assume CDF is approximately monotonic and sort to force it to be
    ind = f.argsort()
    f = f[ind]
    z = z[ind]
    # Drop values where CDF is constant (ie, prob = 0)
    mask = np.ones(nz,dtype='bool')
    for n in range(1,nz):
        if f[n] == f[n-1]:
            mask[n] = False
            
    nVals = mask.sum()
    f[0:nVals] = f[mask]
    z[0:nVals] = z[mask]
    
    if nVals < nz:
        # Pad the remainder of values
        z[nVals:] = z.max()
        f[nVals:] = 2.
    
    return f, z
    
#def setup_r_bins(IC):
#    """
#    Set up radial bins for calculating density (rho) as a function of (r, z)
#    
#    Parameters
#    ----------
#    
#    IC : ICobject
#        An initial conditions object
#    
#    Returns
#    -------
#    
#    R : SimArray
#        Radial bins
#    """
#    
#    nr = IC.settings.rho_calc.nr
#    
#    if nr is None:
#        
#        R = IC.sigma.r_bins
#        
#    else:
#        
#        rmax = IC.sigma.r_bins[[-1]]
#        R = np.linspace(0, rmax, nr)
#        
#    return R
