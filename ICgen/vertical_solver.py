# -*- coding: utf-8 -*-
"""
Defines the vertical_solver class for solving vertical hydrostatic equilibrium
at a fixed radius.


Created on Mon Aug  3 16:01:44 2015

@author: ibackus
"""

# Internal imports
from diskpy.utils import strip_units

# External imports
import numpy as np
import pynbody as pb
SimArray = pb.array.SimArray
from scipy.integrate import cumtrapz, trapz
from scipy.optimize.nonlin import NoConvergence
from scipy.optimize import newton_krylov, fmin

# Constants
G = SimArray(1.0, 'G')
kB = SimArray(1.0, 'k')

class vertical_solver():
    """
    Defines the class which solves vertical hydrostatic equilibrium for a thin
    Keplerian disk orbiting a central star
    
    Examples
    --------
    
    >>> IC = ICgen.load('IC.p')
    >>> R = SimArray(0.5, 'au')
    >>> solver = vertical_solver(IC, R)
    >>> solver.fit(maxiter=100)
    >>> z = solver.results['z']
    >>> rho = solver.results['rho']
    
    This example loads up an already created initial conditions object then
    solves for hydrostatic equilibrium at R=0.5 AU

    """
    
    def __init__(self, IC, R, rescale_rho=True):
        
        # initialize parameters
        z, r, c, zscale, rhoscale = setup(IC, R)
        UNITS = {'length': zscale, 'rho': rhoscale}
        self.UNITS = UNITS
        
        # Initialize z-position stuff
        nz = len(z)
        dz = z[1] - z[0]
        POS = {'z': z, 'nz': nz, 'Z': zscale*z, 'R': R, 'dz': dz}
        self.POS = POS
        
        # Constants used for calculations
        A = r**3 * z * (z**2 + r**2)**-1.5
        #B = 2*c*z        
        CONSTANTS = {'c':c, 'r': r, 'A': A}
        self.CONSTANTS = CONSTANTS
        
        # Boolean flag (rescale rho as part of the residual)
        self.rescale_rho = rescale_rho
        
        # Perform initial guess
        h = self._fit_h()
        self.CONSTANTS['h'] = float(h)        
        self.rho = rho0(z, h)
        self.rho0 = self.rho.copy()
        
        # Initialize empty results
        self.results = {}
        
        
    def _h_residual(self, h):
        """
        Estimates sum of the squares of the residuals for a density of the 
        form exp(-z^2/2h^2) (see rho0)
        """
        z = self.POS['z']
        rho = rho0(z, h)
        res = self.residual(rho)
        return (res**2).sum()
        
    def _fit_h(self, h0=1):
        """
        Fit a functional form to minimize the sum of the squares of the 
        residual for density of form exp(-z^2/2h^2)
        """
        sol = fmin(self._h_residual, h0, disp=False)
        
        return sol
        
    def _drho(self, rho):
        """
        Estimates the derivative of rho with respect to logz
        """
        
        dz = self.POS['dz']
        grad = np.gradient(rho, dz)
        grad[0] = 0
        return grad
        
    def _rhoint(self, rho):
        """
        Estimates (cumulative) integral of rho over z
        (equivalently, the integral of rho*z over logz)
        """
        
        dz = self.POS['dz']
        
        return cumtrapz(rho, dx = dz, initial = 0)
        
    def _rescale(self, rho):
        """
        rescales rho such that the integral of rho from 0 to inf = 1/2
        """
        # Ignore all negative numbers
        pos = (rho >= 0)
        # Rho should integrate to 1 over all space, or 1/2 from 0 to inf
        # (in these units)
        dz = self.POS['dz']
        scale = trapz(pos*rho, dx = dz)/0.5
        rho /= scale
        return rho
        
    def _setup_results(self):
        """
        Sets up the results to have units, etc...
        """
        rho = self.rho
        z = self.POS['z']
        dz = self.POS['dz']
        # estimate scale height
        h2 = trapz(rho * z**2, dx=dz)/trapz(rho, dx=dz)
        h = np.sqrt(h2)
        
        # Handle units
        rho2 = self.UNITS['rho'] * rho
        Z = self.UNITS['length'] * z
        H = self.UNITS['length'] * h
        
        self.results = {'z': Z, 'rho': rho2, 'h': H}
        
    def residual(self, rho):
        """
        Estimates the residual of the differential equation for rho for all
        z
        """
        
        A = self.CONSTANTS['A']
        c = self.CONSTANTS['c']
        
        if self.rescale_rho:
            
            rho = self._rescale(rho)
        
        drho = self._drho(rho)
        res = drho + A*rho + 2*c*rho*self._rhoint(rho)
        # penalize positive gradients
        drhonorm = drho/abs(drho).max()
        penalty = np.ones(len(rho))
        mask = drhonorm > 0
        penalty[mask] = 1 + 10*drhonorm[mask]
        
        # weight the residual
        weight = 1
            
        return weight * res * penalty
        
    def fit(self, **kwargs):
        """
        Attempts to find the root of the residual of the differential equation 
        for rho at all z, given the initial guess self.rho0
        
        Uses a Newton's method with the Krylov approximation for the Jacobian
        (see scipy.optimize.newton_krylov)
        
        
        Parameters
        ----------
        
        **kwargs : keyword arguments
            Additional arguments passed to the root finder
            [see scipy.optimize.newton_krylov]
            
        Returns
        -------
        
        rho : array
            Dimensionless density, saves it to self.rho
        """
        # Set up the default options
        options = {'maxiter': 150, 'verbose': False}
        for k, v in kwargs.iteritems():
            
            options[k] = v
        
        # initial guess
        rho0 = self.rho0
        
        # Find root
        try:
            
            rho = newton_krylov(self.residual, rho0, **options)
            
        except NoConvergence as err:
            # Assume it didn't converge because f_tol was too strict
            # Read exception
            print 'Didnt converge'
            rho = err.args[0]
        
        rho[rho < 0] = 0
        # Normalize rho
        if self.rescale_rho:
            
            rho = self._rescale(rho)
            
        self.rho = rho
        self._setup_results()        
        
    def fitrobust(self, nh=10, scansize=5, **kwargs):
        """
        Attempts to find the root of the residual of the differential equation 
        for rho at all z, given the initial guess generated by different values
        of h.
        
        This routine is slower than fit but can avoid local minima.  If fit
        doesn't seem to be converging to the right values, consider using this.
        
        h is a dimensionless disk height, which is near 1.  This algorithm
        scans over values of h to create initial guesses rho0.  fit() is then
        called, and the results are stored.  The fit value of rho which is
        globally closest to 0 is then used.
        
        Uses a Newton's method with the Krylov approximation for the Jacobian
        (see scipy.optimize.newton_krylov)
        
        
        Parameters
        ----------
        
        nh : int
            (optional) Number of h values to scan over
        scansize : float
            (optional) Sets the range of h values to scan over.  Scans from 
            h0/scansize to h0*scansize
        **kwargs : keyword arguments
            Additional arguments passed to the root finder
            [see scipy.optimize.newton_krylov]
            
        Returns
        -------
        
        rho : array
            Dimensionless density, saves it to self.rho
            
        """
        nz = self.POS['nz']
        z = self.POS['z']
        rhos = np.zeros([nh, nz])
        h = np.zeros(nh)
        residual = np.zeros(nh)
        h0 = self.CONSTANTS['h']
        h[0:-1] = np.linspace(h0/scansize, h0*scansize, nh-1)
        h[-1] = h0
        
        for i in range(nh):
            
            self.rho0 = rho0(z, h[i])
            self.fit(**kwargs)
            rhos[i] = self.rho.copy()
            residual[i] = abs(self.residual(rhos[i])).max()
            
        # Find best solution
        best = residual.argmin()
        rho = rhos[best]
        # Normalize residuals by rhomax
        residual /= rho.max()
        # measure height of the solution
        hmeas = np.sqrt((rho*z**2).sum()/rho.sum())
        print 'h(estimate): {}\nh(best): {}\nh(measured): {}'\
        .format(h[-1], h[best], hmeas)
        print 'residuals normalized by rho:'
        print '\tbest: {}'.format(residual[best])
        print '\test: {}'.format(residual[-1])
        self.rho = rho
        self._setup_results()
        
def dI(I, dz=1):
    """
    Derivative of I with respect to z
    """    
    return np.gradient(I, dz)

def setup(IC, R):
    """
    Initialize values for solving vertical hydrostatic equlibrium.
    
    Parameters
    ----------
    
    IC : IC object
        Initial conditions object
    R : SimArray
        Radius at which to solve hydrostatic equilibrium
    
    Returns
    -------
    
    z : array
        Dimensionless z bins, in units of h
    r : float
        Radius in units of h
    c : float
        Constant in residual equation
    zscale : SimArray
        Amount to scale z (all lengths) by to get dimensions back.
        ie, Z(actual) = z*zscale
    rhoscale : SimArray
        Scale for dimensionless rho
    """
    
    # Load everything from IC
    sigma = IC.sigma(R)
    T = IC.T(R)
    M = IC.settings.physical.M
    m = IC.settings.physical.m
    zmax = IC.settings.rho_calc.zmax
    nz = IC.settings.rho_calc.nz
    
    # Set-up constants
    a = G*M*m/(kB*T)
    b = 2*np.pi*G*m/(kB*T)
    h = np.sqrt(R*R*R/a).in_units('au')
    c = b*h*sigma
    c.convert_units('1')
    zscale = h
    rhoscale = sigma/h
    
    # Set-up z bins
    if zmax is None:
        
        zmax = 6*h
        
    zmax = strip_units((zmax/h).in_units('1'))
    z = np.linspace(0, zmax, nz)
    
    # Make everything dimensionless
    r = strip_units((R/h).in_units('1'))
    c = strip_units(c)
    
    return z, r, c, zscale, rhoscale
    
def rho0(z, h=1):
    """
    Approximate functional form for rho.  All quantities are dimensionless
    
    To first order, h should equal 1.
    
    Parameters
    ----------
    
    z : Array
        z positions to estimate rho at
    h : float
        dimensionless scale height.  Default=1
    
    Returns
    -------
    
    rho : array
        Dimensionless density
    """
    
    return np.exp(-z**2/(2*h**2))/(h*np.sqrt(2*np.pi))