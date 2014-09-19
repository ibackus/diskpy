# -*- coding: utf-8 -*-
"""
Calculates rho(z) to maintain hydrostatic equilibrium in a thin disc.  
Assumes uniform temperature in the disc, and an infinite disc where 
rho can be treated (at least locally) as only a function of z.

Created on Mon Jan 20 12:30:06 2014

@author: ibackus
"""
# ICgen packages
import isaac

# External packages
import numpy as np
import scipy
import scipy.integrate as nInt
import scipy.optimize as opt
from scipy.interpolate import interp1d
import pynbody
from pynbody.array import SimArray
from warnings import warn
import sys

def rho_z(sigma, T, r, settings):
    """ 
    rho,z = rho_z(...)
    
    Calculates rho(z) to maintain hydrostatic equilibrium in a thin disc.  
    Assumes uniform temperature in the disc, and an infinite disc where 
    rho can be treated (locally) as only a function of z.
    
    Only calculates for z>=0, since the density is assumed to be symmetric
    about z=0
    
    The initial guess for rho (a gaussian) only really seems to work for
    Mstar >> Mdisc.  Otherwise the solution can diverge violently.
    
    * NUMERICAL CALCULATION OF RHO(Z) *
    The calculation proceeds using several steps.
    1) Make an initial guess for I, the integral of rho from z to inf.  This
        is an error function
    2) Modify length scale of the initial guess to minimize the residual
        for the differential equation governing I.  Use this as the new
        initial guess.
    3) Find the root I(z) for the differential equation governing I, with
        the boundary condition that I(0) = sigma/2
    4) Set rho = -dI/dz
    5) Find the root rho(z) for the diff. eq. governing rho.
    6) In order to satisfy the BC on I, scale rho so that:
        Integral(rho) = I(0)
    7) Repeat (5) and (6) until rho is rescaled by a factor closer to unity
        than rho_tol
        
        Steps 5-7 are done because the solution for I does not seem to
        satisfy the diff. eq. for rho very well.  But doing it this way 
        allows rho to satisfy the surface density profile
    
    * Arguments *
    
    sigma - The surface density at r
    __stdout__
    T - the temperature at r
    
    r - The radius at which rho is being calculated.  Should have units
    
    settings - ICobj settings (ie, ICobj.settings)
        
    * Output *
    Returns a 1D SimArray (see pynbody) of rho(z) and a 1D SimArray of z,
    with the same units as ICobj.settings.rho_calc.zmax
    
    """
    # Parse settings
    rho_tol = settings.rho_calc.rho_tol
    nz = settings.rho_calc.nz
    zmax = settings.rho_calc.zmax
    
    m = settings.physical.m
    M = settings.physical.M
    
    # Physical constants
    kB = SimArray(1.0,'k')
    G = SimArray(1.0,'G')
    
    # Set up default units
    mass_unit = M.units
    length_unit = zmax.units
    r = (r.in_units(length_unit)).copy()
    
    # Initial conditions/physical parameters
    rho_int = 0.5*sigma.in_units(mass_unit/length_unit**2)   # Integral of rho from 0 to inf
    a = (G*M*m/(kB*T)).in_units(length_unit)
    b = (2*np.pi*G*m/(kB*T)).in_units(length_unit/mass_unit)
    z0guess = np.sqrt(2*r*r*r/a).in_units(length_unit)# Est. scale height of disk
    z0_dummy = (2/(b*sigma)).in_units(length_unit)

    z = np.linspace(0.0,zmax,nz)
    dz = z[[1]]-z[[0]]
    
    # Echo parameters used
    print '***********************************************'
    print '* Calculating rho(z)'
    print '***********************************************'
    print 'sigma          = {0} {1}'.format(sigma,sigma.units)
    print 'zmax           = {0} {1}'.format(zmax,zmax.units)
    print 'r              = {0} {1}'.format(r,r.units)
    print 'molecular mass = {0} {1}'.format(m,m.units)
    print 'Star mass      = {0} {1}'.format(M,M.units)
    print 'Temperature    = {0} {1}'.format(T,T.units)
    print ''
    print 'rho_tol        = {0}'.format(rho_tol)
    print 'nz             = {0}'.format(nz)
    print '***********************************************'
    print 'a              = {0} {1}'.format(a,a.units)
    print 'b              = {0} {1}'.format(b,b.units)
    print 'z0guess        = {0} {1}'.format(z0guess,z0guess.units)
    print '***********************************************'
    
    print 'z0 (from sech^2) = {0} {1}'.format(z0_dummy,z0_dummy.units)
    # --------------------------------------------------------
    # STRIP THE UNITS FROM EVERYTHING!!!
    # This has to be done because many of the scipy/numpy functions used cannot
    # handle pynbody units.  Before returning z, rho, or anything else, the
    # Units must be re-introduced
    # --------------------------------------------------------
    
    rho_int, a, b, z0guess, z0_dummy, z, dz, r, T, sigma \
    = isaac.strip_units([rho_int, a, b, z0guess, z0_dummy, z, dz, r, T, sigma])
    
    # --------------------------------------------------------
    # Check sigma and T
    # --------------------------------------------------------
    if sigma < 1e-100:
        
        warn('Sigma too small.  setting rho = 0')
        rho0 = np.zeros(len(z))
        # Set up units
        rho0 = isaac.set_units(rho0, mass_unit/length_unit**3)
        z = isaac.set_units(z, length_unit)
        
        return rho0, z
        
    if T > 1e100:
        
        warn('Temperature too large.  Setting rho = 0')
        rho0 = np.zeros(len(z))
        # Set up units
        rho0 = isaac.set_units(rho0, mass_unit/length_unit**3)
        z = isaac.set_units(z, length_unit)
        
        return rho0, z  
        
    # -------------------------------------------------------------------
    # FUNCTION DEFINITIONS
    # -------------------------------------------------------------------
    def dI_dz(I_in):
        """
        Finite difference approximation of dI/dz, assuming I is odd around I(0)
        """
        I = I_in.copy()
        dI = np.zeros(len(I))
        # Fourth order center differencing
        dI[0] = (-I[2] + 8*I[1] - 7*I[0])/(6*dz)
        dI[1] = (-I[3] + 8*I[2] - 6*I[0] - I[1])/(12*dz)
        dI[2:-2] = (-I[4:] + 8*I[3:-1] -8*I[1:-3] + I[0:-4])/(12*dz)
        # Second order backward differencing for right edge
        dI[-2:] = (3*I[-2:] -4*I[-3:-1] + I[-4:-2])/(2*dz)
        
        return dI
        
    def d2I_dz2(I_in):
        
        # Finite difference for d2I/dz2 assuming it is 0 at the origin
        I = I_in.copy()
        d2I = np.zeros(len(I))
        # Boundary condition
        d2I[0] = 0
        # Centered 4th order finite difference
        d2I[1] = (-I[3] + 16*I[2] - 30*I[1] + 16*I[0] -(2*I[0] - I[1]))/(12*dz**2)
        d2I[2:-2] = (-I[4:] + 16*I[3:-1] - 30*I[2:-2] + 16*I[1:-3] - I[0:-4])/(12*(dz**2))
        # second order backward difference for right edge
        d2I[-2:] = (-2*I[-2:] + 5*I[-3:-1] -4*I[-4:-2] + I[-5:-3])/dz**2
        
        return d2I
        
    def Ires(I_in):
        """
        Calculate the residual for the differential equation governing I,
        the integral of rho from z to "infinity."
        """
        # DEFINE INITIAL CONDITION:
        I = I_in.copy()
        I[0] = rho_int
        #I[-1] = 0.0
        weight = 1.0

        res = d2I_dz2(I) + dI_dz(I)*(a*z/((z**2 + r**2)**(1.5)) + 2*b*(I[0] - I))
        
        return weight*res
        
    def drho_dz(rho_in):
        """
        Fourth order, centered finite difference for d(rho)/dz, assumes that
        rho is an even function.  The right-hand boundary is done using
        backward differencing
        """
        rho = rho_in.copy()
        drho = np.zeros(len(rho))
        drho[0] = 0.0   # defined by boundary condition, rho[0] = max(rho)
        drho[1] = (-rho[3] + 8*rho[2] - 8*rho[0] + rho[1])/(12*dz)
        drho[2:-2] = (-rho[4:] + 8*rho[3:-1] - 8*rho[1:-3] + rho[0:-4])/(12*dz)
        drho[-2:] = (3*rho[-2:] - 4*rho[-3:-1] + rho[-4:-2])/(2*dz)
        
        return drho
        
    def residual(rho_in):
        """
        Estimate d(rho)/dz
        """
        rho = rho_in.copy()
        # Estimate integral of rho
        I = np.zeros(len(rho))
        I[1:] = nInt.cumtrapz(rho,z)
        # Estimate residual 
        res = drho_dz(rho) + a*rho*z/((z**2 + r**2)**(1.5)) + 2*b*rho*I
        
        return res
        
    def erf_res(scale_size):
        
        testfct = rho_int*(1 - scipy.special.erf(z/scale_size))
        
        return abs(Ires(testfct)).sum()
        

    pass
    # -------------------------------------------------------------------
    # FIND RHO
    # -------------------------------------------------------------------
    
    maxiter = 40
    
    # Estimate the scale length of the error function
    z0 = opt.fminbound(erf_res,z0guess/100.0,5.0*z0guess)
    print 'Length scale guess: {0} {1}'.format(z0guess, length_unit)
    print 'Final length scale: {0} {1}'.format(z0, length_unit)
    
    # Begin by finding I, the integral of rho (from z to inf)
    # Assuming rho is gaussian, I is an error function
    guess = rho_int*(1 - scipy.special.erf(z/z0))
    
    # Find the root of the differential equation for I
    f_tol = rho_int * 6e-6
    try:
        
        Isol = opt.newton_krylov(Ires,guess,maxiter=maxiter,f_tol=f_tol)
        
    except :
        # Assume it didn't converge because f_tol was too strict
        # Read exception
        xepshun = sys.exc_info()
        # Extract rho from the exception
        Isol = xepshun[1][0]
        
    # rho is the negative derivative
    rho0 = -dI_dz(Isol)
    
    # Now apply the diff eq on rho
    for n in range(maxiter):
        
        print 'Iteration {0}'.format(n+1)
        f_tol = rho0.max() * 6e-6
        try:
            
            rho0 = opt.newton_krylov(residual,rho0,maxiter=maxiter, f_tol=f_tol)
            
        except:
            # Assume it didn't converge because f_tol was too strict
            # Read exception
            xepshun = sys.exc_info()
            # Extract rho from the exception
            rho0 = xepshun[1][0]
            
        rho_scale = rho_int/nInt.cumtrapz(rho0,z)[-1]
        print 'Scaling rho by {0}'.format(rho_scale)
        rho0 = rho0*rho_scale
        
        if abs(1-rho_scale) < rho_tol - 1:
            
            break
        
    if n >= maxiter:
        
        print 'Warning: solution to rho did not converge for r = {0}'.format(r)
    
    # Re-introduce units
    rho0 = isaac.set_units(rho0, mass_unit/length_unit**3)
    z = isaac.set_units(z, length_unit)
    
    return SimArray(rho0,'Msol au**-3'), SimArray(z,'au')
    
def cdfinv_z(z,rho):
    """
    Calculates the inverse of the cumulative distribution function for
    probability as a function of z for a given r
    
    *** Arguments ***
    
    * z *   z positions to calculate over.  1D array
            
    * rho *     Density as a function of z.  Treated as an un-normalized
    probability.  1D array
    
    IF Z doesn't have units, units of 'au' are assumed
    
    *** Returns ***
        
    Returns the inverse normalized CDF as 1D spline interpolation
    """
    # Check for units
    if pynbody.units.has_units(z):
        
        zunit = z.units
        
    else:
        
        zunit = pynbody.units.au
        
    # Calculate the CDF from prob
    nz = len(z)
    f = np.zeros(nz)
    f[1:] = nInt.cumtrapz(rho,z)
    if f.max() <= 0.0:
        # The density (rho) is zero here for all z or neg or something.
        # Make all particles go to z = 0.0
        def finv(m_in):
            return m_in*0.0
        return finv
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
    f = f[mask]
    z = z[mask]
    finv_spline = interp1d(f,z,kind='linear')
    
    def finv(m):
        
        return SimArray(finv_spline(m), zunit)
        
    return finv

#def cdfinv_z(z,rho):
#    """
#    Calculates the inverse of the cumulative distribution function for
#    probability as a function of z for a given r
#    
#    *** Arguments ***
#    
#    * z *   z positions to calculate over.  1D array
#            
#    * rho *     Density as a function of z.  Treated as an un-normalized
#    probability.  1D array
#    
#    IF Z doesn't have units, units of 'au' are assumed
#    
#    *** Returns ***
#        
#    Returns the inverse normalized CDF as 1D spline interpolation
#    """
#    # Check for units
#    if pynbody.units.has_units(z):
#        
#        zunit = z.units
#        
#    else:
#        
#        zunit = pynbody.units.au
#        
#    # Calculate the CDF from prob
#    nz = len(z)
#    f = np.zeros(nz)
#    f[1:] = nInt.cumtrapz(rho,z)
#    if f.max() <= 0.0:
#        # The density (rho) is zero here for all z or neg or something.
#        # Make all particles go to z = 0.0
#        def finv(m_in):
#            return m_in*0.0
#        return finv
#    f /= f.max()
#    # Calculate the inverse CDF.
#    # Assume CDF is approximately monotonic and sort to force it to be
#    ind = f.argsort()
#    f = f[ind]
#    z = z[ind]
#    # Drop values where CDF is constant (ie, prob = 0)
#    mask = np.ones(nz,dtype='bool')
#    for n in range(1,nz):
#        if f[n] == f[n-1]:
#            mask[n] = False
#    f = f[mask]
#    z = z[mask]
#    finv_spline = interp1d(f,z,kind='linear')
#    
#    def finv(m):
#        
#        return SimArray(finv_spline(m), zunit)
#        
#    return finv
