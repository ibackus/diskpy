# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 10:17:53 2014

@author: ibackus
"""

# External modules
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray

# diskpy modules
from diskpy.pdmath import smoothstep
from diskpy.utils import match_units


def make_profile(ICobj):
    """
    A wrapper for generating surface density profiles according to the IC object.
    
    Settings for the profile are defined in ICobj.settings.  Which profile gets
    used is defined by ICobj.settings.sigma.kind
    
    Currently available kinds are:
    
    viscous
    powerlaw
    MQWS
    
    **RETURNS**
    
    r : SimArray
        Radii at which sigma is calculated
    sigma : SimArray
        Surface density profile as a function of R
    """
    kind = ICobj.settings.sigma.kind
    
    if kind == 'powerlaw':
        
        r, sigma = powerlaw(ICobj.settings, ICobj.T)
        
    elif (kind == 'mqws') | (kind == 'MQWS'):
        
        r, sigma = MQWS(ICobj.settings, ICobj.T)
        
    elif (kind == 'viscous'):
        
        r, sigma = viscous(ICobj.settings)
        
    else:
        
        raise TypeError, 'Could not make profile for kind {0}'.format(kind)
        
    if hasattr(ICobj.settings.sigma, 'innercut'):
        
        sigma = _applycut(r, sigma, ICobj.settings.sigma.innercut, False)
        
    if hasattr(ICobj.settings.sigma, 'outercut'):
        
        sigma = _applycut(r, sigma, ICobj.settings.sigma.outercut, True)
    
    return r, sigma
    
def _applycut(r, sigma, rcut, outer=True):
    """
    Applies a hard cut to a surface density profile (sigma).  If outer=True,
    sigma = 0 at r > rcut.  Otherwise, sigma = 0 at r < rcut.  If rcut is
    None, inf, or nan no cut is performed.
    """
    
    if rcut is None:
        
        return sigma
        
    elif np.isnan(rcut) or np.isinf(rcut):
        
        return sigma
        
    if outer:
        
        mask = r > rcut
        
    else:
        
        mask = r < rcut
    
    if np.any(mask):
        
        sigma[mask] = 0
        
    return sigma
    
    
def viscous(settings):
    """
    Generates a surface density profile derived from a self-similarity solution
    for a viscous disk, according to:
    
        sigma ~ r^-gamma exp(-r^(2-gamma))
        
    Where r is a dimensionless radius and gamma is a constant less than 2.
    Rd (disk radius) is defined as the radius containing 95% of the disk mass
    
    **ARGUMENTS**
    
    settings : IC settings
        settings like those contained in an IC object (see ICgen_settings.py)
        
    **RETURNS**
    
    R : SimArray
        Radii at which sigma is calculated
    sigma : SimArray
        Surface density profile as a function of R
    """
    Rd = settings.sigma.Rd
    rin = settings.sigma.rin
    rmax = settings.sigma.rmax
    n_points = settings.sigma.n_points
    gamma = settings.sigma.gamma
    m_disk = settings.sigma.m_disk
    
    # Define the fraction of mass contained within Rd
    A = 0.95
    # Normalization for r
    R1 = Rd / (np.log(1/(1-A))**(1/(2-gamma)))
    Rmax = rmax * Rd
    Rin = rin * Rd
    
    R = np.linspace(0, Rmax, n_points)
    r = (R/R1).in_units('1')
    sigma = (r**-gamma) * np.exp(-r**(2-gamma)) * (m_disk/(2*np.pi*R1*R1)) * (2-gamma)   
    # Deal with infinities at the origin with a hard cut off
    sigma[0] = sigma[1]
    
    # Apply interior cutoff
    cut_mask = R < Rin
    if np.any(cut_mask):
        
        sigma[cut_mask] *= smoothstep(r[cut_mask],degree=21,rescale=True)
    
    
    return R, sigma
    
def powerlaw(settings, T = None):
    """
    Generates a surface density profile according to a powerlaw sigma ~ r^p
    with a smooth interior cutoff and smooth exterior exponential cutoff.
    
    **ARGUMENTS**
    
    settings : IC settings
        settings like those contained in an IC object (see ICgen_settings.py)
    T : callable function
        Function that returns temperature of the disk as a function of radius
        IF none, a powerlaw temperature is assumed
    
    **RETURNS**
    
    R : SimArray
        Radii at which sigma is calculated
    sigma : SimArray
        Surface density profile as a function of R
    """
    # Parse settings
    Rd = settings.sigma.Rd
    rin = settings.sigma.rin
    rmax = settings.sigma.rmax
    cutlength = settings.sigma.cutlength
    Mstar = settings.physical.M
    Qmin = settings.sigma.Qmin
    n_points = settings.sigma.n_points
    m = settings.physical.m
    power = settings.sigma.power    

    if T is None:
        # If no callable object to calculate Temperature(R) is provided, 
        # default to a powerlaw T ~ R^-q
        
        T0 = SimArray([129.0],'K') # Temperature at 1 AU
        R0 = SimArray([1.0],'au')
        q = 0.59
        def T(x):
            
            return T0 * np.power((x/R0).in_units('1'),-q)
        
    Rd = match_units(pynbody.units.au, Rd)[1]
    Mstar = match_units(pynbody.units.Msol, Mstar)[1]
    # Molecular weight
    m = match_units(m, pynbody.units.m_p)[0]
    # Maximum R to calculate sigma at (needed for the exponential cutoff region)
    Rmax = rmax*Rd
    
    # Q calculation parameters:
    G = SimArray([1.0],'G')
    kB = SimArray([1.0],'k')
    
    # Initialize stuff
    A = SimArray(1.0,'Msol')/(2*np.pi*np.power(Rd,2))
    R = np.linspace(0,Rmax,n_points)
    r = np.array((R/Rd).in_units('1'))
    
    # Calculate sigma
    # Powerlaw
    #sigma = A/r
    #dflemin3 edit 06/10/2015: Try powerlaw of the form sigma ~ r^power
    sigma = A*np.power(r,power)
    sigma[0] = 0.0
    # Exterior cutoff
    sigma[r>1] *= np.exp(-(r[r>1] - 1)**2 / (2*cutlength**2))
    # Interior cutoff
    sigma[r<rin] *= smoothstep(r[r<rin],degree=21,rescale=True)
    
    # Calculate Q
    Q = np.sqrt(Mstar*kB*T(R)/(G*m*R**3))/(np.pi*sigma)
    Q.convert_units('1')
    
    # Rescale sigma to meet the minimum Q requirement
    sigma *= Q.min()/Qmin
    
    # Calculate Q
    Q = np.sqrt(Mstar*kB*T(R)/(G*m*R**3))/(np.pi*sigma)
    Q.convert_units('1')
    return R, sigma
    
#def MQWS(n_points=1000, rin=4.0, rout=20.0, rmax = None, m_disk=0.1):
def MQWS(settings, T):
    """
    Generates a surface density profile as the per method used in Mayer, Quinn,
    Wadsley, and Stadel 2004
    
    ** ARGUMENTS **
    NOTE: if units are not supplied, assumed units are AU, Msol
    
    settings : IC settings
        settings like those contained in an IC object (see ICgen_settings.py)
        
    T : callable
        A function to calculate temperature as a function of radius
        
    ** RETURNS **
    
    r : SimArray
        Radii at which sigma is calculated
    sigma : SimArray
        Surface density profile as a function of R
    """
    # Q calculation parameters:
    G = SimArray([1.0],'G')
    kB = SimArray([1.0],'k')
    
    # Load in settings
    n_points = settings.sigma.n_points
    rin = settings.sigma.rin
    rout = settings.sigma.rout
    rmax = settings.sigma.rmax
    Qmin = settings.sigma.Qmin
    m = settings.physical.m
    Mstar = settings.physical.M
    #m_disk = settings.sigma.m_disk
    
    rin = match_units(pynbody.units.au, rin)[1]
    rout = match_units(pynbody.units.au, rout)[1]
    #m_disk = match_units(pynbody.units.Msol, m_disk)[1]
    
    if rmax is None:
        
        rmax = 2.5 * rout
        
    else:
        
        rmax = match_units(pynbody.units.au, rmax)[1]
        
    r = np.linspace(0, rmax, n_points)
    
    a = (rin/r).in_units('1')
    b = (r/rout).in_units('1')
    sigma = (np.exp(-a**2 - b**2)/r) * Mstar.units/r.units
    
    # Calculate Q
    Q = np.sqrt(Mstar*kB*T(r)/(G*m*r**3))/(np.pi*sigma)
    Q.convert_units('1')

    sigma *= np.nanmin(Q)/Qmin
    
    # Remove all nans
    sigma[np.isnan(sigma)] = 0.0
    
    
    return r, sigma
