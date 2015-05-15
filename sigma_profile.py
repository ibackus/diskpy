# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 10:17:53 2014

@author: ibackus
"""

import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import isaac

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
    
    return r, sigma
    
def viscous(settings):
    """
    Generates a surface density profile derived from a self-similarity solution
    for a viscous disk, according to:
    
        sigma ~ r^-gamma exp(-r^(2-gamma))
        
    Where r = R/Rdisk is a dimensionless radius and gamma is a constant less
    than 2.
    
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
    #Mstar = settings.physical.M
    n_points = settings.sigma.n_points
    gamma = settings.sigma.gamma
    m_disk = settings.sigma.m_disk
    
    Rmax = rmax * Rd
    
    R = np.linspace(0, Rmax, n_points)
    r = np.linspace(0, rmax, n_points)
    sigma = (r**-gamma) * np.exp(-r**(2-gamma)) * (m_disk/(2*np.pi*Rd*Rd)) * (2-gamma)   
    
    # Apply interior cutoff
    cut_mask = r < rin
    if np.any(cut_mask):
        
        sigma[r<rin] *= isaac.smoothstep(r[r<rin],degree=21,rescale=True)
    
    
    return R, sigma
    
#def powerlaw(Rd=SimArray(1.0,'au'), rin=0.5, rmax=2.3, cutlength=0.3, \
#Mstar=SimArray(1.0/3.0,'Msol'), Qmin=1.5, n_points=1000, m=2.0, T=None):
def powerlaw(settings, T = None):
    """
    Generates a surface density profile according to a powerlaw sigma ~ 1/r
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

    if T is None:
        # If no callable object to calculate Temperature(R) is provided, 
        # default to a powerlaw T ~ R^-q
        
        T0 = SimArray([129.0],'K') # Temperature at 1 AU
        R0 = SimArray([1.0],'au')
        q = 0.59
        def T(x):
            
            return T0 * np.power((x/R0).in_units('1'),-q)
        
    Rd = isaac.match_units(pynbody.units.au, Rd)[1]
    Mstar = isaac.match_units(pynbody.units.Msol, Mstar)[1]
    # Molecular weight
    m = isaac.match_units(m, pynbody.units.m_p)[0]
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
    sigma = A/r
    sigma[0] = 0.0
    # Interior cutoff
    sigma[r>1] *= np.exp(-(r[r>1] - 1)**2 / (2*cutlength**2))
    # Exterior cutoff
    sigma[r<rin] *= isaac.smoothstep(r[r<rin],degree=21,rescale=True)
    
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
    
    rin = isaac.match_units(pynbody.units.au, rin)[1]
    rout = isaac.match_units(pynbody.units.au, rout)[1]
    #m_disk = isaac.match_units(pynbody.units.Msol, m_disk)[1]
    
    if rmax is None:
        
        rmax = 2.5 * rout
        
    else:
        
        rmax = isaac.match_units(pynbody.units.au, rmax)[1]
        
    r = np.linspace(0, rmax, n_points)
    
    #A = m_disk * np.exp(2 * (rin/rout).in_units('1'))/(rout * np.pi**1.5)
    
    a = (rin/r).in_units('1')
    b = (r/rout).in_units('1')
    #sigma = A * np.exp(-a**2 - b**2)/r
    sigma = (np.exp(-a**2 - b**2)/r) * Mstar.units/r.units
    
    # Calculate Q
    Q = np.sqrt(Mstar*kB*T(r)/(G*m*r**3))/(np.pi*sigma)
    Q.convert_units('1')

    sigma *= np.nanmin(Q)/Qmin
    
    # Remove all nans
    sigma[np.isnan(sigma)] = 0.0
    
    
    return r, sigma
