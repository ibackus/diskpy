# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 10:17:53 2014

@author: ibackus
"""

import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import isaac

def powerlaw(Rd=SimArray(1.0,'au'), rin=0.5, rmax=2.3, cutlength=0.3, \
Mstar=SimArray(1.0/3.0,'Msol'), Qmin=1.5, n_points=1000, m=2.0, T=None):
    """
    Generates a surface density profile according to a powerlaw sigma ~ 1/r
    with a smooth interior cutoff and smooth exterior exponential cutoff.
    
    ** ARGUMENTS **
    
    Rd : float or SimArray
        Disk radius at which the exponential cutoff begins.  If float, it is
        assumed to be in units of AU (pynbody.units.au)
    rin : float
        Fraction Rd of at which the interior cutoff happens.  i.e., the radius
        of the interior cutoff in units of Rd
    rmax : float
        Largest radius to calculate sigma at in units of Rd.  This is needed 
        for a smooth exponential cutoff.  Default of 2.3 should work
    cutlength : float
        Length scale for the exponential cut off in units of Rd
    Mstar : float or SimArray
        Mass of the central star.  If a float, it is assumed to be in units of
        Msol
    Qmin : float
        Minimum Toomre Q for the disk.  The surface density is scaled to match
        Qmin.
    n_points : int
        Number of bins to calculate sigma, R at.
    m : float or SimArray
        mean molecular weight.  If float, assumed to be in units of m_p
    T : callable
        Temperature as a function of radius.  If None, it is assumed to be
        a powerlaw T ~ r^-q
    
    
    ** RETURNS **
    
    R : SimArray
        Radii at which sigma is calculated
    sigma : SimArray
        Surface density profile as a function of R
    """

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