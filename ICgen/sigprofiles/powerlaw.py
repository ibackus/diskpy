# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:23:50 2015

@author: ibackus
"""


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
    power = settings.sigma.power    

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
    #sigma = A/r
    #dflemin3 edit 06/10/2015: Try powerlaw of the form sigma ~ r^power
    sigma = A*np.power(r,power)
    sigma[0] = 0.0
    # Exterior cutoff
    sigma[r>1] *= np.exp(-(r[r>1] - 1)**2 / (2*cutlength**2))
    # Interior cutoff
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