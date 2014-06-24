# -*- coding: utf-8 -*-
"""
Calculates the temperature at a given r for a power law:
    T = T0(r/r0)^Tpower
If r_in has no units, au are assumed

Created on Wed Jan 29 13:20:04 2014

@author: ibackus
"""
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import isaac

class T:
    """
    Calculates temperature as a function of radius.  Updating the settings 
    contained in the ICobject will automatically cause T to be calculated
    accordingly.  ICobject should contain the attribute 'settings'
    
    USAGE:
    
    T = calc_temp.T(ICobject)       # Generate function to calculate T
    
    T(r)        # Calculate T at r
    
    """
    def __init__(self, ICobj):
        
        self._parent = ICobj
        
    def __call__(self, r):
        
        # Load settings
        params = self._parent.settings.physical
        
        if not hasattr(params, 'kind'):
            # Add this check for backwards compatibility.  Previous versions
            # only had one kind of temperature profile
            params.kind = 'powerlaw'
            
        T0 = params.T0
        Tmin = params.Tmin
        r0 = params.r0
        kind = params.kind
        
        # Calculate T(r)
        r = isaac.match_units(r, r0)[0]
        a = (r/r0)
        a = isaac.match_units(a, '1')[0]
        
        if kind.lower() == 'powerlaw':
            
            Tpower = params.Tpower
            Tout = T0 * np.power(a, Tpower)
            Tout[Tout < Tmin] = Tmin
            
        elif kind.lower() == 'mqws':
            
            # NOTE: I'm not sure how exactly they generate the temperature
            # profile.  The quoted equation doesn't match their figures
            Tout = T0 * np.exp(-a)**1.5 + Tmin
        
        if hasattr(params, 'Tmax'):
            Tmax = params.Tmax
            Tout[Tout > Tmax] = Tmax
        
        return Tout