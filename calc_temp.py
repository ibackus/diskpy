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
    
    REQUIRED SETTINGS:
    
    physical.T0
    physical.r0
    physical.Tpower
    physical.Tmin
    
    """
    def __init__(self, ICobj):
        
        self._parent = ICobj
        
    def __call__(self, r):
        
        # Load settings
        params = self._parent.settings.physical
        T0 = params.T0
        Tpower = params.Tpower
        r0 = params.r0
        Tmin = params.Tmin
        
        if hasattr(params, 'Tmax'):
            Tmax = params.Tmax
        # Calculate T(r)
        r = isaac.match_units(r, r0)[0]
        a = (r/r0)
        a = isaac.match_units(a, '1')[0]
        Tout = T0 * np.power(a, Tpower)
        Tout[Tout < Tmin] = Tmin
        if hasattr(params, 'Tmax'):
            Tout[Tout > Tmax] = Tmax
        
        return Tout