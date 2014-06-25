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
    accordingly.  ICobject should contain the attribute 'settings'.  Settings
    for temperature are contained in settings.physical
    
    USAGE:
    
    T = calc_temp.T(ICobject)       # Generate function to calculate T
    
    T(r)        # Calculate T at r
    
    There are multiple kinds of available temperature profiles.  They are set
    in ICobject.settings.physical.kind.  They are:
    
    'powerlaw' : (default)
        Follows a power law, T = T0 * (r/r0)**Tpower
        Settings needed:
            T0
            r0
            Tpower
            Tmin (optional)
            Tmax (optional)
            
    'MQWS'
        Mimics Mayer, Quinn et. all 2004 ICs.  Settings needed:
            T0
            r0 (same as r_in in the paper)
            Tmin
            Tmax (optional)
    
    """
    def __init__(self, ICobj):
        
        self._parent = ICobj
        self.kind = ICobj.settings.physical.kind
        
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
        
        # Powerlaw temperature (default)
        if kind == 'powerlaw':
            
            Tpower = params.Tpower
            Tout = T0 * np.power(a, Tpower)
            Tout[Tout < Tmin] = Tmin
            
        # MQWS temperature profile
        elif (kind == 'mqws') | (kind == 'MQWS'):
            
            # NOTE: I'm not sure how exactly they generate the temperature
            # profile.  The quoted equation doesn't match their figures
            Tout = T0 * np.exp(-3*a/2) + Tmin
            
        else:
            
            raise TypeError, 'Could not find temperature kind {}'.format(kind)
        
        if hasattr(params, 'Tmax'):
            Tmax = params.Tmax
            Tout[Tout > Tmax] = Tmax
            
        
        return Tout        
        