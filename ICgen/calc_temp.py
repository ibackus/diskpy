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

from diskpy.utils import strip_units, match_units

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
        self.adiabatic_ready = False
        self.kind = ICobj.settings.physical.kind
        
    def __call__(self, r):
        
        # Load settings
        params = self._parent.settings.physical
        
        # Calculate the basic temperature profile
        Tout = self.T_nocut(r)
        
        # Handle an adiabatic equation of state
        # The interior cutoff affects the choice of temperature profile
        if (params.eos == 'adiabatic') and self.adiabatic_ready:
            
            if not isinstance(r, np.ndarray):
                
                r = np.array(r)
                
            nr = np.product(r.shape)
            mask = r < self.r_a
            
            if np.any(mask) and (nr > 1):
                
                Tout[mask] = self.T_adiabatic(r[mask])
                
            elif np.any(mask):
                
                Tout = self.T_adiabatic(r)
        
        return Tout
        
    def _est_r_adiabatic(self):
        """
        Estimates the radius at which the temperature profile should be made
        adiabatic (for adiabatic disks only!).  This is defined the be the
        first radius at which the entropy gradient switches from negative
        to positive
        
        Returns
        -------
        
        r_adiabatic : SimArray
            radius
        
        Notes
        -----
        
        Entropy gradients are estimated assuming:
        
        * Star gravity dominates disk self-gravity
        
        * Disk is thin
        
        * No vertical temperature gradient
        
        """
        ICobj = self._parent
        
        if not hasattr(ICobj, 'sigma'):
            
            raise RuntimeError, 'Sigma not found in initial conditions object'
            
        gamma = ICobj.settings.physical.gamma
        p = (gamma - 1.)/(gamma + 1.)
        a = 1/p        
        r = ICobj.sigma.r_bins
        temp = strip_units(self.T_nocut(r))
        sigma = strip_units(ICobj.sigma(r))
        r1 = strip_units(r)
        # Entropy/mass (up to a multiplicative and an additive constant)
        s = 0.5*a*np.log(temp) + 1.5*np.log(r1) - np.log(sigma)
        # Entropy gradient (up to a multiplicative constant)
        ds = np.gradient(s)
        # Find all negative to positive zero crossings
        neg = ds < 0
        ind = np.where(neg[0:-1] & ~neg[1:])[0] + 1
        
        if len(ind) > 1:
            
            print 'WARNING: Multiple r_adiabatics found.  Selecting the lowest'
            r_adiabatic = r[[ind[0]]]
            
        elif len(ind) < 1:
            
            print 'WARNING: could not find r_adiabatic'
            r_adiabatic = SimArray([0.],'au')
            
        else:
            
            r_adiabatic = r[ind]
            
        return r_adiabatic
        
    def _est_Tscale(self):
        """
        Estimates the temperature scaling factor for the adiabatic temperature
        profile so that the adiabatic temperature matches the user-set temp
        at r_adiabatic.
        """
        
        self.Tscale = 1.
        
        T0 = self.T_nocut(self.r_adiabatic)
        T1 = self.T_adiabatic(self.r_adiabatic)
        
        A = T0/T1
#        p = self.p
#        
#        r = self.r_adiabatic.in_units('au')
#        T0 = self.T_nocut(r)
#        r = strip_units(r)
#        sigma = self._parent.sigma(r)
#        sigma.convert_units('Msol au**-2')
#        sigma = strip_units(sigma)
#        
#        A = T0 * (sigma**2/r**3)**-p
        return A
        
            
    def T_nocut(self, r):
        """
        Calculates the temperature as a function of r ignoring any interior
        cutoff.  The interior cutoff would affect T for adiabatic disks
        
        Paramters
        ---------
        
        r : array, SimArray, or float
            radius (or radii) at which to calculate the temperature
        
        Returns
        -------
        
        T : SimArray
            temperature
        """
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
        r = match_units(r, r0)[0]
        a = (r/r0)
        a = match_units(a, '1')[0]
        
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
            
            raise TypeError, 'Could not find temperature kind {0}'.format(kind)
            
        # Apply max Temp cutoff
        if hasattr(params, 'Tmax'):
            Tmax = params.Tmax
            Tout[Tout > Tmax] = Tmax
        
        return Tout
    
    def setup_interior(self):
        """
        Sets up the interior temperature profile.
        
        For an adiabatic disk, this is an adiabatic profile inside of the
        cutoff.  Otherwise, nothing is done.
        
        Requires that sigma already be calculated
        """
        
        ICobj = self._parent
        
        if ICobj.settings.physical.eos == 'adiabatic':
            
            if not hasattr(ICobj, 'sigma'):
                
                raise RuntimeError('Cannot setup interior temp profile: sigma '
                'needs to be calculated first')
                
#            gamma = ICobj.settings.physical.gamma
#            p = (gamma - 1.)/(gamma + 1.)
#            
#            self.p = p
#            self.gamma = gamma
#            self.r_adiabatic = self._est_r_adiabatic()
#            self.Tscale = self._est_Tscale()
            
            r_a = self._est_r_adiabatic()
            self.sigma_a = ICobj.sigma(r_a)
            self.T_a = self.T_nocut(r_a)
            self.r_a = r_a
            self.adiabatic_ready = True
            
        else:
            
            pass
        
    def T_adiabatic(self, r):
        """
        Estimates the adiabatic temperature profile as a function of r
        
        setup_interior must be run first
        """
#        if not self.adiabatic_ready:
#            
#            # Try to setup the adiabatic profile
#            self.setup_interior()
            
#        p = self.p
#        A = self.Tscale
#        print A
#        
#        r = match_units(r, 'au')[0]
#        sigma = self._parent.sigma(r)
#        sigma.convert_units('Msol au**-2')
#        sigma = strip_units(sigma)
#        r = strip_units(r)
#        return A * ((sigma**2)/(r**3))**p
            
        b = 1.5
        
        sigma = self._parent.sigma(r)
        r_a = self.r_a
        r = match_units(r, r_a.units)[0]
        x = (r/r_a).in_units('1')
        sigma_a = self.sigma_a
        y = (sigma/sigma_a).in_units('1')
        
        x = strip_units(x)
        y = strip_units(y)
        
        gamma = self._parent.settings.physical.gamma
        p = (gamma-1.)/(gamma+1.)
        
        T_a = self.T_a
        
        T = T_a * (y**(2*p)) * (x**-b)
        
        return T
        
#        r = strip_units(match_units(r,'au')[0])
#        
#        return A * (r**2)