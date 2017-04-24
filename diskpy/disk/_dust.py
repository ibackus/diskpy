#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:18:34 2017

@author: ibackus
"""

import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import warnings

from diskpy.pychanga import units_from_param
from diskpy.utils import as_simarray

def setDustfrac(f, dustFrac=0.01, scale_mass=False):
    """
    Set the dust fraction to a SimSnap or sub-snap, and optionally scale the
    mass to keep the original mass (e.g. star mass or gas mass) constant.
    
    Parameters
    ----------
    f : SimSnap
        Snapshot or sub-snap
    dustFrac : number or array-like
        Dust fraction to apply
    scale_mass : bool
        If true, the mass will be scaled to keep the non-dust mass constant.
    """
    if np.any((dustFrac == 1.0 )) and scale_mass:
        
        raise ValueError, "Cannot scale mass and set dust fraction to 1.0"\
            " (divide by zero)"
        
    if len(f) > 0 and dustFrac is not None:
        
        dustFrac0 = f.get('dustFrac', np.array(0.0)).copy()
        f['dustFrac'] = dustFrac
        
        if scale_mass:
            
            mScale = (1.0 - dustFrac0)/(1.0 - dustFrac)
            
            if (np.size(mScale) > 1) and np.any(mScale != mScale[0]):
                
                warnings.warn("Scaling mass by non-uniform values")
                
            f['mass'] *= mScale
    
    return

def add_dust_sim(f, param, gasDustFrac=0.01, starDustFrac=0., 
                 scaleGasMass=True, scaleStarMass=True, grainSize='1 mm',
                 grainDensity='3 g cm**-3'):
    """
    Add a single size of dust to a simulation.  Changes are applied in-place
    to the snapshot and the param dict.
    
    Parameters
    ----------
    f : SimSnap
        Simulation snapshot
    param : dict
        Param dictionary (see diskpy.utils.configparser).
    gasDustFrac, starDustFrac : number or array-like
        Sets the gas fraction for the gas, star.  If None, no dust fraction
        is set.
    scaleGasMass, scaleStarMass : bool
        Scale the gas, star mass to keep it constant after adding dust.
        For example, this can be useful to keep the gas pressure constant.
    grainSize : number, SimArray, str
        The grain size (i.e. radius).  Assumed units (if none are supplied) are
        'mm'
    grainDensity : number, SimArray, str
        Intrinsic grain density (e.g. for silica).  Assumed units if non are
        supplied are g/cm^3
    """
    units = units_from_param(param)
    grainSize = as_simarray(grainSize, 'mm')
    grainDensity = as_simarray(grainDensity, 'g cm**-3')
    grainSize.convert_units(units['l_unit'])
    grainDensity.convert_units(units['rho_unit'])
    if len(f.g) > 0:
        setDustfrac(f.g, gasDustFrac, scaleGasMass)
    if len(f.s) > 0:
        setDustfrac(f.s, starDustFrac, scaleStarMass)
    param['dDustSize'] = float(grainSize)
    param['dDustGrainDensity'] = float(grainDensity)
