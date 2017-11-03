#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 17:01:21 2017

@author: ibackus
"""
import pynbody
SimArray = pynbody.array.SimArray
import diskpy

# ----------------------------------------------------------------------
# Utility function
# ----------------------------------------------------------------------
def get_parent_snap(f):
    """
    Find the parent (last ancestor) of SubSnap f
    """
    parent = f.ancestor
    while parent != f:
        f = parent
        parent = f.ancestor
    
    return parent

def get_parent_param_dict(f):
    """
    Get the param dict from parent of f.  If it doesn't exist, try to infer
    it from parent._paramfile
    """
    parent = get_parent_snap(f)
    if not hasattr(parent, 'param'):
        parent.param = diskpy.utils.snap_param(parent)
        
    return parent.param

def get_all_units(f):
    """
    Retrieve a dict of units from a tipsy SimSnap
    """
    parent = get_parent_snap(f)
    if not hasattr(parent, 'units'):
        parent_dict = get_parent_param_dict(parent)
        parent.units = diskpy.pychanga.units_from_param(parent_dict)
    
    return parent.units

def get_snap_unit(f, unit):
    """
    Attempts to get the simulation units for unit
    
    To see available units, try get_all_units(f)
    """
    units = get_all_units(f)
    return units[unit]

def get_snap_param(f, key, use_defaults=False):
    """
    Get the runtime param from a snapshot.  IF use_defaults, keys not present
    in the .param file will be substituted with defaults
    """
    param = get_parent_param_dict(f)
    if use_defaults:
        return diskpy.pychanga.getpar(key, param)
    else:
        return param[key]


def get_snap_gamma(f):
    """
    Retrieve the adiabatic index for a simulation
    """
    gamma = get_snap_param(f, 'dConstGamma', use_defaults=True)
    return gamma

def get_snap_mu(f):
    """
    Retrieve the mean molecular weight for a tipsy SimSnap
    """
    m = get_snap_param(f, 'dMeanMolWeight', use_defaults=True)
    return SimArray(m, 'm_p', dtype=float)

def is_isothermal(f):
    """
    Check to see if a SimSnap is isothermal
    """
    param = get_parent_param_dict(f)
    if 'bGasIsothermal' in param:
        isothermal = (param['bGasIsothermal'] == 1)
    elif 'bGasAdiabatic' in f.param:
        isothermal = (param['bGasAdiabatic'] != 0)
    else:
        isothermal = False
    return isothermal

def polar_phi_hat(f, array):
    """
    Gets the phi-hat (in polar coords) component of a 2D or 3D cartesian array.
    """
    val = f[array]
    if val.shape[1] not in (2, 3):
        raise ValueError, "array to get phi-hat component of must be 2D or 3D"
    return (f['x']*val[:,1] - f['y'] * val[:,0])/f['rxy']

def polar_r_hat(f, array):
    """
    Gets the radial (in polar coords) component of a 2D or 3D cartesian array.
    """
    val = f[array]
    if val.shape[1] not in (2, 3):
        raise ValueError, "array to get phi-hat component of must be 2D or 3D"
    return (f['x'] * val[:,0] + f['y'] * val[:,1])/f['rxy']
