#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This defines derived arrays and utility functions to extend the utility of
tipsy SimSnaps.  Note that this has only been tested on tipsy formatted
snapshots.

The ONLY functions that should be defined OR defined as an import in this file
should be deriveda arrays.  This is to make generating documentation on the
derived arrays simple.

i.e. don't do from module import func.  Just do import module and use 
module.func

Created on Sat Sep 16 16:17:20 2017

@author: ibackus
"""


import numpy as np
import pynbody
SimArray = pynbody.array.SimArray

import warnings

import _simsnaputils as sutil

kB = SimArray(1.0,'k')
G = SimArray(1.0, 'G')


# -------------------------------------------------------------------
# Array wrappers to handle loading aux arrays with units
# -------------------------------------------------------------------

@pynbody.derived_array
def u_dustFrac(f):
    """
    Returns the auxilliary array of the same name (excluding the 'u_' prefix)
    with units inferred from the runtime parameters.
    """
    array_name = 'dustFrac'
    units = '1'
    if not pynbody.units.has_units(f[array_name]):
        f[array_name].units = units
    return f[array_name]

@pynbody.derived_array
def u_dustFracDot(f):
    """
    Returns the auxilliary array of the same name (excluding the 'u_' prefix)
    with units inferred from the runtime parameters.
    """
    array_name = 'dustFracDot'
    
    if not pynbody.units.has_units(f[array_name]):
        # Infer units
        units = 1/sutil.get_snap_unit(f, 't_unit')
        f[array_name].units = units
    return f[array_name]

@pynbody.derived_array
def u_smoothlength(f):
    """
    Returns the auxilliary array of the same name (excluding the 'u_' prefix)
    with units inferred from the runtime parameters.
    """
    array_name = 'smoothlength'
    
    if not pynbody.units.has_units(f[array_name]):
        # Infer units
        units = sutil.get_snap_unit(f, 'l_unit')
        f[array_name].units = units
    return f[array_name]

@pynbody.derived_array
def u_dustVel(f):
    """
    Returns the auxilliary array of the same name (excluding the 'u_' prefix)
    with units inferred from the runtime parameters.
    """
    array_name = 'dustVel'
    
    if not pynbody.units.has_units(f[array_name]):
        # Infer units
        units = sutil.get_snap_unit(f, 'v_unit')
        f[array_name].units = units
    
    # This is a necessary check -- currently the master branch of pynbody
    # cannot load 3D arrays.  This has been implemented on branch issue-379
    # of ibackus' fork, but the pull request has not been merged
    if (np.ndim(f[array_name]) != 2) or (f[array_name].shape[1] != 3):
        raise RuntimeError, "dustVel is a 3D array but was not read as one."\
        "  This is problem with pynbody"
    return f[array_name]

# -------------------------------------------------------------------
# Derived arrays
# -------------------------------------------------------------------
@pynbody.derived_array
def roche_rho(f):
    r"""
    Roche density
    
    .. math::
        \rho_{roche} \equiv M_{star}/r^3
    
    where :math:`r` is the distance to the center of mass of the star(s).  The
    total star mass is used.
    """
    parent = sutil.get_parent_snap(f)
    if parent != f:
        warnings.warn("Accessing stars from parent snapshot")
    
    star_mass = parent.s['mass'].sum()
    star_pos = pynbody.analysis.halo.center_of_mass(parent.s)
    r_star = np.sqrt(((f['pos'] - star_pos)**2).sum(1))
    rho = 16. * star_mass/r_star**3
    rho_unit = sutil.get_snap_unit(parent, 'rho_unit')
    return rho.in_units(rho_unit)

@pynbody.derived_array
def roche_rho_fac(f):
    """
    Gas density divided by roche density
    """
    return (f['rho']/f['roche_rho']).in_units('1')

@pynbody.derived_array
def cs(f):
    """
    Sound speed
    """
    if sutil.is_isothermal(f):
        gamma = 1.
    else:
        gamma = sutil.get_snap_gamma(f)
    m = sutil.get_snap_mu(f)
    v_unit = sutil.get_snap_unit(f, 'v_unit')
    c = np.sqrt(kB * f['temp']*gamma/m).in_units(v_unit)
    return c

@pynbody.derived_array
def P(f):
    """
    Pressure
    """
    pres = f['rho'] * f['cs']**2
    if 'dustFrac' in f.all_keys():
        pres *= (1 - f['dustFrac'])
    pres_unit = sutil.get_snap_unit(f, 'pres_unit')
    return pres.in_units(pres_unit)

@pynbody.derived_array
def tstop(f):
    """
    Dust stopping time
    """
    units = sutil.get_all_units(f)
    grainSize = SimArray(sutil.get_snap_param(f, 'dDustSize'), units['l_unit'])
    grainDensity = SimArray(sutil.get_snap_param(f, 'dDustGrainDensity'), units['rho_unit'])
    
    if sutil.is_isothermal(f):
        gamma = 1.
    else:
        gamma = sutil.get_snap_gamma(f)
        
    t = ((grainSize*grainDensity)/(f['rho'] * f['cs'])) * np.sqrt(np.pi*gamma/8.)
    return t.in_units(units['t_unit'])

@pynbody.derived_array
def stokes(f):
    """
    Dust stokes parameter, defined as the dust stopping time times the 
    Keplerian orbital angular velocity
    
    R is calculated relative to the origin!  NOT THE STARS!
    """
    parent = sutil.get_parent_snap(f)
    Mstar = parent.s['mass'].sum()
    tstop = f['tstop']
    omega = np.sqrt(G*Mstar/f['rxy']**3)
    
    if pynbody.units.has_units(tstop) and pynbody.units.has_units(omega):
        return (tstop*omega).in_units('1')
    else:
        return tstop*omega
    
@pynbody.derived_array
def rho_dust(f):
    """
    Dust density
    """
    return f['u_dustFrac'] * f['rho']


# Velocities ---------------------
@pynbody.derived_array
def dust_delta_v(f):
    """
    Dust delta v is defined as v_dust - v_gas
    """
    # Note: I add a small number to the denominator to avoid divide by zero
    # as dustFrac -> 1, the gas velocity is not well defined and delta v is
    # therefore poorly defined
    return (f['u_dustVel'] - f['vel'])/(1 - f['u_dustFrac'][:,None] + 1e-15)

@pynbody.derived_array
def dust_delta_v_phi(f):
    """
    Cylindrical tangential component of the dust delta v velocity in the x-y 
    plane.  
    Dust delta v is defined as v_dust - v_gas
    """
    return sutil.polar_phi_hat(f, 'dust_delta_v')

@pynbody.derived_array
def dust_delta_v_rxy(f):
    """
    Radial dust delta v in x-y plane.
    Dust delta v is defined as v_dust - v_gas
    """
    return sutil.polar_r_hat(f, 'dust_delta_v')


@pynbody.derived_array
def dustVel_phi(f):
    """
    Cylindrical tangential component of the dust velocity in the x-y plane
    """
    return sutil.polar_phi_hat(f, 'u_dustVel')

@pynbody.derived_array
def dustVel_rxy(f):
    """
    Radial dust velocity in the x-y plane
    """
    return sutil.polar_r_hat(f, 'u_dustVel')

@pynbody.derived_array
def gas_vel(f):
    """
    Gas velocity
    """
    v = f['vel']
    if 'dustFrac' in f.all_keys():
        v -= f['u_dustFrac'][:, None] * f['u_dustVel']
        v /= (1 - f['u_dustFrac'][:,None] + 1e-15)
    return v

@pynbody.derived_array
def gas_v_rxy(f):
    """
    Gas cylindrical radial velocity
    """
    return sutil.polar_r_hat(f, 'gas_vel')

@pynbody.derived_array
def gas_v_phi(f):
    """
    Cylindrical tangential component of the gas velocity in the x-y plane
    """
    return sutil.polar_phi_hat(f, 'gas_vel')

@pynbody.derived_array
def dustmass(f):
    """
    dustFrac * mass
    """
    return f['u_dustFrac'] * f['mass']

# -------------------------------------------
# analysis arrays
# -------------------------------------------

@pynbody.derived_array
def abs_z(f):
    """
    Absolute value of z
    """
    return abs(f['z'])

@pynbody.derived_array
def dust_midplane_v(f):
    """
    component of delta v (v_dust - v_gas) in the direction of the midplane
    """
    sign = np.sign(f['z'])
    sign.units = '1'
    return -sign * f['dust_delta_v'][:,2]

@pynbody.derived_array
def scaleheight(f):
    """
    |z|/r (in cylindrical coords)
    """
    return (abs(f['z'])/f['rxy']).in_units('1')

@pynbody.derived_array
def abs_dust_rho_dot(f):
    """
    abs(rho * dustFracDot)
    """
    return abs(dust_rho_dot(f))

@pynbody.derived_array
def dust_rho_dot(f):
    """
    dustFracDot * rho
    """
    return f['u_dustFracDot'] * f['rho']
