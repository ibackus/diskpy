# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 15:39:28 2014

@author: ibackus
"""

__version__ = "$Revision: 1 $"
# $Source$

import numpy as np
import pynbody
SimArray = pynbody.array.SimArray

import isaac
import ICgen_utils

import os
import glob

def v_xy(f, param, changbin=None, nr=50, min_per_bin=100, changa_preset=None,r=None):
    """
    Attempts to calculate the circular velocities for particles in a thin
    (not flat) keplerian disk.  Requires ChaNGa
    
    **ARGUMENTS**
    
    f : tipsy snapshot
        For a gaseous disk
    param : dict
        a dictionary containing params for changa. (see isaac.configparser)
    changbin : str  (OPTIONAL)  
        If set, should be the full path to the ChaNGa executable.  If None, 
        an attempt to find ChaNGa is made
    nr : int (optional)
        number of radial bins to use when averaging over accelerations
    min_per_bin : int (optional)
        The minimum number of particles to be in each bin.  If there are too
        few particles in a bin, it is merged with an adjacent bin.  Thus,
        actual number of radial bins may be less than nr.
    changa_preset : str
        Which ChaNGa execution preset to use (ie 'mpi', 'local', ...).  See
        ICgen_utils.changa_command
    r : SimArray (optional)
        To save on memory, if the cylindrical radius r has already been calculated
        it can be passed to this function
        
    **RETURNS**
    
    vel : SimArray
        An N by 3 SimArray of gas particle velocities.
    """
    
    # Load stuff from the snapshot
    if r is None:
        
        r = f.g['rxy']
        
    cosine = (f.g['x']/r).in_units('1')
    sine = (f.g['y']/r).in_units('1')
    z = f.g['z']
    vel0 = f.g['vel'].copy()
    vel = f.g['vel'].copy()
    
    # Temporary filenames for running ChaNGa
    f_prefix = str(np.random.randint(0, 2**32))
    f_name = f_prefix + '.std'
    p_name = f_prefix + '.param'
    
    # Update parameters
    p_temp = param.copy()
    p_temp['achInFile'] = f_name
    p_temp['achOutName'] = f_prefix
    p_temp['dDelta'] = 1e-10
    if 'dDumpFrameTime' in p_temp: p_temp.pop('dDumpFrameTime')
    if 'dDumpFrameStep' in p_temp: p_temp.pop('dDumpFrameStep')
    
    # --------------------------------------------
    # Estimate velocity from gravity only
    # --------------------------------------------
    # Note, accelerations due to gravity are calculated twice to be extra careful
    # This is so that any velocity dependent effects are properly accounted for
    # (although, ideally, there should be none)
    # The second calculation uses the updated velocities from the first
    for iGrav in range(2):
        # Save files
        f.write(filename=f_name, fmt = pynbody.tipsy.TipsySnap)
        isaac.configsave(p_temp, p_name, ftype='param')
        
        # Run ChaNGa, only calculating gravity
        command = ICgen_utils.changa_command(p_name, changa_preset, changbin, '-gas +n 0')
        print command
        p = ICgen_utils.changa_run(command)
        p.wait()            
    
        # Load accelerations
        acc_name = f_prefix + '.000000.acc2'
        a = isaac.load_acc(acc_name)
        
        # Clean-up
        for fname in glob.glob(f_prefix + '*'): os.remove(fname)
        
        # Calculate cos(theta) where theta is angle above x-y plane
        cos = (r/np.sqrt(r**2 + z**2)).in_units('1')
        # Calculate radial acceleration times r^2
        ar2 = (a[:,0]*cosine + a[:,1]*sine)*r**2
        
        # Bin the data
        r_edges = np.linspace(r.min(), (1+np.spacing(2))*r.max(), nr + 1)
        ind, r_edges = isaac.digitize_threshold(r, min_per_bin, r_edges)
        ind -= 1
        nr = len(r_edges) - 1
        
        r_bins, ar2_mean, err = isaac.binned_mean(r, ar2, binedges=r_edges, \
        weighted_bins=True)
        
        # Fit lines to ar2 vs cos for each radial bin
        m = np.zeros(nr)
        b = np.zeros(nr)    
        
        for i in range(nr):
            
            mask = (ind == i)
            p = np.polyfit(cos[mask], ar2[mask], 1)
            m[i] = p[0]
            b[i] = p[1]
            
        # Interpolate the line fits
        m_spline = isaac.extrap1d(r_bins, m)
        b_spline = isaac.extrap1d(r_bins, b)
        
        # Calculate circular velocity
        ar2 = SimArray(m_spline(r)*cos + b_spline(r), ar2.units)
        v_calc = (np.sqrt(abs(ar2)/r)).in_units(vel.units)
        vel[:,0] = -v_calc*sine
        vel[:,1] = v_calc*cosine
        
        # Assign to f
        f.g['vel'] = vel
        
    # --------------------------------------------
    # Estimate pressure/gas dynamics accelerations
    # --------------------------------------------
    a_grav = a
    ar2_calc_grav = ar2
    
    # Save files
    f.write(filename=f_name, fmt = pynbody.tipsy.TipsySnap)
    isaac.configsave(p_temp, p_name, ftype='param')
    
    # Run ChaNGa, including SPH
    command = ICgen_utils.changa_command(p_name, changa_preset, changbin, '+gas -n 0')
    p = ICgen_utils.changa_run(command)
    p.wait()
        
    # Load accelerations
    acc_name = f_prefix + '.000000.acc2'
    a_total = isaac.load_acc(acc_name)
    
    # Clean-up
    for fname in glob.glob(f_prefix + '*'): os.remove(fname)
    
    # Estimate the accelerations due to pressure gradients/gas dynamics
    a_gas = a_total - a_grav
    ar2_gas = (a_gas[:,0]*cosine + a_gas[:,1]*sine)*r**2
    
    logr_bins, ratio, err = isaac.binned_mean(np.log(r), ar2_gas/ar2_calc_grav, nbins=nr,\
    weighted_bins=True)
    r_bins = np.exp(logr_bins)
    ratio_spline = isaac.extrap1d(r_bins, ratio)
    
    ar2_calc = ar2_calc_grav*(1 + ratio_spline(r))
    
    v = (np.sqrt(abs(ar2_calc)/r)).in_units(vel0.units)
    
    vel[:,0] = -v*sine
    vel[:,1] = v*cosine
    
    # more cleanup
    f.g['vel'] = vel0
    
    return vel