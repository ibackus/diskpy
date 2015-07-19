# -*- coding: utf-8 -*-
"""
Same as calc_velocity.py, but calls mpi with changa to allow many nodes
NOTE.  mpirrun must be already loaded.  Also, should do export MX_RCACHE=0
before loading python

Created on Wed Apr  9 15:39:28 2014

@author: ibackus
"""

import numpy as np
import pynbody
SimArray = pynbody.array.SimArray

import isaac

import subprocess
import os
import glob

import time

def v_xy(f, param, changbin=None, nr=50, min_per_bin=100):
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
        
    **RETURNS**
    
    vel : SimArray
        An N by 3 SimArray of gas particle velocities.
    """
    
    if changbin is None:
        # Try to find the ChaNGa binary full path
        changbin = os.popen('which ChaNGa_uw_mpi').read().strip()
        
    # Load up mpi
    
    # Load stuff from the snapshot
    x = f.g['x']
    y = f.g['y']
    z = f.g['z']
    r = f.g['rxy']
    vel0 = f.g['vel'].copy()
    
    # Remove units from all quantities
    r = isaac.strip_units(r)
    x = isaac.strip_units(x)
    y = isaac.strip_units(y)
    z = isaac.strip_units(z)
    
    # Temporary filenames for running ChaNGa
    f_prefix = str(np.random.randint(0, 2**32))
    f_name = f_prefix + '.std'
    p_name = f_prefix + '.param'
    
    # Update parameters
    p_temp = param.copy()
    p_temp['achInFile'] = f_name
    p_temp['achOutName'] = f_prefix
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
        command = 'mpirun --mca mtl mx --mca pml cm ' + changbin + ' -gas -n 0 ' + p_name
        #command = 'charmrun ++local ' + changbin + ' -gas -n 0 ' + p_name
        
        p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        
        while p.poll() is None:
            
            time.sleep(0.1)
            
    
        # Load accelerations
        acc_name = f_prefix + '.000000.acc2'
        a = isaac.load_acc(acc_name)
        
        # Clean-up
        for fname in glob.glob(f_prefix + '*'): os.remove(fname)
        
        # If a is not a vector, calculate radial acceleration.  Otherwise, assume
        # a is the radial acceleration                
        a_r = a[:,0]*x/r + a[:,1]*y/r
        
        # Make sure the units are correct then remove them
        a_r = isaac.match_units(a_r, a)[0]
        a_r = isaac.strip_units(a_r)
        
        # Calculate cos(theta) where theta is angle above x-y plane
        cos = r/np.sqrt(r**2 + z**2)
        ar2 = a_r*r**2
        
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
        ar2_calc = m_spline(r)*cos + b_spline(r)
        v_calc = np.sqrt(abs(ar2_calc)/r)
        vel = f.g['vel'].copy()
        v_calc = isaac.match_units(v_calc,vel)[0]
        vel[:,0] = -v_calc*y/r
        vel[:,1] = v_calc*x/r
        
        # Assign to f
        f.g['vel'] = vel
        
    # --------------------------------------------
    # Estimate pressure/gas dynamics accelerations
    # --------------------------------------------
    a_grav = a
    ar2_calc_grav = ar2_calc
    
    # Save files
    f.write(filename=f_name, fmt = pynbody.tipsy.TipsySnap)
    isaac.configsave(p_temp, p_name, ftype='param')
    
    # Run ChaNGa, including SPH
    command = 'mpirun --mca mtl mx --mca pml cm ' + changbin + ' +gas -n 0 ' + p_name
    
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    
    while p.poll() is None:
        
        time.sleep(0.1)
        
    # Load accelerations
    acc_name = f_prefix + '.000000.acc2'
    a_total = isaac.load_acc(acc_name)
    
    # Clean-up
    for fname in glob.glob(f_prefix + '*'): os.remove(fname)
    
    # Estimate the accelerations due to pressure gradients/gas dynamics
    a_gas = a_total - a_grav
    ar_gas = a_gas[:,0]*x/r + a_gas[:,1]*y/r
    ar_gas = isaac.strip_units(ar_gas)
    ar2_gas = ar_gas*r**2
    
    logr_bins, ratio, err = isaac.binned_mean(np.log(r), ar2_gas/ar2_calc_grav, nbins=nr,\
    weighted_bins=True)
    r_bins = np.exp(logr_bins)
    ratio_spline = isaac.extrap1d(r_bins, ratio)
    
    ar2_calc = ar2_calc_grav*(1 + ratio_spline(r))
    a_calc = ar2_calc/r**2
    
    v = np.sqrt(r*abs(a_calc))
    v = isaac.match_units(v, vel0.units)[0]
    vel = vel0.copy()
    vel[:,0] = -v*y/r
    vel[:,1] = v*x/r
    
    # more cleanup
    f.g['vel'] = vel0
    
    return vel