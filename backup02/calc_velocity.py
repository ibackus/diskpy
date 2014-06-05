# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 15:39:28 2014

@author: ibackus
"""

import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import matplotlib.pyplot as plt
import scipy.interpolate as interp

import isaac

import subprocess
import os
import glob

import time

def v_xy(snapshot, param, changbin=None):
    """
    Attempts to calculate the circular velocities for particles in a thin
    (not flat) keplerian disk.  Requires ChaNGa
    
    ARGUMENTS:
    
    snapshot:   a tipsy snapshot for a gaseous disk
    param:      a dictionary containing params for changa.
                (see isaac.configparser)
    changbin:   (OPTIONAL)  If set, should be the full path to the ChaNGa
                executable.  If None, an attempt to find ChaNGa is made
    """
    
    if changbin is None:
        # Try to find the ChaNGa binary full path
        changbin = os.popen('which ChaNGa').read().strip()
    # --------------------------------------------
    # Initialize
    # --------------------------------------------
    # Load things from snapshot
    x = snapshot.g['x']
    y = snapshot.g['y']
    z = snapshot.g['z']
    r = snapshot.g['rxy']
    v_initial = snapshot.g['vel'].copy()
    
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
    # Save files
    snapshot.write(filename=f_name, fmt = pynbody.tipsy.TipsySnap)
    isaac.configsave(p_temp, p_name, ftype='param')
    
    # Run ChaNGa, only calculating gravity
    command = 'charmrun ++local ' + changbin + ' -gas -n 0 ' + p_name
    
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    
    while p.poll() is None:
        
        time.sleep(0.1)
        
    # Load accelerations
    acc_name = f_prefix + '.000000.acc2'
    a = isaac.load_acc(acc_name)
    
    # Calculate radial acceleration
    a_r = a[:,0]*x/r + a[:,1]*y/r
    
    # circular velocity
    v = np.sqrt(abs(r*a_r))
    
    # Assign
    snapshot.g['vel'][:,0] = -v*y/r
    snapshot.g['vel'][:,1] = v*x/r
    
    # --------------------------------------------
    # Estimate velocity from gravity & sph (assuming bDoGas=1 in param)
    # --------------------------------------------
    # Write file
    snapshot.write(filename=f_name, fmt = pynbody.tipsy.TipsySnap)
    
    # Run ChaNGa, only calculating gravity
    changbin = os.popen('which ChaNGa').read().strip()
    #command = 'charmrun ++local ' + changbin + ' +gas -n 0 ' + p_name
    command = 'charmrun ++local ' + changbin + ' -n 0 ' + p_name
    
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    
    while p.poll() is None:
        
        time.sleep(0.1)
        
    # Load accelerations
    acc_name = f_prefix + '.000000.acc2'
    a = isaac.load_acc(acc_name)
    
    # Calculate radial acceleration
    a_r = a[:,0]*x/r + a[:,1]*y/r
    
    # circular velocity
    v = np.sqrt(abs(r*a_r))
    
    logv = np.log(v)
    logr = np.log(r)
    
    logr_bins, logv_mean, err = isaac.binned_mean(logr, logv, 50, \
    weighted_bins=True)
    
    logv_spline = isaac.extrap1d(logr_bins, logv_mean)
    
    v_calc = np.exp(logv_spline(logr))
    v_calc = isaac.match_units(v_calc, v_initial.units)[0]
    v_out = v_initial * 0.0
    v_out[:,0] = -v_calc * y/r
    v_out[:,1] = v_calc * x/r
    
    # Clean-up
    for fname in glob.glob(f_prefix + '*'): os.remove(fname)
    
    # Re-assign velocity to snapshot
    snapshot.g['vel'] = v_initial
    
    return v_out