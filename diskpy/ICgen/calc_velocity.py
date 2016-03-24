# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 15:39:28 2014

@author: ibackus
"""

__version__ = "$Revision: 2 $"
# $Source$

# External packages
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import os
import glob
import gc

# diskpy packages
import ICgen_utils
from diskpy.utils import configsave
from diskpy.pdmath import extrap1d, digitize_threshold, binned_mean
from diskpy.pychanga import load_acc

def v_xy(f, param, changbin=None, nr=50, min_per_bin=100, changa_preset=None, \
max_particles=None, est_eps=True):
    """
    Attempts to calculate the circular velocities for particles in a thin
    (not flat) keplerian disk.  Also estimates gravitational softening (eps)
    for the gas particles and a reasonable time step (dDelta)
    
    Requires ChaNGa
    
    Note that this will change the velocities IN f
    
    Parameters
    ----------
    f : tipsy snapshot
        For a gaseous disk
    param : dict
        a dictionary containing params for changa. (see configparser)
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
    max_particles : int or None
        Specifies the maximum number of particles to use for calculating
        accelerations and velocities.  Setting a smaller number can speed up
        computation and save on memory but can yield noisier results.
        If None, max is unlimited.
    est_eps : bool
        Estimate eps (gravitational softening length).  Default is True.
        If False, it is assumed eps has already been estimated
        
    Returns
    -------
    dDelta : float
        A reasonable time step for the simulation (in code units).  Velocties
        and eps are updated in the snapshot.
    """
    # If the snapshot has too many particles, randomly select gas particles
    # To use for calculating velocity and make a view of the snapshot
    n_gas = len(f) - 1
    subview = (n_gas > max_particles) and (max_particles is not None)
    if subview:
        
        max_particles = int(max_particles)
        mask = np.zeros(n_gas + 1, dtype=bool)
        mask[-1] = True # Use the star particle always
        # randomly select particles to use
        m = np.random.rand(n_gas)
        ind = m.argsort()[0:max_particles]
        mask[ind] = True
        # Make a subview and create a reference to the complete snapshot
        complete_snapshot = f
        f = complete_snapshot[mask]
        # Scale gas mass
        m_scale = float(n_gas)/float(max_particles)
        f.g['mass'] *= m_scale
        
        if not est_eps:
            
            f.g['eps'] *= m_scale**(1.0/3)
        
    # Load stuff from the snapshot
    r = f.g['rxy'].astype(np.float32)
        
    cosine = (f.g['x']/r).in_units('1').astype(np.float32)
    sine = (f.g['y']/r).in_units('1').astype(np.float32)
    z = f.g['z']
    vel = f.g['vel']
    a = None # arbitrary initialization
    
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
    for iGrav in range(2):
        # Save files
        f.write(filename=f_name, fmt = pynbody.tipsy.TipsySnap)
        configsave(p_temp, p_name, ftype='param')
        
        if iGrav == 0:
            # Run ChaNGa calculating all forces (for initial run)
            command = ICgen_utils.changa_command(p_name, changa_preset, changbin, '+gas +n 0')
        else:
            # Run ChaNGa, only calculating gravity (on second run)
            command = ICgen_utils.changa_command(p_name, changa_preset, changbin, '-gas +n 0')
            
        print command
        p = ICgen_utils.changa_run(command)
        p.wait()
        
        if (iGrav == 0) and est_eps:
            # Estimate the gravitational softening length on the first iteration
            smoothlength_file = f_prefix + '.000000.smoothlength'
            eps = ICgen_utils.est_eps(smoothlength_file)
            f.g['eps'] = eps
    
        # Load accelerations
        acc_name = f_prefix + '.000000.acc2'
        del a
        gc.collect()
        a = load_acc(acc_name, low_mem=True)
        gc.collect()
        
        # Clean-up
        for fname in glob.glob(f_prefix + '*'): os.remove(fname)
        
        # Calculate cos(theta) where theta is angle above x-y plane
        cos = (r/np.sqrt(r**2 + z**2)).in_units('1').astype(np.float32)
        # Calculate radial acceleration times r^2
        ar2 = (a[:,0]*cosine + a[:,1]*sine)*r**2
        
        # Bin the data
        r_edges = np.linspace(r.min(), (1+np.spacing(2))*r.max(), nr + 1)
        ind, r_edges = digitize_threshold(r, min_per_bin, r_edges)
        ind -= 1
        nr = len(r_edges) - 1
        
        r_bins, ar2_mean, err = binned_mean(r, ar2, binedges=r_edges, \
        weighted_bins=True)
        
        gc.collect()
        
        # Fit lines to ar2 vs cos for each radial bin
        m = np.zeros(nr)
        b = np.zeros(nr)    
        
        for i in range(nr):
            
            mask = (ind == i)
            p = np.polyfit(cos[mask], ar2[mask], 1)
            m[i] = p[0]
            b[i] = p[1]
            
        # Interpolate the line fits
        m_spline = extrap1d(r_bins, m)
        b_spline = extrap1d(r_bins, b)
        
        # Calculate circular velocity
        ar2 = SimArray(m_spline(r)*cos + b_spline(r), ar2.units)
        gc.collect()
        v_calc = (np.sqrt(abs(ar2)/r)).in_units(vel.units)
        gc.collect()
        vel[:,0] = -v_calc*sine
        vel[:,1] = v_calc*cosine
        del v_calc
        gc.collect()
        
    # --------------------------------------------
    # Estimate pressure/gas dynamics accelerations
    # --------------------------------------------
    
    # Save files
    f.write(filename=f_name, fmt = pynbody.tipsy.TipsySnap)
    configsave(p_temp, p_name, ftype='param')
    
    # Run ChaNGa, including SPH
    command = ICgen_utils.changa_command(p_name, changa_preset, changbin, '+gas -n 0')
    p = ICgen_utils.changa_run(command)
    p.wait()
        
    # Load accelerations
    acc_name = f_prefix + '.000000.acc2'
    a_total = load_acc(acc_name, low_mem=True)
    gc.collect()
    
    # Estimate the accelerations due to pressure gradients/gas dynamics
    a_gas = a_total - a
    absa = np.sqrt((a_total**2).sum(1)) # magnitude of the acceleration
    del a_total, a
    gc.collect()
    ar2_gas = (a_gas[:,0]*cosine + a_gas[:,1]*sine)*r**2
    del a_gas
    gc.collect()
    
    logr_bins, ratio, err = binned_mean(np.log(r), ar2_gas/ar2, nbins=nr,\
    weighted_bins=True)
    r_bins = np.exp(logr_bins)
    del ar2_gas
    gc.collect()
    ratio_spline = extrap1d(r_bins, ratio)
    
    # Calculate time stepping parameters
    f0 = pynbody.load(f_prefix + '.000000')
    etaGrav = param.get('dEta', 0.2)
    dtGrav = etaGrav * np.sqrt(f0.g['eps']/absa)
    etaCourant = param.get('dEtaCourant', 0.4)
    mumax = f0.g['mumax']
    mumax[mumax < 0] = 0
    h = f0.g['smoothlength']
    dtCourant = etaCourant * h/(1.6*f0.g['c'] + 1.2*mumax)
    
    
    # If not all the particles were used for calculating velocity,
    # Make sure to use them now
    if subview:
        
        # Re-scale mass back to normal
        f.g['mass'] /= m_scale
        # Scale eps appropriately
        f.g['eps'] /= m_scale**(1.0/3)
        complete_snapshot.g['eps'] = f.g['eps'][[0]]
        # Re-scale time steps
        dtGrav /= m_scale**(1.0/6)
        dtCourant /= m_scale**(1.0/3)
        
        # Rename complete snapshot
        f = complete_snapshot
        # Calculate stuff for all particles
        r = f.g['rxy']
        z = f.g['z']
        cos = (r/np.sqrt(r**2 + z**2)).in_units('1').astype(np.float32)
        ar2 = SimArray(m_spline(r)*cos + b_spline(r), ar2.units)
        cosine = (f.g['x']/r).in_units('1').astype(np.float32)
        sine = (f.g['y']/r).in_units('1').astype(np.float32)
        vel = f.g['vel']
    
    dt = np.array([dtGrav, dtCourant]).min(0)
    dDelta = np.median(dt)    
    
    ar2_calc = ar2*(1 + ratio_spline(r))
    del ar2
    gc.collect()
    
    # Calculate velocity
    v = (np.sqrt(abs(ar2_calc)/r)).in_units(f.g['vel'].units)
    del ar2_calc
    gc.collect()
    
    vel[:,0] = -v*sine
    vel[:,1] = v*cosine
    
    # Clean-up
    for fname in glob.glob(f_prefix + '*'): os.remove(fname)
    
    return dDelta
