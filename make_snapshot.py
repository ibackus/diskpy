# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 15:11:31 2014

@author: ibackus
"""

__version__ = "$Revision: 1 $"
# $Source$

import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import gc
import os

import isaac
import calc_velocity
import ICgen_utils
import ICglobal_settings
global_settings = ICglobal_settings.global_settings

def snapshot_gen(ICobj):
    """
    Generates a tipsy snapshot from the initial conditions object ICobj.
    
    Returns snapshot, param
    
        snapshot: tipsy snapshot
        param: dictionary containing info for a .param file
    """
    
    print 'Generating snapshot...'
    # Constants
    G = SimArray(1.0,'G')
    # ------------------------------------
    # Load in things from ICobj
    # ------------------------------------
    print 'Accessing data from ICs'
    settings = ICobj.settings
    # snapshot file name
    snapshotName = settings.filenames.snapshotName
    # particle positions
    r = ICobj.pos.r
    xyz = ICobj.pos.xyz
    # Number of particles
    nParticles = ICobj.pos.nParticles
    # molecular mass
    m = settings.physical.m
    # star mass
    m_star = settings.physical.M.copy()
    # disk mass
    m_disk = ICobj.sigma.m_disk.copy()
    m_disk = isaac.match_units(m_disk, m_star)[0]
    # mass of the gas particles
    m_particles = m_disk / float(nParticles)
    # re-scale the particles (allows making of lo-mass disk)
    m_particles *= settings.snapshot.mScale
    
    # -------------------------------------------------
    # Assign output
    # -------------------------------------------------
    print 'Assigning data to snapshot'
    # Get units all set up
    m_unit = m_star.units
    pos_unit = r.units
    
    if xyz.units != r.units:
        
        xyz.convert_units(pos_unit)
        
    # time units are sqrt(L^3/GM)
    t_unit = np.sqrt((pos_unit**3)*np.power((G*m_unit), -1)).units
    # velocity units are L/t
    v_unit = (pos_unit/t_unit).ratio('km s**-1')
    # Make it a unit
    v_unit = pynbody.units.Unit('{} km s**-1'.format(v_unit))
    
    # Other settings
    metals = settings.snapshot.metals
    star_metals = metals
    
    # Generate snapshot
    # Note that empty pos, vel, and mass arrays are created in the snapshot
    snapshot = pynbody.new(star=1,gas=nParticles)
    snapshot['vel'].units = v_unit
    snapshot['eps'] = 0.01*SimArray(np.ones(nParticles+1, dtype=np.float32), pos_unit)
    snapshot['metals'] = SimArray(np.zeros(nParticles+1, dtype=np.float32))
    snapshot['rho'] = SimArray(np.zeros(nParticles+1, dtype=np.float32))
    
    snapshot.gas['pos'] = xyz
    snapshot.gas['temp'] = ICobj.T(r)
    snapshot.gas['mass'] = m_particles
    snapshot.gas['metals'] = metals
    
    snapshot.star['pos'] = SimArray([[ 0.,  0.,  0.]],pos_unit)
    snapshot.star['vel'] = SimArray([[ 0.,  0.,  0.]], v_unit)
    snapshot.star['mass'] = m_star
    snapshot.star['metals'] = SimArray(star_metals)
    # Estimate the star's softening length as the closest particle distance
    snapshot.star['eps'] = r.min()
    
    # Make param file
    param = isaac.make_param(snapshot, snapshotName)
    param['dMeanMolWeight'] = m
       
    gc.collect()
    
    # CALCULATE VELOCITY USING calc_velocity.py.  This also estimates the 
    # gravitational softening length eps
    print 'Calculating circular velocity'
    preset = settings.changa_run.preset
    max_particles = global_settings['misc']['max_particles']
    calc_velocity.v_xy(snapshot, param, changa_preset=preset, max_particles=max_particles)
    
    gc.collect()
    
    print 'Wrapping up'
    # Now set the star particle's tform to a negative number.  This allows
    # UW ChaNGa treat it as a sink particle.
    snapshot.star['tform'] = -1.0
    
    # Update params
    r_sink = isaac.strip_units(r.min())
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    param['dSinkMassMin'] = 0.9 * isaac.strip_units(m_star)
    param['bDoSinks'] = 1
    
    return snapshot, param
    
def make_director(ICobj, res=1200):
    
    director = {}
    director['render'] = 'tsc'
    director['FOV'] = 45.0
    director['clip'] = [0.0001, 500]
    director['up'] = [1, 0, 0]
    director['project'] = 'ortho'
    director['softgassph'] = 'softgassph'
    director['physical'] = 'physical'
    director['size'] = [res, res]
    
    sig_set = ICobj.settings.sigma
    mScale = ICobj.settings.snapshot.mScale
    snapshot_name = ICobj.settings.filenames.snapshotName
    f_prefix = os.path.splitext(os.path.basename(snapshot_name))[0]
    
    director['file'] = f_prefix
    
    
    if sig_set.kind == 'MQWS':
        
        rmax = sig_set.rout + 3*sig_set.rin
        zmax = float(rmax)
        director['eye'] = [0, 0, zmax]
        vmin = float(ICobj.rho(0, rmax))
        vmax = float(ICobj.rho.rho_binned[0,:].max())
        vmax *= mScale
        director['logscale'] = [vmin, 10*vmax]
        director['colgas'] = [1, 1, 1]
        
    return director
        