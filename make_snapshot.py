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

import isaac
import calc_velocity
import ICgen_utils

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
    
#    # ------------------------------------
#    # Initial calculations
#    # ------------------------------------
#    print 'Calculating temperature'
#    T = ICobj.T(r)
    
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
    
#    # 3-D velocity
#    vel = SimArray(np.zeros([nParticles,3]),v_unit)
    
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
    
#    snapshot.gas['vel'] = vel
    snapshot.gas['pos'] = xyz
    snapshot.gas['temp'] = ICobj.T(r)
    snapshot.gas['mass'] = m_particles
    snapshot.gas['metals'] = metals
    # Initial eps...totally arbitrary (it gets estimated below)
    #snapshot.gas['eps'] = 0.01
    #snapshot.gas['rho'] = 0
    
    snapshot.star['pos'] = SimArray([[ 0.,  0.,  0.]],pos_unit)
    snapshot.star['vel'] = SimArray([[ 0.,  0.,  0.]], v_unit)
    snapshot.star['mass'] = m_star
    snapshot.star['metals'] = SimArray(star_metals)
    # Estimate the star's softening length as the closest particle distance
    snapshot.star['eps'] = r.min()
    #snapshot.star['rho'] = 0
    
    # Make param file
    param = isaac.make_param(snapshot, snapshotName)
    param['dMeanMolWeight'] = m
    
    # Estimate reasonable gravitational softening for the gas
    preset = settings.changa_run.preset
    
    gc.collect()
    
    # CALCULATE VELOCITY USING calc_velocity.py
    print 'Calculating circular velocity'
    calc_velocity.v_xy(snapshot, param, changa_preset=preset,r=r,calc_eps=True)
    
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