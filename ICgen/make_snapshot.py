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

# diskpy packages
from diskpy.utils import match_units, configsave, strip_units
from diskpy.pychanga import make_director, make_param
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
    # filenames
    snapshotName = settings.filenames.snapshotName
    paramName = settings.filenames.paramName
        
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
    m_disk = match_units(m_disk, m_star)[0]
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
    v_unit = pynbody.units.Unit('{0} km s**-1'.format(v_unit))
    
    # Other settings
    metals = settings.snapshot.metals
    star_metals = metals
    
    # -------------------------------------------------
    # Initialize snapshot
    # -------------------------------------------------
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
    param = make_param(snapshot, snapshotName)
    param['dMeanMolWeight'] = m
    eos = (settings.physical.eos).lower()
    
    if eos == 'adiabatic':
        
        param['bGasAdiabatic'] = 1
        param['bGasIsothermal'] = 0
        
    param['dConstGamma']
       
    gc.collect()
    
    # -------------------------------------------------
    # CALCULATE VELOCITY USING calc_velocity.py.  This also estimates the 
    # gravitational softening length eps
    # -------------------------------------------------
    print 'Calculating circular velocity'
    preset = settings.changa_run.preset
    max_particles = global_settings['misc']['max_particles']
    calc_velocity.v_xy(snapshot, param, changa_preset=preset, max_particles=max_particles)
    
    gc.collect()
    
    # -------------------------------------------------
    # Estimate time step for changa to use
    # -------------------------------------------------
    # Save param file
    configsave(param, paramName, 'param')
    # Save snapshot
    snapshot.write(filename=snapshotName, fmt=pynbody.tipsy.TipsySnap)
    # est dDelta
    dDelta = ICgen_utils.est_time_step(paramName, preset)
    param['dDelta'] = dDelta
    
    # -------------------------------------------------
    # Create director file
    # -------------------------------------------------
    # largest radius to plot
    r_director = float(0.9 * r.max())
    # Maximum surface density
    sigma_min = float(ICobj.sigma(r_director))
    # surface density at largest radius
    sigma_max = float(ICobj.sigma.input_dict['sigma'].max())
    # Create director dict
    director = make_director(sigma_min, sigma_max, r_director, filename=param['achOutName'])
    ## Save .director file
    #configsave(director, directorName, 'director')
    
    # -------------------------------------------------
    # Wrap up
    # -------------------------------------------------
    print 'Wrapping up'
    # Now set the star particle's tform to a negative number.  This allows
    # UW ChaNGa treat it as a sink particle.
    snapshot.star['tform'] = -1.0
    
    # Update params
    r_sink = strip_units(r.min())
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    param['dSinkMassMin'] = 0.9 * strip_units(m_star)
    param['bDoSinks'] = 1
    
    return snapshot, param, director
