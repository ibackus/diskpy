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
    kB = SimArray(1.0,'k')
    # ------------------------------------
    # Load in things from ICobj
    # ------------------------------------
    print 'Accessing data from ICs'
    settings = ICobj.settings
    # snapshot file name
    snapshotName = settings.filenames.snapshotName
    # particle positions
    theta = ICobj.pos.theta
    r = ICobj.pos.r
    x = ICobj.pos.x
    y = ICobj.pos.y
    z = ICobj.pos.z
    # Number of particles
    nParticles = ICobj.pos.nParticles
    dr = (ICobj.sigma.r_bins[[1]] - ICobj.sigma.r_bins[[0]])/10.0
    # molecular mass
    m = settings.physical.m
    # star mass
    m_star = settings.physical.M.copy()
    # disk mass
    m_disk = ICobj.sigma.m_disk.copy()
    m_disk = isaac.match_units(m_disk, m_star)[0]
    # mass of the gas particles
    m_particles = np.ones(nParticles) * m_disk / float(nParticles)
    # re-scale the particles (allows making of lo-mass disk)
    m_particles *= settings.snapshot.mScale
    
    # ------------------------------------
    # Initial calculations
    # ------------------------------------
    print 'Running initial calculations'
    # Find total mass interior to every particle
#    N_interior = np.array(r.argsort().argsort())
#    m_int = m_particles[[0]]*N_interior + m_star
#    # Retrieve rho (density) at each position
#    rho = ICobj.rho(z,r)
#    # Retrieve radial derivative at each position
#    drho_dr = ICobj.rho.drho_dr(z,r)
    # Get temperature at each position
    T = ICobj.T(r)
    
#    # ------------------------------------
#    # Calculate particle velocities
#    # ------------------------------------
#    print 'Calculating initial guess for particle velocities'
#    # Find keperlerian velocity squared due to gravity
#    v2grav = G*m_int/r
#    # Find contribution from density gradient
#    v2dens = (kB*T/m)*(r*drho_dr/rho)
#    #       ignore nans and infs
#    v2dens[(np.isnan(v2dens)) | (np.isinf(v2dens))] = 0.0
#    # Find contribution from temperature gradient
#    
#    dT_dr = (ICobj.T(r+dr) - ICobj.T(r-dr))/(2*dr)
#    v2temp = r * dT_dr * kB/m
#    #v2temp = (kB*T/m)*Tpower
#    # Now find velocity from all contributions
#    v = np.sqrt(v2grav + v2dens + v2temp)
#    # Sometimes, at large r, the velocities due to the pressure and temp
#    # Gradients become negative.  If this is the case, set them to 0
#    # Also, the temperature gradient can become infinite at r=0
#    nanind = np.isnan(v) | np.isinf(v)
#    v[nanind] = 0.0
    
    # -------------------------------------------------
    # Assign output
    # -------------------------------------------------
    print 'Assigning data to snapshot'
    # Get units all set up
    m_unit = m_star.units
    pos_unit = r.units
    # time units are sqrt(L^3/GM)
    t_unit = np.sqrt((pos_unit**3)*np.power((G*m_unit), -1)).units
    # velocity units are L/t
    v_unit = (pos_unit/t_unit).ratio('km s**-1')
    # Make it a unit
    v_unit = pynbody.units.Unit('{} km s**-1'.format(v_unit))
    x.convert_units(pos_unit)
    y.convert_units(pos_unit)
    z.convert_units(pos_unit)
    
    # 3-D velocity
    vel = SimArray(np.zeros([nParticles,3]),v_unit)
#    vel[:,0] = -np.sin(theta)*v
#    vel[:,1] = np.cos(theta)*v
    
    # Generate positions
    xyz = SimArray(np.zeros([nParticles,3]),pos_unit)
    xyz[:,0] = x
    xyz[:,1] = y
    xyz[:,2] = z
    
    # Other settings
    metals = settings.snapshot.metals
    star_metals = metals
    metals *= SimArray(np.ones(nParticles))
    
    # Generate snapshot
    snapshot = pynbody.new(star=1,gas=nParticles)
    snapshot.gas['vel'] = vel
    snapshot.gas['pos'] = xyz
    snapshot.gas['temp'] = T
    snapshot.gas['mass'] = m_particles
    snapshot.gas['metals'] = metals
    # Initial eps...totally arbitrary (it gets estimated below)
    snapshot.gas['eps'] = 0.01
    snapshot.gas['rho'] = 0
    
    snapshot.star['pos'] = SimArray([[ 0.,  0.,  0.]],pos_unit)
    snapshot.star['vel'] = SimArray([[ 0.,  0.,  0.]], v_unit)
    snapshot.star['mass'] = m_star
    snapshot.star['metals'] = SimArray(star_metals)
    # Estimate the star's softening length as the closest particle distance
    snapshot.star['eps'] = r.min()
    snapshot.star['rho'] = 0
    
    # Make param file
    param = isaac.make_param(snapshot, snapshotName)
    
    # Estimate reasonable gravitational softening for the gas
    print 'Estimating gas gravitational softening'
    preset = settings.changa_run.preset
    snapshot.g['eps'] = ICgen_utils.est_eps(snapshot, preset)
    
    # CALCULATE VELOCITY USING calc_velocity.py
    print 'Calculating circular velocity'
    vel = calc_velocity.v_xy(snapshot, param, changa_preset=preset)
    snapshot.gas['vel'] = vel
    
    print 'Wrapping up'
    # Now set the star particle's tform to a negative number.  This allows
    # UW ChaNGa treat it as a sink particle.
    snapshot.star['tform'] = -1.0
    
    # Update params
    r_sink = isaac.strip_units(snapshot.g['rxy'].min())
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    param['dSinkMassMin'] = 0.9 * isaac.strip_units(m_star)
    param['bDoSinks'] = 1
    
    return snapshot, param