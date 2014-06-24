# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 15:11:31 2014

@author: ibackus
"""

import pynbody
SimArray = pynbody.array.SimArray
import numpy as np

import isaac
import calc_velocity

def snapshot_gen(ICobj):
    """
    Generates a tipsy snapshot from the initial conditions object ICobj.
    
    Returns snapshot, param
    
        snapshot: tipsy snapshot
        param: dictionary containing info for a .param file
    """
    
    # Constants
    G = SimArray(1.0,'G')
    kB = SimArray(1.0,'k')
    # ------------------------------------
    # Load in things from ICobj
    # ------------------------------------
    # snapshot file name
    snapshotName = ICobj.settings.filenames.snapshotName
    # particle positions
    theta = ICobj.pos.theta
    r = ICobj.pos.r
    x = ICobj.pos.x
    y = ICobj.pos.y
    z = ICobj.pos.z
    # Number of particles
    nParticles = ICobj.pos.nParticles
    # Temperature power law (used for pressure gradient)
    Tpower = ICobj.settings.physical.Tpower
    dr = (ICobj.sigma.r_bins[[1]] - ICobj.sigma.r_bins[[0]])/10.0
    # molecular mass
    m = ICobj.settings.physical.m
    # star mass
    m_star = ICobj.settings.physical.M.copy()
    # disk mass
    m_disk = ICobj.sigma.m_disk.copy()
    m_disk = isaac.match_units(m_disk, m_star)[0]
    # mass of the gas particles
    m_particles = np.ones(nParticles) * m_disk / float(nParticles)
    # re-scale the particles (allows making of lo-mass disk)
    m_particles *= ICobj.settings.snapshot.mScale
    
    # ------------------------------------
    # Initial calculations
    # ------------------------------------
    # Find total mass interior to every particle
    N_interior = np.array(r.argsort().argsort())
    m_int = m_particles[[0]]*N_interior + m_star
    # Retrieve rho (density) at each position
    rho = ICobj.rho(z,r)
    # Retrieve radial derivative at each position
    drho_dr = ICobj.rho.drho_dr(z,r)
    # Get temperature at each position
    T = ICobj.T(r)
    
    # ------------------------------------
    # Calculate particle velocities
    # ------------------------------------
    # Find keperlerian velocity squared due to gravity
    v2grav = G*m_int/r
    # Find contribution from density gradient
    v2dens = (kB*T/m)*(r*drho_dr/rho)
    #       ignore nans and infs
    v2dens[(np.isnan(v2dens)) | (np.isinf(v2dens))] = 0.0
    # Find contribution from temperature gradient
    
    dT_dr = (ICobj.T(r+dr) - ICobj.T(r-dr))/(2*dr)
    v2temp = r * dT_dr * kB/m
    #v2temp = (kB*T/m)*Tpower
    # Now find velocity from all contributions
    v = np.sqrt(v2grav + v2dens + v2temp)
    # Sometimes, at large r, the velocities due to the pressure and temp
    # Gradients become negative.  If this is the case, set them to 0
    # Also, the temperature gradient can become infinite at r=0
    nanind = np.isnan(v) | np.isinf(v)
    v[nanind] = 0.0
    
    # -------------------------------------------------
    # Assign output
    # -------------------------------------------------
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
    vel[:,0] = -np.sin(theta)*v
    vel[:,1] = np.cos(theta)*v
    
    # Generate positions
    xyz = SimArray(np.zeros([nParticles,3]),pos_unit)
    xyz[:,0] = x
    xyz[:,1] = y
    xyz[:,2] = z
    
    # Other settings
    eps = ICobj.settings.snapshot.eps
    star_eps = eps
    eps *= SimArray(np.ones(nParticles), pos_unit)
    metals = ICobj.settings.snapshot.metals
    star_metals = metals
    metals *= SimArray(np.ones(nParticles))
    
    # Generate snapshot
    snapshot = pynbody.new(star=1,gas=nParticles)
    snapshot.gas['vel'] = vel
    snapshot.gas['pos'] = xyz
    snapshot.gas['temp'] = T
    snapshot.gas['mass'] = m_particles
    snapshot.gas['metals'] = metals
    snapshot.gas['eps'] = eps
    snapshot.gas['mu'].derived = False
    snapshot.gas['mu'] = float(m.in_units('m_p'))
    snapshot.gas['rho'] = 0
    
    snapshot.star['pos'] = SimArray([[ 0.,  0.,  0.]],pos_unit)
    snapshot.star['vel'] = SimArray([[ 0.,  0.,  0.]], v_unit)
    snapshot.star['mass'] = m_star
    snapshot.star['metals'] = SimArray(star_metals)
    snapshot.star['eps'] = SimArray(star_eps, pos_unit)
    snapshot.star['rho'] = 0
    
    # Make param file
    param = isaac.make_param(snapshot, snapshotName)
    
    # CALCULATE VELOCITY USING calc_velocity.py
    vel = calc_velocity.v_xy(snapshot, param)
    snapshot.gas['vel'] = vel
    
    # Now add sinks to the params.  This has to be done after calculating
    # velocity
    r_sink = isaac.strip_units(snapshot.g['rxy'].min())
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    
    # Now set the star particle's tform to a negative number.  This allows
    # UW ChaNGa treat it as a sink particle.
    snapshot.star['tform'] = -1.0
    
    return snapshot, param