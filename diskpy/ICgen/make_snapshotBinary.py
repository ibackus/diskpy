# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 15:11:31 2014

@author: ibackus
@editor: dfleming
-Note: indentation is 4 spaces in this file, not a tab!

"""

__version__ = "$Revision: 1 $"
# $Source$

import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import math
import gc
import AddBinary
import calc_velocity
import ICgen_utils
from diskpy import global_settings

from diskpy.utils import match_units, strip_units, configsave
from diskpy.pychanga import make_param, make_director

def snapshot_gen(ICobj):
    """
    Generates a tipsy snapshot from the initial conditions object ICobj.
    
    Returns snapshot, param
    
        snapshot: tipsy snapshot
        param: dictionary containing info for a .param file
    Note: Code has been edited (dflemin3) such that now it returns a snapshot for a circumbinary disk
    where initial conditions generated assuming star at origin of mass M.  After gas initialized, replaced
    star at origin with binary system who's center of mass lies at the origin and who's mass m1 +m2 = M
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
    
    # re-scale the particles (allows making of low-mass disk)
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
    # Make it a unit, save value for future conversion
    v_unit_vel = v_unit
    #Ensure v_unit_vel is the same as what I assume it is.
    assert(np.fabs(AddBinary.VEL_UNIT-v_unit_vel)<AddBinary.SMALL),"VEL_UNIT not equal to ChaNGa unit! Why??"			
	
    v_unit = pynbody.units.Unit('{0} km s**-1'.format(v_unit))
    
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
    #snapshot.star['eps'] = r.min()
    
    # Make param file
    param = make_param(snapshot, snapshotName)
    param['dMeanMolWeight'] = m
       
    gc.collect()
    
    # CALCULATE VELOCITY USING calc_velocity.py.  This also estimates the 
    # gravitational softening length eps
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

    #Now that velocities and everything are all initialized for gas particles, create new snapshot to return in which
    #single star particle is replaced by 2, same units as above
    snapshotBinary = pynbody.new(star=2,gas=nParticles)
    snapshotBinary['eps'] = 0.01*SimArray(np.ones(nParticles+2, dtype=np.float32), pos_unit)
    snapshotBinary['metals'] = SimArray(np.zeros(nParticles+2, dtype=np.float32))
    snapshotBinary['vel'].units = v_unit
    snapshotBinary['pos'].units = pos_unit
    snapshotBinary['mass'].units = snapshot['mass'].units
    snapshotBinary['rho'] = SimArray(np.zeros(nParticles+2, dtype=np.float32))

    #Assign gas particles with calculated/given values from above
    snapshotBinary.gas['pos'] = snapshot.gas['pos']
    snapshotBinary.gas['vel'] = snapshot.gas['vel']
    snapshotBinary.gas['temp'] = snapshot.gas['temp']
    snapshotBinary.gas['rho'] = snapshot.gas['rho']
    snapshotBinary.gas['eps'] = snapshot.gas['eps']
    snapshotBinary.gas['mass'] = snapshot.gas['mass']
    snapshotBinary.gas['metals'] = snapshot.gas['metals']

    #Load Binary system obj to initialize system
    binsys = ICobj.settings.physical.binsys
    
    x1,x2,v1,v2 = binsys.generateICs()

    #Put velocity in sim units
    #!!! Note: v_unit_vel will always be 29.785598165 km/s when m_unit = Msol and r_unit = 1 AU in kpc!!!
    #conv = v_unit_vel #km/s in sim units
    #v1 /= conv
    #v2 /= conv

    #Assign position, velocity assuming CCW orbit

    snapshotBinary.star[0]['pos'] = SimArray(x1,pos_unit)
    snapshotBinary.star[0]['vel'] = SimArray(v1,v_unit)
    snapshotBinary.star[1]['pos'] = SimArray(x2,pos_unit)
    snapshotBinary.star[1]['vel'] = SimArray(v2,v_unit)

    #Set stellar masses
    #Set Mass units
    #Create simArray for mass, convert units to simulation mass units
    priMass = SimArray(binsys.m1,m_unit)
    secMass = SimArray(binsys.m2,m_unit)

    snapshotBinary.star[0]['mass'] = priMass
    snapshotBinary.star[1]['mass'] = secMass
    snapshotBinary.star['metals'] = SimArray(star_metals)

    #Estimate stars' softening length as fraction of distance to COM
    d = np.sqrt(AddBinary.dotProduct(x1-x2,x1-x2))

    snapshotBinary.star[0]['eps'] = SimArray(math.fabs(d)/4.0,pos_unit)
    snapshotBinary.star[1]['eps'] = SimArray(math.fabs(d)/4.0,pos_unit)
 
    print 'Wrapping up'
    # Now set the star particle's tform to a negative number.  This allows
    # UW ChaNGa treat it as a sink particle.
    snapshotBinary.star['tform'] = -1.0
    
    #Set Sink Radius to be mass-weighted average of Roche lobes of two stars
    r1 = AddBinary.calcRocheLobe(binsys.m1/binsys.m2,binsys.a) 
    r2 = AddBinary.calcRocheLobe(binsys.m2/binsys.m1,binsys.a)
    p = strip_units(binsys.m1/(binsys.m1 + binsys.m2))

    r_sink = (r1*p) + (r2*(1.0-p))
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    param['dSinkMassMin'] = 0.9 * strip_units(secMass)
    param['bDoSinks'] = 1
    
    return snapshotBinary, param, director
    
        
