# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 15:11:31 2014

@author: ibackus
@editor: dflemin3
-Note: indentation is 4 spaces in this file, not a tab!

This module initializes an S-type binary system in which the gas disk is around
the primary, not both stars!  Assumes a_bin >> r_disk such that the disk's
velocity is dominated by the influence of the primary.
"""

__version__ = "$Revision: 1 $"
# $Source$

import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import gc
import AddBinary
import ICgen_utils

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
 
    #Load user supplied snapshot (assumed to be in cwd)
    path = "/astro/store/scratch/tmp/dflemin3/nbodyshare/9au-Q1.05-129K/"
    snapshot = pynbody.load(path + snapshotName)
 
    # particle positions
    r = snapshot.gas['r']
    xyz = snapshot.gas['pos']
    
    # Number of particles
    nParticles = len(snapshot.gas)
    
    # molecular mass
    m = settings.physical.m
    
    #Pull star mass from user-supplied snapshot
    ICobj.settings.physical.M = snapshot.star['mass'] #Total stellar mass in solar masses
    m_star = ICobj.settings.physical.M
    
    # disk mass
    m_disk = np.sum(snapshot.gas['mass'])
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
    
    # Estimate the star's softening length as the closest particle distance
    eps = r.min()
    
    # Make param file
    param = make_param(snapshot, snapshotName)
    param['dMeanMolWeight'] = m
       
    gc.collect()
    
    # CALCULATE VELOCITY USING calc_velocity.py.  This also estimates the 
    # gravitational softening length eps

    preset = settings.changa_run.preset
  
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
    
    """
    Now that the gas disk is initializes around the primary (M=m1), add in the
    second star as specified by the user.
    """    
    
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
    m_disk = strip_units(np.sum(snapshotBinary.gas['mass']))
    binsys.m1 = strip_units(m_star)
    binsys.m1 = binsys.m1 + m_disk  
    #Recompute cartesian coords considering primary as m1+m_disk    
    binsys.computeCartesian()
    
    x1,x2,v1,v2 = binsys.generateICs()

    #Assign position, velocity assuming CCW orbit
    snapshotBinary.star[0]['pos'] = SimArray(x1,pos_unit)
    snapshotBinary.star[0]['vel'] = SimArray(v1,v_unit)
    snapshotBinary.star[1]['pos'] = SimArray(x2,pos_unit)
    snapshotBinary.star[1]['vel'] = SimArray(v2,v_unit)

    """
    We have the binary positions about their center of mass, (0,0,0), so 
    shift the position, velocity of the gas disk to be around the primary.
    """
    snapshotBinary.gas['pos'] += snapshotBinary.star[0]['pos']
    snapshotBinary.gas['vel'] += snapshotBinary.star[0]['vel']  
    
    #Set stellar masses: Create simArray for mass, convert units to simulation mass units
    snapshotBinary.star[0]['mass'] = SimArray(binsys.m1-m_disk,m_unit)
    snapshotBinary.star[1]['mass'] = SimArray(binsys.m2,m_unit)
    snapshotBinary.star['metals'] = SimArray(star_metals)
 
    print 'Wrapping up'
    # Now set the star particle's tform to a negative number.  This allows
    # UW ChaNGa treat it as a sink particle.
    snapshotBinary.star['tform'] = -1.0
    
    #Set sink radius, stellar smoothing length as fraction of distance
    #from primary to inner edge of the disk
    r_sink = eps
    snapshotBinary.star[0]['eps'] = SimArray(r_sink/2.0,pos_unit)
    snapshotBinary.star[1]['eps'] = SimArray(r_sink/2.0,pos_unit)
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    param['dSinkMassMin'] = 0.9 * binsys.m2
    param['bDoSinks'] = 1
    
    return snapshotBinary, param, director
