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
import os
import AddBinary
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
    m_disk = isaac.match_units(m_disk, m_star)[0]
    
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
<<<<<<< HEAD
    #Ensure v_unit_vel is the same as what I assume it is.
    assert(np.fabs(AddBinary.VEL_UNIT-v_unit_vel)<AddBinary.SMALL),"VEL_UNIT not equal to ChaNGa unit! Why??"			
	
=======
>>>>>>> master
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
  
	# -------------------------------------------------
    # Estimate time step for changa to use
    # -------------------------------------------------
    # Save param file
    isaac.configsave(param, paramName, 'param')
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
    director = isaac.make_director(sigma_min, sigma_max, r_director, filename=param['achOutName'])
    ## Save .director file
    #isaac.configsave(director, directorName, 'director')

<<<<<<< HEAD
    #Now that velocities and everything are all initialized for gas particles, create new snapshot to return in which
=======
    #Now that velocities are all initialized for gas particles, create new snapshot to return in which
>>>>>>> master
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

<<<<<<< HEAD
    #Load Binary system obj to initialize system
    binsys = ICobj.settings.physical.binsys
    
    x1,x2,v1,v2 = binsys.generateICs()

    #Put velocity in sim units
    #!!! Note: v_unit_vel will always be 29.785598165 km/s when m_unit = Msol and r_unit = 1 AU in kpc!!!
    #conv = v_unit_vel #km/s in sim units
    #v1 /= conv
    #v2 /= conv

    #Assign position, velocity assuming CCW orbit
=======
    #Calculate binary star parameters
    #***Consider binary system where m_star = m1 + m2***
    #1 = primary (more massive) star, 2 = secondary

    #Calculate Binary Star masses
    p = ICobj.settings.physical.priMassPerc
    m1 = isaac.strip_units(m_star)*p
    m2 = isaac.strip_units(m_star)*(1.0-p)

    #Calculate Binary Semimajor axis (in au) given period
    period = ICobj.settings.physical.period #days    
    e = ICobj.settings.physical.ecc #eccentricity
    a = AddBinary.pToA(period,m1+m2)
    
    #Get other Orbital Parameters to use later
    w = ICobj.settings.physical.w #Argument of Pericenter
    Omega = ICobj.settings.physical.Omega #Longitude of Ascending Node
    inc = ICobj.settings.physical.inc
    MA = ICobj.settings.physical.M

    x1, x2, v1, v2 = AddBinary.initializeBinary(a,e,inc,Omega,w,MA,m1,m2,angleFlag=True,scaleFlag=False)

    #Calculate Binary Star x positions, perihelion velocities
    #x1, x2 = AddBinary.calcPositions(m1+m2,a,e,p)
    #v1, v2 = AddBinary.calcV(m1,m2,a,e)

    #Put velocity in sim units
    #!!! Note: v_unit_vel will always be 29.785598165 km/s when m_unit = Msol and r_unit = 1 AU in kpc!!!
    conv = v_unit_vel #km/s in sim units
    #v1 /= conv
    #v2 /= conv

    #Assign stellar parameters for CCW rotation
    #snapshotBinary.star[0]['pos'] = SimArray([[x1, 0., 0.]],pos_unit)
    #snapshotBinary.star[0]['vel'] = SimArray([[0., v1, 0.]],v_unit)
    
>>>>>>> master
    snapshotBinary.star[0]['pos'] = SimArray(x1,pos_unit)
    snapshotBinary.star[0]['vel'] = SimArray(v1,v_unit)
    snapshotBinary.star[1]['pos'] = SimArray(x2,pos_unit)
    snapshotBinary.star[1]['vel'] = SimArray(v2,v_unit)

<<<<<<< HEAD
    #Set stellar masses
    #Set Mass units
    #Create simArray for mass, convert units to simulation mass units
    priMass = SimArray(binsys.m1,m_unit)
    secMass = SimArray(binsys.m2,m_unit)
=======
    #snapshotBinary.star[1]['pos'] = SimArray([[x2, 0., 0.]],pos_unit)
    #snapshotBinary.star[1]['vel'] = SimArray([[0., v2, 0.]],v_unit)

    #Set stellar masses
    #Set Mass units
    #Create simArray for mass, convert units to simulation mass units
    priMass = SimArray(m1,m_unit)
    secMass = SimArray(m2,m_unit)
>>>>>>> master
    snapshotBinary.star[0]['mass'] = priMass
    snapshotBinary.star[1]['mass'] = secMass
    snapshotBinary.star['metals'] = SimArray(star_metals)

<<<<<<< HEAD
    #Estimate stars' softening length as fraction of distance to COM
    d = np.linalg.norm(x1) - np.linalg.norm(x2)
=======
    d = np.linalg.norm(x1) - np.linalg.norm(x2)
    #Estimate stars' softening length as fraction of distance to COM
>>>>>>> master
    snapshotBinary.star[0]['eps'] = SimArray(math.fabs(d)/4.0,pos_unit)
    snapshotBinary.star[1]['eps'] = SimArray(math.fabs(d)/4.0,pos_unit)
 
    print 'Wrapping up'
    # Now set the star particle's tform to a negative number.  This allows
    # UW ChaNGa treat it as a sink particle.
    snapshotBinary.star['tform'] = -1.0
    
<<<<<<< HEAD
    #Set Sink Radius to be mass-weighted average of Roche lobes of two stars
    r1 = AddBinary.calcRocheLobe(binsys.m1/binsys.m2,binsys.a) 
    r2 = AddBinary.calcRocheLobe(binsys.m2/binsys.m1,binsys.a)
    p = isaac.strip_units(binsys.m1/(binsys.m1 + binsys.m2))
=======
    #Set Sink Radius to be weighted average of Roche lobes of two stars
    r1 = AddBinary.calcRocheLobe(m1/m2,a) 
    r2 = AddBinary.calcRocheLobe(m2/m1,a)
>>>>>>> master
    r_sink = (r1*p) + (r2*(1.0-p))
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    param['dSinkMassMin'] = 0.9 * isaac.strip_units(secMass)
    param['bDoSinks'] = 1
    
    return snapshotBinary, param, director
    
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
        
