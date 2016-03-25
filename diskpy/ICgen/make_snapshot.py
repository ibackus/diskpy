# -*- coding: utf-8 -*-
"""
Contains functions for genrerating a simulation snapshot from an IC object 
which already has particle positions generated.  Snapshot are made with
snapshot_gen(IC)

Created on Fri Mar 21 15:11:31 2014

@author: ibackus
"""

__version__ = "$Revision: 2 $"
# $Source$

import pynbody
SimArray = pynbody.array.SimArray
import gc
import numpy as np

# diskpy packages
import diskpy
from diskpy import global_settings
from diskpy.utils import match_units, strip_units
from diskpy.pychanga import make_director, make_param, setup_units
import calc_velocity
import AddBinary

# Constants
G = SimArray(1.0,'G')

def snapshot_gen(IC):
    """
    Generates a tipsy snapshot from the initial conditions object IC.  Includes
    a routine to calculate the velocity.
    
    Parameters
    ----------
    IC : ICobj
    
    Returns
    -------
    snapshot : SimSnap
        Simulation snapshot with the velocity calculated
    param : dict
        dictionary containing info for a .param file
    """
    
    print 'Generating snapshot...'   
    # Initialize snapshot
    snapshot = init_snapshot(IC)
    # Make param file
    param = init_param(IC, snapshot)
       
    # -------------------------------------------------
    # CALCULATE VELOCITY USING calc_velocity.py.  This also estimates the 
    # gravitational softening length eps and a good timestep
    # -------------------------------------------------
    print 'Calculating circular velocity'
    preset = IC.settings.changa_run.preset
    max_particles = global_settings['misc']['max_particles']
    dDelta = calc_velocity.v_xy(snapshot, param, \
    changa_preset=preset, max_particles=max_particles)
    param['dDelta'] = dDelta    
    print 'Calculated time step.  dDelta = ', dDelta
    gc.collect()
    
    # Create director file
    director = init_director(IC, param)
    # Handle binaries
    starMode = IC.settings.physical.starMode.lower()
    if starMode == 'binary':
        
        snapshot = make_binary(IC, snapshot)
        
    # Finalize
    print 'Wrapping up'
    setup_sinks(IC, snapshot, param)
    
    return snapshot, param, director
    
    
def make_binary(IC, snapshot):
    """
    Turns a snapshot for a single star into a snapshot of a binary system
    
    Parameters
    ----------
    IC : IC object
    snapshot : SimSnap
        Single star system to turn into a binary
    
    Returns
    -------
    snapshotBinary : SimSnap
        A binary version of the simulation snapshot
    """
    # Initialize snapshot
    snapshotBinary = pynbody.new(star=2, gas=len(snapshot.g))
    # Copy gas particles over
    for key in snapshot.gas.keys():
        
        snapshotBinary.gas[key] = snapshot.gas[key]
        
    # Load Binary system obj to initialize system
    starMode = IC.settings.physical.starMode.lower()
    binsys = IC.settings.physical.binsys
    
    if starMode == 'stype':
        
        # Treate the primary as a star of mass mStar + mDisk
        m_disk = strip_units(np.sum(snapshotBinary.gas['mass']))
        binsys.m1 += m_disk
        binsys.computeCartesian()
        
    x1,x2,v1,v2 = binsys.generateICs()
    
    #Assign star parameters assuming CCW orbit
    snapshotBinary.star[0]['pos'] = x1
    snapshotBinary.star[0]['vel'] = v1
    snapshotBinary.star[1]['pos'] = x2
    snapshotBinary.star[1]['vel'] = v2
    
    #Set stellar masses
    priMass = binsys.m1
    secMass = binsys.m2    
    snapshotBinary.star[0]['mass'] = priMass
    snapshotBinary.star[1]['mass'] = secMass    
    snapshotBinary.star['metals'] = snapshot.s['metals']
    
    if starMode == 'stype':
        # We have the binary positions about their center of mass, (0,0,0), so 
        # shift the position, velocity of the gas disk to be around the primary.
        snapshotBinary.gas['pos'] += snapshotBinary.star[0]['pos']
        snapshotBinary.gas['vel'] += snapshotBinary.star[0]['vel']  
        # Remove disk mass from the effective star mass
        snapshotBinary[0]['mass'] -= m_disk
        binsys.m1 -= m_disk
        # Star smoothing
        snapshotBinary.star['eps'] = snapshot.star['eps']

        
    elif starMode == 'binary':
        
        # Estimate stars' softening length as fraction of distance to COM
        d = np.sqrt(AddBinary.dotProduct(x1-x2,x1-x2))
        pos_unit = snapshotBinary['pos'].units
        snapshotBinary.star['eps'] = SimArray(abs(d)/4.0,pos_unit)
        
    return snapshotBinary

def init_director(IC, param=None):
    """
    Creates a director dict for a PPD.
    
    Parameters
    ----------
    IC : ICobj
    param : .param dict
        If not supplied, param = IC.param
    """    
    if param is None:
        
        param = IC.snapshot_param
        
    # largest radius to plot
    r_director = float(0.9 * IC.pos.r.max())
    # Maximum surface density
    sigma_min = float(IC.sigma(r_director))
    # surface density at largest radius
    sigma_max = float(IC.sigma.input_dict['sigma'].max())
    # Create director dict
    director = make_director(sigma_min, sigma_max, r_director, \
    filename=param['achOutName'])
    
    return director
    
def sink_radius(IC):
    """
    Determine a reasonable sink radius for the star particles depending on 
    the star system type (e.g., single star, binary, etc...)
    
    Parameters
    ----------
    IC : IC object
    
    Returns
    -------
    r_sink : SimArray
        Sink radius for star particles
    """
    
    # Set up the sink radius
    starMode = IC.settings.physical.starMode.lower()
    
    if starMode == 'binary':
        
        binsys = IC.settings.physical.binsys
        #Set Sink Radius to be mass-weighted average of Roche lobes of two stars
        r1 = AddBinary.calcRocheLobe(binsys.m1/binsys.m2,binsys.a) 
        r2 = AddBinary.calcRocheLobe(binsys.m2/binsys.m1,binsys.a)
        p = strip_units(binsys.m1/(binsys.m1 + binsys.m2))
        r_sink = (r1*p) + (r2*(1.0-p))
    
    else:
    
        r_sink = IC.pos.r.min()
        
    return r_sink
    
def setup_sinks(IC, snapshot, param):
    """
    Sets up snapshot and param for stars that are sinks.
    
    Parameters
    ----------
    IC : IC obj
    snapshot : SimSnap
    param : param dict
    
    Returns
    -------
    None
    """
    units = diskpy.pychanga.units_from_param(param)
    # Set the star tforms to a negative number.  This allows UW ChaNGa treat 
    # stars as sink particles
    snapshot.star['tform'] = -1.0   
    # Set sink radius for stars
    r_sink = sink_radius(IC)
    r_sink = float(strip_units(r_sink))
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    # Set sink mass to be 90% of the smallest star
    Mstar = snapshot.s['mass'].min()
    Mstar.convert_units(units['m_unit'])
    param['dSinkMassMin'] = 0.9 * strip_units(Mstar)
    # Turn sinks on
    param['bDoSinks'] = 1

def init_param(IC, snapshot=None):
    """
    Initializes a ChaNGa param dict (see also diskpy.pychanga.make_param) for
    an IC object.
    
    Parameters
    ----------
    IC : ICobject
        Initial conditions object containing required settings for generating
        the param dict
    snapshot : SimSnap
        Snapshot to create param for.  If None, IC.snapshot is used
        
    Returns
    -------
    param : dict
        Dict containing the param.  Can be saved with diskpy.utils.configsave
    """
    if snapshot is None:
        
        snapshot = IC.snapshot
        
    # Make param file
    settings = IC.settings
    snapshotName = settings.filenames.snapshotName
    param = make_param(snapshot, snapshotName)
    param['dMeanMolWeight'] = settings.physical.m
    eos = (settings.physical.eos).lower()
    
    if eos == 'adiabatic':
        
        param['bGasAdiabatic'] = 1
        param['bGasIsothermal'] = 0
        
    return param

def init_snapshot(IC):
    """
    Initialize a snapshot for the IC object.  Requires that positions have
    been created.  Also sets:
     * pos
     * metals
     * temp
     * mass
     * star eps
    
    Parameters
    ----------
    IC : ICobj
    
    Returns
    -------
    snapshot : SimSnap
    """
    # Get required settings from IC
    settings = IC.settings
    # particle positions
    r = IC.pos.r
    xyz = IC.pos.xyz
    nParticles = IC.pos.nParticles
    m_star = settings.physical.M
    m_disk = IC.sigma.m_disk
    m_disk = match_units(m_disk, m_star)[0]
    m_particles = m_disk / float(nParticles)
    metals = settings.snapshot.metals
    # re-scale the particles (allows making of lo-mass disk)
    m_particles *= settings.snapshot.mScale
    
    # Handle units
    units = setup_units(m_star, r)
    
    if xyz.units != r.units:
        
        xyz.convert_units(units['x'])
    
    # Initialize arrays
    snapshot = pynbody.new(star=1,gas=nParticles)
    snapshot['vel'].units = units['v']
    snapshot['eps'] = SimArray(0.01, units['x'])
    snapshot['rho'] = 0.
    snapshot['metals'] = metals
    # Assign array values
    snapshot.gas['pos'] = xyz
    snapshot.gas['temp'] = IC.T(r)
    snapshot.gas['mass'] = m_particles
    
    snapshot.star['pos'] = 0.
    snapshot.star['mass'] = m_star
    # Estimate the star's softening length as the closest particle distance/2
    snapshot.star['eps'] = r.min()/2.
    
    return snapshot
