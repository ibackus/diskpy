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

# diskpy packages
import diskpy
from diskpy import global_settings
from diskpy.utils import match_units, strip_units
from diskpy.pychanga import make_director, make_param, setup_units
import calc_velocity

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
    # Finalize
    print 'Wrapping up'
    setup_sinks(snapshot, param, r_sink = IC.pos.r.min())
    
    return snapshot, param, director

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
    
def setup_sinks(snapshot, param, r_sink=None):
    """
    Sets up snapshot and param for stars that are sinks.
    
    Parameters
    ----------
    snapshot : SimSnap
    param : param dict
    r_sink : SimArray-like
        If provide, this will be the sink radius.  Otherwise it's taken to be
        the minimum gas cylindrical radius
    
    Returns
    -------
    None
    """
    units = diskpy.pychanga.units_from_param(param)
    
    if r_sink is None:
        
        r_sink = snapshot.g['rxy'].min()
        r_sink.convert_units(units['l_unit'])
        
    r_sink = strip_units(r_sink)
    # Set the star tforms to a negative number.  This allows UW ChaNGa treat 
    # stars as sink particles
    snapshot.star['tform'] = -1.0    
    # Update params
    Mstar = snapshot.s['mass'].min()
    Mstar.convert_units(units['m_unit'])
    param['dSinkBoundOrbitRadius'] = r_sink
    param['dSinkRadius'] = r_sink
    param['dSinkMassMin'] = 0.9 * strip_units(Mstar)
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
    # Estimate the star's softening length as the closest particle distance
    snapshot.star['eps'] = r.min()
    
    return snapshot
