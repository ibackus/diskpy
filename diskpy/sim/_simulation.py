# -*- coding: utf-8 -*-
"""
Created on Wed May 11 15:52:34 2016

@author: ibackus
"""
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
Profile = pynbody.analysis.profile.Profile
import diskpy
import os


class Simulation():
    
    def __init__(self, directory, paramname=None, logfile=None, 
                 notefile='notes.txt', nbins=40):
        
        self.directory = None
        self.param = None
        self.paramname = None
        self.IC = None
        self.snaps = None
        self.profs = None
        self.nbins = nbins
        self.t = None
        self._redges = None
        
        
    def __delitem__(self, key):
        
        raise RuntimeError, "Cannot delete"
        
    def __getitem__(self, key):
        
        return self.snaps[key]
        
    def __setitem__(self, key, value):
        
        self.snaps[key] = value
        
    def __len__(self):
        
        return len(self.snaps)
        
        
    def set_nbins(self, nbins):
        
        self.nbins = nbins
        self.profs, self._redges = setupProfiles(self.snaps, self.nbins)
        
            
def loadSim(paramname, useAbsPath=True, nbins=40):
    
    sim = Simulation(nbins)
    # Get the directory of the simulation
    sim.directory = os.path.dirname(paramname)    
    if useAbsPath:
        
        sim.directory = os.path.abspath(sim.directory)
        
    # Load the param file
    sim.param = diskpy.utils.configparser(paramname, 'param')
    sim.paramname = paramname
    # Get the units
    sim.units = diskpy.pychanga.units_from_param(sim.param)
    # Load the simulation list
    sim.IC, sim.snaps, sim.filenames = loadSimList(sim.directory, sim.param, paramname, useAbsPath)     
    # Get the times of each snapshot
    print 'Loading sim times'
    sim.t = loadTimes(sim.snaps, timeUnits=sim.units['t_unit'])
    # Get length of each snapshot
    sim.nparticles = np.zeros(len(sim.snaps), dtype=int)
    for i, snap in enumerate(sim.snaps):
        
        sim.nparticles[i] = len(snap)
    
    # initialize radial profiles
    sim.profs, sim._redges = setupProfiles(sim.snaps, sim.nbins)
    
    # set-up molecular weight for the simulation
    molecularWeight = sim.param.get('dMeanMolWeight', 1)
    molecularWeight = SimArray(molecularWeight, 'm_p')
    sim.mu = molecularWeight
    
    return sim
    
def setupProfiles(simlist, nbins=40):
    
    snap = simlist[0]
    redges = diskpy.pdmath.setupbins(snap['r'], bins=nbins)
    profs = []
    
    for snap in simlist:
        
        print 'Setting up profile: ', snap
        prof = Profile(snap, ndim=3, bins=redges)
        profs.append(prof)
        
    return profs, redges
    
def loadSimList(directory, param, paramname, useAbsPath=True):
    
    
    # Load initial conditions
    ICname = os.path.join(directory, param['achInFile'])
    IC = pynbody.load(ICname, paramname=paramname)
    # Load simulation outputs
    fprefix = param['achOutName']
    fnames = diskpy.pychanga.get_fnames(fprefix, directory)
    snapshots = []
    
    for fname in fnames:
        
        print 'Loading sim: ', fname
        snapshot = pynbody.load(fname, paramname=paramname)
        snapshots.append(snapshot)
        
    return IC, snapshots, fnames

def loadTimes(simlist, timeUnits='yr'):
    
    t = SimArray(np.zeros(len(simlist)), timeUnits)
    for i, sim in enumerate(simlist):
        
        t[i] = diskpy.pychanga.snapshot_time(sim).in_units(timeUnits)
        
    return t