# -*- coding: utf-8 -*-
"""
Defines all the data classes for storing, managing, organizing, and accessing
clumps.

Created on Thu May  7 21:38:26 2015

@author: ibackus
"""

import pynbody as pb
SimArray = pb.array.SimArray
SimSnap = pb.snapshot.SimSnap
import numpy as np
from warnings import warn


def newSimgroup(simlist):
    """
    Generates a new simgroup object from a list of clumplists.
    """    
    nSims = len(simlist)
    
    simgrp = simgroup(nSims)
    
    for i, s in enumerate(simlist):
        
        simgrp[i] = newSim(s)
        
    return simgrp

def newSim(clumplist):
    """
    Generates a sim object from a list of clump dictionaries for an entire
    simulation.
    """
    if len(clumplist) < 1:
        
        return sim(0)
        
    nClumps = len(clumplist)
    s = sim(nClumps)
    
    for i, cl in enumerate(clumplist):
        
        s[i] = newClump(cl)
    
    return s

def newClump(clumpDict):
    """
    Generates a new clump object from a clump dictionary (see clumps.blank_clump
    and clumps.build_clumps).
    """
    # Number of timesteps in this simulation
    nt = len(clumpDict['pos'])
    
    # ---Compatability---
    # Deal with an unfortunate name convention.  pynbody uses 'vel' by default
    # to refer to 3D velocity, whereas clumpDicts may use 'v'
    if 'v' in clumpDict:
        
        clumpDict['vel'] = clumpDict.pop('v')
        
    if 'm' in clumpDict:
        
        clumpDict['mass'] = clumpDict.pop('m')
        
    # NaN values for N are unnatural--make them zeros
    if 'N' in clumpDict:
        
        clumpDict['N'][np.isnan(clumpDict['N'])] = 0
    
    # Initialize a blank clump
    cl = clump(nt)
    
    # Load all the arrays that are created by default in clump.__init__
    unloaded_keys = set(clumpDict.keys())
    for k in cl.keys():
        
        if k in clumpDict.keys():
            
            cl[k] = clumpDict[k]
            unloaded_keys.remove(k)
    
    # Try to load any new arrays which are present in clumpDict
    for k in clumpDict.keys():
        
        if k in unloaded_keys:
            
            unloaded_keys.remove(k)
            v = clumpDict[k]
            if np.ndim(v) > 1:
                
                ndim = v.shape[1]
                
            else:
                
                ndim = 1
                
            cl._create_array(k, ndim, v.dtype)
            cl[k] = v
            
    if len(unloaded_keys) > 0:
        
        warn('Could not load {0} into clump'.format(unloaded_keys))
    
    return cl
    
def clump(nt):
    """
    Generates a clump object, which is basically just a pynbody SimSnap where
    time steps are represented by particles.
    
    **ARGUMENTS**
    
    nt : int
        Number of time steps in the simulation
    """
    
    cl = pb.new(s=nt)
    cl._create_arrays(['r_clump', 'time', 'm', 'T', 'rho'], ndim=1,\
    dtype=float)
    cl._create_arrays(['L'], ndim=3, dtype=float)
    cl._create_array('N', dtype=int)
    
    return cl
    
class sim(np.ndarray):
    """
    The sim class.  sims are basically just numpy arrays containing clump
    objects.  It stores data for all the clumps in a given simulation (for all
    available time steps).  Additionally, slicing can be done in different ways
    (see below)
    
    To initialize a blank sim object:
    
        >>> s = sim(nClumps)
        
    To generate a sim object from a list of clumpDicts for a simulation (ie
    the output of :func:`clumps.clump_tracker`)
    
        >>> s = newSim(clumplist)
        
    Where nClumps is the number of clumps in the simulation
    
    Accessing different quantites:
    
        >>> s = newSim(clumplist) # create sim object from clumplist
        >>> x = s('x', t=slice(None)) # get x for all clumps at all times
        # Get the mass for the first time step at which the clumps exist
        >>> m0 = s('mass',t=0,clumpexist=True)
        # Get velocity for the final timestep at which a clump exists
        >>> vf = s('vel', t=-1, clumpexist=True)
        # Get masses of all clumps at end of simulation
        >>> all_masses = s('mass', t=-1)
        >>> mask = ~np.isnan(all_masses)
        >>> mf = all_masses[mask]
        
    Print available keys:
    
        >>> s.keys()
        
    Select a sub-group of clumps by normal array slicing:
    
        >>> s1 = s[0:10]
        
    **CALL FUNCTION ARGUMENTS**
    s(key, t=None, clumpexist=False, clump=None)
    
    key : str
        Clump key to access
    t : slicer
        Anything that can be used to slice a numpy array, eg [1,2,4] or
        slice(0,-1).  Slices according to the simulation time
    clumpexist : bool
        Filter out all quantities by whether the clump exists at that time.
        e.g., this makes it easy to select the first timestep when a clump
        exists
    clump : slicer
        Anything that can be used to slice a numpy array, eg [1,2,4] or
        slice(0,-1). Slices by clump number
    """
    
    def __new__(subtype, shape, buffer=None, offset=0, \
    strides=None, order=None):
        
        dtype = object
        
        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, \
        strides, order)
        
        return obj
        
    def __array_finalize__(self, obj):
        
        if obj is None: 
            newSim
            return
            
    def __call__(self, key, t=None, clumpexist=False, clump=None):
        """
        """
        # ---------------------------------------
        # Initialize
        # ---------------------------------------
        s = _simslice(self, clump)
        test_array = s[0][0][key]
        units = _getunits(test_array)
        dim = _getdim(test_array)
        nt = _getnt(s, t)
        
        if nt <= 0:
            # Nothing matches t
            return None
        
        nClumps = len(s)
        
        dummy = s[0][0][key]
        dtype = dummy.dtype
        
        outarray = self._init_array(nClumps, nt, dim, dtype, units) 
        usemask = np.ones(nClumps, dtype=bool)
        
        # ---------------------------------------
        # Access the data
        # ---------------------------------------
        for iClump, cl in enumerate(s):
            
            if clumpexist:
                
                # only look at the times where the clump exists
                cl = cl[cl['exists']]
                
            if t is not None:
                
                try:
                    
                    cl = cl[t]
                    
                except IndexError:
                    
                    cl = []
                
            nt_use = len(cl)
            if nt_use > 0:
                
                if nt > 1:
                    
                    outarray[iClump, 0:nt_use] = cl[key]
                    
                else:
                    
                    outarray[iClump] = cl[key]
                    
            else:
                
                usemask[iClump] = False
                
        # ---------------------------------------
        # Filter/return data
        # ---------------------------------------
        if clumpexist:
            
            return outarray[usemask]
            
        else:
            
            return outarray
            
    def __repr__(self):
        
        printstr = "sim object.  {0} clumps, {1} time steps".format(\
        self.nClumps(), self.nt())        
        return printstr
        
    def __str__(self):
        
        return self.__repr__()
        
            
    def _init_array(self, nClumps, nt=0, dim=0, dtype=float, units=None):
        """
        Initialize a blank array for multi-clump slicing.  Default fill values
        are NaN for floats, -1 for ints, and 0 otherwise.  This is useful
        for automatically flagging when clumps do not exist
        
        **ARGUMENTS**
        
        nClumps : int
            Number of clumps
        nt : int
            Number of time steps
        dim : int
            Dimension of the array.  IE 3 for 'vel' or 1 for 'x'
        dtype : dtype
            dtype of the array
        units : (see pynbody.units)
            Units for the array
            
        **RETURNS**
        
        An array of a shape suitable for nClumps, nt, and dim
        """
        shape = [nClumps]
        
        if nt > 1:
            
            shape.append(nt)
            
        if dim > 1:
            
            shape.append(dim)
            
        if np.issubdtype(dtype, float):
            
            fill_val = np.nan
            
        elif np.issubdtype(dtype, int):
            
            fill_val = -1
            
        else:
            
            fill_val = 0
            
        outarray = SimArray(fill_val*np.ones(shape, dtype=dtype), units)
        
        return outarray
        
    def keys(self):
        """
        Return keys present in all clumps
        """
        return _keys(self)
    
    def nClumps(self):
        """
        Returns the number of clumps in the simulation
        """
        
        return len(self)
        
    def nt(self):
        """
        Returns the number of timesteps in the simulation
        """
        if self.nClumps() > 0:
            
            return len(self[0])
            
        else:
            
            return 0
        

        
class simgroup(np.ndarray):
    """
    The simgroup class.  Basically just an array containing a bunch of sim
    objects.  Meant to contain information for clumps in a suite of simulations
    
    To initialize a blank sim group:
    
        >>> sgrp = simgroup(nSims)
        
    To initialize from a list of simulations (ie a list of the outputs from
    :func:`clumps.clumptracker`):
    
        >>> sgrp = newSimgroup(simlist)
        
    Accessing different quantites:
    
        >>> s = newSimgroup(simlist) # create sim object from simlist
        >>> x = s('x', t=slice(None)) # get x for all clumps at all times
        # Get the mass for the first time step at which the clumps exist
        >>> m0 = s('mass',t=0,clumpexist=True)
        # Get velocity for the final timestep at which a clump exists
        >>> vf = s('vel', t=-1, clumpexist=True)
        # Get masses of all clumps at end of simulation
        >>> all_masses = s('mass', t=-1)
        >>> mask = ~np.isnan(all_masses)
        >>> mf = all_masses[mask]
        
    Print available keys:
    
        >>> s.keys()
        
    Select a sub-group of simulations by normal array slicing:
    
        >>> s1 = s[0:10]
        
    **CALL FUNCTION ARGUMENTS**
    s(key, t=None, clumpexist=False, clump=None, sims=None)
    
    key : str
        Clump key to access
    t : slicer
        Anything that can be used to slice a numpy array, eg [1,2,4] or
        slice(0,-1).  Slices according to the simulation time
    clumpexist : bool
        Filter out all quantities by whether the clump exists at that time.
        e.g., this makes it easy to select the first timestep when a clump
        exists
    clump : slicer
        Anything that can be used to slice a numpy array, eg [1,2,4] or
        slice(0,-1). Slices by clump number
    sims : slicer
        Anything that can be used to slice a numpy array, eg [1,2,4] or
        slice(0,-1). Slices by simulation number
    """
    
    def __new__(subtype, shape, buffer=None, offset=0, \
    strides=None, order=None):
        
        dtype = sim
        
        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, \
        strides, order)
        
        return obj
        
    def __array_finalize__(self, obj):
        
        return
            
    def __call__(self, key, t=None, clumpexist=False, clump=None, sims=None):
        """
        Call documentation for simgroups
        """
        # Slice the simulation group is needed
        simgrp = _simgroupslice(self, sims)
        
        # No simulations were selected
        if len(simgrp) < 1:
            
            return None
            
        # Loop through all the simulations at generate the request arrays
        outlist = []
        
        for iSim, s in enumerate(simgrp):
            
            if len(s) > 0:
                # The simulation has a clump.  Now generate the requested array
                val = s(key,t,clumpexist,clump)
                
            else:
                # The simulation has no clumps
                val = None
            
            # Append output (if it's not None)
            if val is not None:
                
                if len(val) > 0:
                    
                    outlist.append(val)
                    
        # Concatenate the list of output arrays into a single SimArray
        outarray = arraylistcat(outlist)
        
        return outarray
        
    def __repr__(self):
        
        printstr = 'simgroup object. {0} simulations, {1} clumps'.format(\
        self.nsim(), self.nClumps())
        
        return printstr
    
    def __str__(self):
        
        return self.__repr__()
        
    def keys(self):
        """
        Return keys present in all clumps
        """
        return _keys(self)
        
    def nsim(self):
        """
        Return the number of simulations present here
        """
        return len(self)
        
    def nClumps(self):
        
        n = 0
        for simulation in self:
            
            n += simulation.nClumps()
            
        return n
        
    
def _keys(obj):
    """
    Return the clump keys present in all things here
    """
    
    k = []
    
    if len(obj) > 0:
        # There is at least one thing
        for x in obj:
            # Make sure keys() is defined
            if hasattr(x, 'keys'):
            
                k.extend(x.keys())
        
        k = list(set(k))
        k.sort()
            
    return k
    
@pb.derived_array
def exists(sim):
    """
    Defines a pynbody derived array which determines the time steps a clump 
    exists at
    """
    return (sim['N'] > 0)
    
def len2(x):
    """
    A modified version of len() to work with numbers.  Numbers have a length
    of 1
    """
    
    if hasattr(x, '__len__'):
        
        length = len(x)
        
    elif isinstance(x, (int,float,long,complex)):
        
        length = 1
        
    return length
    
def arraylistcat(arraylist):
    """
    Concatenate a list of like arrays (or SimArrays) into a single array.
    Concatenates along the first dimension
    Returns None for an empty list
    """
    if len(arraylist) < 1:
        
        return None
        
    nx = 0
    for x in arraylist:
        
        nx += len(x)
        
    dummy = arraylist[0]
    shape = list(dummy.shape)
    shape[0] = nx
    units = _getunits(dummy)
    
    outarray = SimArray(np.zeros(shape), units)
    counter = 0
    
    for array in arraylist:
        
        outarray[counter:counter+len(array)] = array
        counter += len(array)
        
    return outarray
    
def _simslice(simulation, clump=None):
    """
    A method for slicing a sim, guaranteeing a sim object is returned.
    clump is anything that can be used to slice an array
    """
    if clump is None:
        # No slicing required
        s = simulation
        
    else:
        # Slice according to clump
        s = simulation[clump]
    
    if not isinstance(s, sim):
        # Cast s as a sim object
        dummy = sim(1)
        dummy[0] = s
        s = dummy
        
    return s

def _simgroupslice(simgrp, sims=None):
    """
    A method for slicing a simgroup, guaranteeing a simgroup object is returned.
    sims is anything that can be used to slice an array
    """
    if sims is None:
        # No slicing required
        s = simgrp
        
    else:
        # Slice according to sims
        s = simgrp[sims]
    
    if not isinstance(s, simgroup):
        # Cast s as a sim object
        dummy = simgroup(1)
        dummy[0] = s
        s = dummy
        
    return s
    
def _getunits(x):
    """
    Attempts to get the units of x.  If x has not units, None is returned
    """
    if pb.units.has_units(x):
        
        units = x.units
        
    else:
        
        units = None
        
    return units
    
def _getdim(x):
    """
    Get the dimension of an array x.  IE, for 'vel' dim=3, for 'z', dim=1
    
    For x shape (N,)                dim = 0
    For x shape (N,m)               dim = m
    For x shape (N1, N2, ...,m)     dim = m
    """
    
    if np.ndim(x) > 1:
        
        dim = x.shape[-1]
        
    else:
        
        dim = 0
    
    return dim
    
def _getnt(simulation, t=None):
    """
    Get the total number of time steps the slicer t will create for a simulation
    """
    nt_sim = simulation.nt()
    
    if t is not None:
        
        dummy = np.zeros(nt_sim)
        nt = len2(dummy[t])
            
    else:
        
        nt = nt_sim
        
    return nt