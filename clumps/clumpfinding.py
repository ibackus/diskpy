# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 12:10:08 2014

@author: ibackus
"""

# Multiprocessing modules
from multiprocessing import Pool, cpu_count

# Generic pacakges
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import subprocess
import glob
import os
import matplotlib.pyplot as plt
import re

# 'Internal' packages
import isaac
from diskpy import utils

def clump_tracker(fprefix, param=None, directory=None, nsmooth=32, verbose=True):
    """
    Finds and tracks clumps over a simulation with multiple time steps and
    calculates various physical properties of the clumps.
    
    Runs all the steps necessary to find/track clumps, these are:
    
    get_fnames
    pFind_clumps
    pClump_properties
    pLink2
    multilink
    build_clumps
    
    If the iord property is not found, the linking will only work if the number
    of particles remains constant through the simulation
    
    **ARGUMENTS**
    
    fprefix : str
        Prefix of the simulation outputs
    param : str (recommended)
        Filename of a .param file for the simulation
    directory : str (optional)
        Directory to search through.  Default is current working directory
    nsmooth : int (optional)
        Number of nearest neighbors used for particle smoothing in the
        simulation.  This is used in the definition of a density threshold
        for clump finding.
    verbose : bool (optional)
        Verbosity flag.  Default is True
        
    **RETURNS**
    
    clump_list : list
        A list containing dictionaries for all clumps foujohn obryan fiddlend in the simulation
        See clump_properties for a list of the properties calculated for clumps
    """
    
    # Get a list of all snapshot files
    fnames = get_fnames(fprefix, directory)
    nfiles = len(fnames)
    
    # Run the clump (halo) finder
    if verbose: print "\n\nRunning clump finder on {} files\n\n".format(nfiles)
    clumpnum_list = pFind_clumps(fnames, nsmooth, param, verbose=verbose)
    nclumps = np.zeros(nfiles, dtype=int)
    
    for i, clumpnums in enumerate(clumpnum_list):
        
        nclumps[i] = clumpnums.max()
        
    if nclumps.max() <= 0:
        
        if verbose: print 'No clumps found'
        return []
    
    # Calculate the physical properties of the clumps
    if verbose: print "\n\nCalculating the physical of properties of clumps\n\n"
    properties = pClump_properties(fnames, clumpnum_list)
    
    # Link clumps on consecutive time-steps
    if verbose: print "\n\nLinking Clumps\n\n"
    link_list = pLink2(properties)
    # Link on multiple time-steps
    multilink_list = multilink(link_list)
    # Build the clumps
    clump_list = build_clumps(multilink_list, properties, fnames, param)
    
    return clump_list
    

def get_fnames(fprefix, directory=None):
    """
    Finds the filenames of ChaNGa simulation outputs.  They are formatted as:
        fprefix.000000
    i.e., fprefix followed by '.' and a 6 digit number
    
    **ARGUMENTS**
    
    fprefix : str
        Prefix of the simulation outputs
    directory : str (optional)
        Directory to search through.  Default is current working directory
        
    **RETURNS**
    
    fnames : list
        A list containing the matching filenames
    """
    
    fnames = []
    
    if directory is not None:
        
        fprefix = os.path.join(directory, fprefix)
        
    repattern = '^' + fprefix + '.(?:(?<!\d)\d{6}(?!\d))$'
        
    for fname in glob.glob(fprefix + '.*'):
        
        if re.search(repattern, fname) is not None:
            
            fnames.append(fname)
            
    fnames.sort()
    
    return fnames

def blank_clump(clump_pars, nt=1):
    """
    Generates a blank clump dictionary, using clump_pars
    """
    
    ignore_keys = ['iord', 'ids']
    keys = clump_pars.keys()
    
    # Remove the ignored keys
    for ignore_key in ignore_keys:
        
        if ignore_key in keys:
            
            keys.remove(ignore_key)
            
    clump = {}
    
    for key in keys:
        
        # Load the current value
        val = clump_pars[key]
        # All values should be nd arrays or nd simarrays
        shape = list(val.shape)
        shape[0] = nt
        # Initialize a blank array with the same dtype
        val_array = np.ones(shape) * np.nan
        # Check if there are units
        if pynbody.units.has_units(val):
            
            val_array = SimArray(val_array, val.units)
            
        clump[key] = val_array
        
    return clump
        

def build_clumps(multilink_list, clump_pars_list, fnames=None, param=None):
    """
    Builds a list of clump dictionaries.  The clump dictionaries contain various
    properties as a function of time for the different clumps.
    
    **ARGUMENTS**
    
    multilink_list : list
        A list of clump-link arrays (output of multilink)
    clump_pars_list : list
        List of dictionaries containing clump properties for all clumps on a
        given time step (see pClump_properties)
    fnames : list (recommended)
        A list of the simulation snapshot filenames.  If provided, the time
        steps are included in the clumps
    param : str (recommended)
        Filename of a .param file for the simulation.  Only used if fnames is
        provided
        
    **RETURNS**
    
    clump_list : list
        A list of clump dictionaries.  Each clump dictionary gives various
        physical properties for a clump as a function of time.  A value of NaN
        means the clump was not present at that time step.
    """
    
    # Initialize list to contain all the clump objects (dictionaries)
    clump_list = []
    
    nt = len(clump_pars_list) # number of timesteps
    nclumps = len(multilink_list) # number of clumps
    
    if nclumps < 1:
        
        return clump_list
    
    # Find the time step of each simulation
    if fnames is not None:
        
        t = SimArray(np.zeros(nt), 'yr')
        
        for i, fname in enumerate(fnames):
            
            f = pynbody.load(fname, paramname=param)
            t_unit = SimArray(1, f.infer_original_units('yr'))
            # Note, for ChaNGa outputs, t0 = t_unit (ie, 1 in simulation units)
            # To correct for this, we subtract off one time unit from the 
            # snapshot's time
            t[i] = f.properties['time'] - t_unit
            
    
    # Start by finding the first non-zero clump_pars (ie, first timestep with
    # clumps)
    for i, clump_pars in enumerate(clump_pars_list):
        
        if clump_pars is not None:
            
            iFirstClump = i
            break
        
    
    # Now fill in the clumps
    for iClump in range(nclumps):
        print iClump
        
        # Initialize a blank clump
        clump = blank_clump(clump_pars_list[iFirstClump], nt)
        
        for iord, iStep in multilink_list[iClump]:
            
            clump_pars = clump_pars_list[iStep]
            
            for key in clump.keys():
                
                clump[key][iStep] = clump_pars[key][iord]
                
        clump['time'] = t.copy()
        clump_list.append(clump)
        
    return clump_list
        

def multilink(link_list):
    """
    Links clumps on multiple time steps.
    
    Given the output of link2 or pLink2 for multiple time-steps, this determines
    every clump's ID as a function of time-step.
    
    **ARGUMENTS**
    
    link_list : list
        A list of link arrays.  link_list[i] contains the links between time
        step i and i+1.
        ie: link_list[i] = link2(clump_pars_list[i], clump_pars_list[i+1])
        Same as the output from pLink2
        
    **RETURNS**
    
    clumpid_list : list
        A list of 2D arrays.  Each array gives pairs of (clumpID, time-step)
    """
    
    n_links = len(link_list)
    
    clump_list = []
    
    iStart = 0
    
    while iStart < n_links:
        
        if link_list[iStart] is not None:
        
            pairs0 = link_list[iStart]
            new_mask = pairs0[:,0] == -1
            new_iord = pairs0[new_mask,1]
        
            for iord0 in new_iord:
                
                t = iStart + 1
                iPair = iStart + 1
                # Initialize a new clump
                clump = [ [iord0, t] ]
                
                while iPair < n_links:
                    
                    pairs1 = link_list[iPair]
                    
                    if pairs1 is None:
                        # There are no clumps at this time step
                        break
                    
                    t = iPair + 1
                    iord = clump[-1][0]
                    
                    new_ind = np.nonzero(pairs1[:,0] == iord)[0]
                    if len(new_ind) > 0:
                        # The clump links to something in the next timestep
                        new_ind = int(new_ind)
                        # Add the new index to the clump
                        #clump.append([pairs1[new_ind, 1], t])
                        clump.append([new_ind, t])
                        # Increment the time step
                        iPair += 1
                        
                    else:
                        # The clump links to nothing.  It has died
                        break
                    
                clump_list.append(np.array(clump))
                
        iStart += 1
                
    return clump_list
    
def _parallel_link2(args):
    """
    A parallel wrapper for link2
    """
    
    return link2(*args)
    
def pLink2(clump_pars_list):
    """
    A parallel (batch) implementation of links2 for linking clumps in
    consecutive time-steps.  
    
    **ARGUMENTS**
    
    clump_pars_list : list
        A list of dictionaries containing the properties of clumps.  Each
        element of the list is a dictio:nary for a single time step.
        (see pClump_properties)
        
    **RETURNS**
    
    link_list : list
        A list of link arrays.  link_list[i] contains the links between time
        step i and i+1.
        ie: link_list[i] = link2(clump_pars_list[i], clump_pars_list[i+1])
    """
    
    arg_list = zip(clump_pars_list[0:-1], clump_pars_list[1:])
    nproc = cpu_count()
    
    pool = Pool(nproc)
    link_list = pool.map(_parallel_link2, arg_list)
    pool.close()
    pool.join()
    
    return link_list

def link2(clump_pars1, clump_pars2, link_thresh = 0.2):
    """
    'Links' the clumps in two consecutive timesteps.  i.e., a clump numbered
    10 in one time step might be numbered 15 in the next.  Requires the particle
    iord (ie, particle ID)
    
    **ARGUMENTS**
    
    clump_pars1 : dict
        Dict containing the properties of the clumps in time-step 1 
        (see clump_properties)
    clump_pars2 : dict
        Dict containing the properties of the clumps in time-step 2 
        (see clump_properties)
    link_thresh : int (optional)
        The minimum fraction of particles in clump i at time 1 which must end
        up in clump j at time 2 to consider the clumps 'linked'
        
    **RETURNS**
    
    clump_pairs : array
        2D numpy array organized according to (parent index, clump number)
        where clump number is the number of a clump in time 2 and parent index
        is the number of it's "parent" in time-step 1 (ie, that same clump's
        ID in time-step 1)
        A parent index of -1 corresponds to no parent (a new clump)
    """
    
    if clump_pars2 is None:
        # there are no clumps in the second time step.  Any clumps in the first
        # time step must have died.  Return None
        return
    
    if clump_pars1 is None:
        # Any clumps in the second time step are new.  This is handled by 
        # saying their parents are -1
        n2 = len(clump_pars2['iord'])
        clump_pairs = np.zeros([n2, 2], dtype=int)
        clump_pairs[:,0] = -1
        clump_pairs[:,1] = np.arange(n2)
        
        return clump_pairs
    
    # number of clumps in each time step
    n1 = len(clump_pars1['iord'])
    n2 = len(clump_pars2['iord'])
    
    iord_list1 = list(clump_pars1['iord'])
    iord_list2 = list(clump_pars2['iord'])
    
    # Used to store how many particles clumps have in common
    connections = np.zeros([n1, n2], dtype=int)
    
    # Calculate the number of connections common to the clumps in clump_pars1
    # and the clumps in clump_pars2
    # Loop over the first set of clumps
    for i, iord1 in enumerate(iord_list1):
        # Loop over the second set of clumps
        iord1.sort()
        npart1 = len(iord1)
        for j, iord2 in enumerate(iord_list2):
            
            # Find which particles are common to clump[i] in pars1 and clump[j]
            # in pars2
            intersect = np.intersect1d(iord1, iord2, assume_unique=True)
            # Save how many are shared in the 2 clumps
            connections[i,j] = len(intersect)
            # Now only retain particles that are not common to the clumps
            # IE, we know where these particles end up, we can stop checking
            # them
            iord1 = np.setdiff1d(iord1, intersect, assume_unique=True)
            iord2 = np.setdiff1d(iord2, intersect, assume_unique=True)
            iord_list2[j] = iord2
            
            if len(iord1) < 1:
                # There are no more particles to look at in the original clump                
                break
            
        # Now ignore any connections where number of particles shared between
        # clumps less than link_thresh * (num. part. in clump 1)
        thresh_mask = connections[i,:] < link_thresh * npart1
        connections[i, thresh_mask] = 0
    
    # Find the clump in clump_pars2 that shares the most number of members for
    # a clump in clump_pars1.  This gives us the children of the clumps in
    # clump_pars1
    col_ind = connections.argmax(1)
    # Set all others to 0
    mask = np.zeros([n1,n2], dtype=bool)
    row_ind = np.arange(n1)
    mask[row_ind, col_ind] = True
    connections[~mask] = 0
    
    # The clumps in clump_pars2 may have multiple parents.  Select the parent
    # which shares the most particles in common (note, if this isn't the
    # child of any clump, parent_index will always be 0)
    parent_index = connections.argmax(0)
    # and find the number of particles inherited from the parent
    n_inherit = connections.max(0)
    
    # Demand to inherit at least 1 particle from the parent
    parent_index[n_inherit < 1] = -1
    
    # Now create the clump pairs
    clump_pairs = np.array( [parent_index, np.arange(n2)] ).T

    return clump_pairs
        

def clump_im(f, clump_array, width, qty='rho', resolution=1200, clim=None, clump_min=None):
    """
    Plots an sph image from f with particles not in clumps colored red and 
    particles in clumps colored green.  Uses pynbody for a backend.
    
    **ARGUMENTS**
    
    f : TipsySnapshot (see pynbody) or str
        The snapshot to plot OR the filename of a snapshot to plot
    clump_array : numpy array
        A array (same length as f) such that 0 entries correspond to particles
        not in clumps and entries > 0 correspond to particles in clumps
    width : str, simarray
        See pynbody.plot.sph.image.  Width of the plot, ie '3 au'
    resolution : int
        Resolution in pixels of the plot.  The plot will be res by res pixels
    clim : tuple,list,array
        Density limits arranged as [low, high].  Any pixels below low are mapped
        to 0, any pixels above high are mapped to 1.
    clump_min : float
        Used to set a cutoff for masking the clumps.  Not needed
        
    **RETURNS**
    
    image : numpy nd-array
        Returns an NxNx3 numpy array of the color image plotted.
    """
    # Check if the input is a filename
    if isinstance(f, str):
        
        f = pynbody.load(f)
    
    # Get the current state for matplotlib (this still needs work, since an
    # extra window will in general be created)
    current_fig = plt.gcf()
    interactive_flag = plt.isinteractive()
    plt.ioff()
    
    # Intermediate figure, required for rendering the plots
    fig1 = plt.figure()
    
    # Render a grayscale image of all the particles
    im_all = pynbody.plot.sph.image(f.g, width=width,resolution=resolution, cmap='gray', qty=qty)
    xlim = plt.xlim()
    ylim = plt.ylim()
    extent = [xlim[0], xlim[1], ylim[0], ylim[1]]
    
    fig1.clf()
    
    # Initialize the color image
    im_color = np.zeros([resolution, resolution, 3])
    # Set the red channel of the color image to be the plot of all particles
    im_color[:,:,0] = np.log(im_all)
    
    # Check to see that there is at least one clump
    clump_flag = (clump_array.max() > 0)
    
    if clump_flag:
        
        # Get a sub-snap of just particles in a clump
        mask = clump_array > 0
        f2 = f[mask]
        # Render an image of just particles in a clump
        im_clump = pynbody.plot.sph.image(f2.g, width=width,resolution=resolution, cmap='gray',qty=qty)
        # Set the clump image as the green channel
        im_color[:,:,1] = np.log(im_clump)
        
        plt.clf()
        
        # Most of the clump image should be blank space: igore everything
        # below a threshold
        if clump_min is None:
        
            clump_min = im_clump.mean()
        
        mask2 = im_clump > clump_min
    
    # Set the color scaling
    if clim is None:
        
        clim = [im_all.min(), im_all.max()]
        
    log_clim = [np.log(clim[0]), np.log(clim[1])]
            
    im_color -= log_clim[0]
    im_color /= (log_clim[1] - log_clim[0])
    im_color[im_color < 0] = 0
    im_color[im_color > 1] = 1
    
    if clump_flag:
        
        # Set all pixels outside a clump (in the clump image) to 0
        im_color[~mask2,1] = 0
        # Set all pixels inside a clump (in the overall image) to 0
        im_color[mask2,0] = 0
    
    else:
            
        im_color[:,:,1] = 0
        
    im_color[:,:,2] = 0
    
    # Plot
    plt.figure(current_fig.number)
    
    if interactive_flag:
        plt.ion()
    
    plt.imshow(im_color, extent=extent, interpolation='none', aspect='equal')
    
    # Echo the color limits used
    print 'clims used: {}'.format(clim)
    
    plt.close(fig1)
        
    return im_color
    
def _parallel_clump_pars(args):
    """
    A wrapper to parallelize clump_properties
    """
    
    return clump_properties(*args)
    
def pClump_properties(flist, clumpnum_list):
    """
    A parallel (batch) implementation of clump_properties.  Calculates the 
    physical properties of clumps in a list of snapshots.
    
    **ARGUMENTS**
    
    flist : list
        A list of tipsy snapshots or filenames pointing to snapshots.
    clumpnum_list : list
        A list of arrays (one per snapshot) which define the clump number
        each particle belongs to (see pFind_clumps)
    
    **RETURNS**
    
    properties : list
        A list of dictionaries which contain the clump properties for every
        simulation (see clump_properties)
    """
    
    nproc = cpu_count()
    
    arg_list = zip(flist, clumpnum_list)
    
    pool = Pool(nproc)
    
    properties = pool.map(_parallel_clump_pars, arg_list)
    pool.close()
    pool.join()
    
    return properties

def clump_properties(f, clump_nums):
    """
    Calculates the physical properties of clumps in a snapshot.
    
    **ARGUMENTS**
    
    f : str -OR- tipsy snapshot
        Either a tipsy snapshot or a filename pointing to a snapshot
    clump_nums : array like
        Clump number that each particle belongs to (see clumps.find_clumps).
        0 corresponds to not being in a clump.
    
    **RETURNS**
    
    properties : dict
        A dictionary containing the physical properties of the clumps.
        Keys are:
        
        'm'         mass
        'N'         Number of particles
        'pos'       xyz position
        'r'         cylindrical radial position
        'v'         center of mass velocity
        'L'         Angular momentum relative to clump center of mass
        'T'         Average temperature
        'rho'       Average density
        'r_clump'   Clump radius.  Sqrt of mass averaged particle distance squared
                    (from the center of mass).  IE: r = sqrt( sum(mr^2)/sum(m))
        'ids'       particle IDs in the clump (first particle in simulation is
                    0, second is 1, etc...)
        'iord'      Particle iord (a particle's ID for the whole simulation)
        
    """
    
    if clump_nums.max() < 1:
        # Return none if there are no clumps
        return
    
    if isinstance(f, str):        
        f = pynbody.load(f)
        
    try:
        
        iorder = f['iord']
        
    except KeyError:
        
        print 'Warning.  iorder not found.  Assuming 0,1,2,3...'
        iorder = np.arange(len(f))
        
    particle_nums = np.arange(len(f))
        
    # Only include particles in a clump AND that are not star particles
    mask1 = clump_nums > 0
    n_star = len(f.s)
    mask1[-(n_star+1):-1] = False
    clump_nums1 = clump_nums[mask1]
    f1 = f[mask1]
    iorder1 = iorder[mask1]
    particle_nums1 = particle_nums[mask1]
    
    # Get units set up
    m_unit = f1['mass'].units
    l_unit = f1['pos'].units
    v_unit = f1['vel'].units
    rho_unit = f1['rho'].units
    
    # Get arrays of pointers to the required quantities
    f_mass = f1['mass']
    f_pos = f1['pos']
    f_v = f1['vel']
    f_T = f1['temp']
    f_rho = f1['rho']
    
    # Initialize arrays
    n_clumps = clump_nums1.max()
    
    m = SimArray(np.zeros(n_clumps), m_unit) # clump mass
    N = np.zeros(n_clumps, dtype=int) # Number of particles/clump
    pos = SimArray(np.zeros([n_clumps,3]), l_unit) # center of mass
    r = SimArray(np.zeros(n_clumps), l_unit) # center of mass radial position
    v = SimArray(np.zeros([n_clumps, 3]), v_unit) # center of mass velocity
    # Angular momentum around the center of mass rest frame
    L = SimArray(np.zeros([n_clumps, 3]), m_unit*l_unit*v_unit)
    T = SimArray(np.zeros(n_clumps), 'K') # mass averaged temperature
    rho = SimArray(np.zeros(n_clumps), rho_unit) # density
    r_clump = SimArray(np.zeros(n_clumps), l_unit) # clump radius (size)
    
    # index of each particle (in this file)
    particle_ids = []
    # universal identity of each particle
    particle_iord = []
    
    # loop over the clumps
    for i in range(n_clumps):
        
        mask2 = (clump_nums1 == i+1)
        
        # Mask the input arrays to look at only the current clump
        p_mass = f_mass[mask2]
        p_pos = f_pos[mask2]
        p_v = f_v[mask2]
        p_T = f_T[mask2]
        p_rho = f_rho[mask2]
        
        # Calculate properties of the clump
        N[i] = mask2.sum()
        m[i] = p_mass.sum()
        pos[i] = np.dot(p_pos.T, p_mass[:,None]).flatten()
        pos[i] /= float(m[i])
        r[i] = np.sqrt((pos[i]**2).sum())
        v[i] = np.dot(p_v.T, p_mass[:,None]).flatten()
        v[i] /= float(m[i])
        
        # position of all particles relative to center of mass
        cm_pos = p_pos - pos[i]
        # velocity of all particles relative to center of mass
        cm_v = p_v - v[i]
        # angular momentum of all particles relative to center of mass
        cm_momentum = (cm_v * p_mass[:,None])
        p_L = np.cross(cm_pos, cm_momentum)
        # Total angular momentum relative to center of mass
        L[i] = p_L.sum(0)
        
        T[i] = p_T.sum()/N[i]
        rho[i] = p_rho.sum()/N[i]
        
        # Clump radius
        try:
            r_clump[i] = np.sqrt((p_mass*( (cm_pos**2).sum(1) )).sum()/m[[i]])
        except pynbody.units.UnitsException:
            print 'i is: {}'.format(i)
            return p_mass, cm_pos, m
        
        particle_ids.append(particle_nums1[mask2])
        particle_iord.append(iorder1[mask2])
        
    properties = {'m':m, 'N':N, 'pos':pos, 'r':r, 'v':v, 'L':L, 'T':T, 'rho':rho,\
    'r_clump': r_clump, 'ids': particle_ids, 'iord': particle_iord}
    
    return properties

def _parallel_find_clumps(args):
    """
    A wrapper to parallelize find_clumps()
    """    
    return find_clumps(*args)
    
def pFind_clumps(f_list, n_smooth=32, param=None, arg_string=None, verbose=True):
    """
    A parallel implementation of find_clumps.  Since SKID is not parallelized
    this can be used to run find_clumps on a set of snapshots from one
    simulation.
    
    **ARGUMENTS**
    
    f_list : list
        A list containing the filenames of snapshots OR the tipsy snapshots
    n_smooth : int (optional)
        Number of nearest neighbors used for particle smoothing in the
        simulation.  This is used in the definition of a density threshold
        for clump finding.
    param : str (optional)
        filename for a tipsy .param file
    arg_string : str (optional)
        Additional arguments to be passed to SKID.  Cannot use -tau, -d, -m, -s, -o
    verbose : bool
        Verbosity flag.  Default is True
        
        
    **RETURNS**
    
    clumpnum_list : list
        A list containing the particle clump assignment for every snapshot in 
        f_list.  clumps[i][j] gives the clump number for particle j in
        snapshot i.
    """
    # Number of processes to create = number of cores
    n_proc = cpu_count()
    
    # Set up the arguments for calls to find_clumps
    arg_list = []
    
    for i, f_name in enumerate(f_list):
        
        arg_list.append([f_name, n_smooth, param, arg_string, i, verbose])
        
    print arg_list
    
    # Set up the pool
    pool = Pool(n_proc)
    
    # Run the job in parallel
    results = pool.map(_parallel_find_clumps, arg_list, chunksize=1)
    pool.close()
    pool.join()
    
    return results

def find_clumps(f, n_smooth=32, param=None, arg_string=None, seed=None, verbose=True):
    """
    Uses skid (https://github.com/N-BodyShop/skid) to find clumps in a gaseous
    protoplanetary disk.  
    
    The linking length used is equal to the gravitational softening length of
    the gas particles.
    
    The density cut-off comes from the criterion that there are n_smooth particles
    within the Hill sphere of a particle.  This is formulated mathematically as:
    
        rho_min = 3*n_smooth*Mstar/R^3
        
    where R is the distance from the star.  The trick used here is to multiply
    all particle masses by R^3 before running skid so the density cut-off is:
    
        rho_min = 3*n_smooth*Mstar
        
    **ARGUMENTS**
    
    *f* : TipsySnap, or str
        A tipsy snapshot loaded/created by pynbody -OR- a filename pointing to a
        snapshot.
    
    *n_smooth* : int (optional)
        Number of particles used in SPH calculations.  Should be the same as used
        in the simulation.  Default = 32
    
    *param* : str (optional)
        filename for a .param file for the simulation
    
    *arg_string* : str (optional)
        Additional arguments to be passed to skid.  Cannot use -tau, -d, -m, -s, -o
    
    *seed* : int
        An integer used to seed the random filename generation for temporary
        files.  Necessary for multiprocessing and should be unique for each
        thread.
        
    *verbose* : bool
        Verbosity flag.  Default is True
    
    **RETURNS**
    
    *clumps* : array, int-like
        Array containing the group number each particle belongs to, with star
        particles coming after gas particles.  A zero means the particle belongs
        to no groups
    """
    # Parse areguments
    if isinstance(f, str):
        
        f = pynbody.load(f, paramfile=param)
        
    if seed is not None:
        
        np.random.seed(seed)
        
    # Estimate the linking length as the gravitational softening length
    tau = f.g['eps'][0]
    
    # Calculate minimum density
    rho_min = 3*n_smooth*f.s['mass'][0]
    
    # Center on star.  This is done because R used in hill-sphere calculations
    # is relative to the star
    star_pos = f.s['pos'].copy()
    f['pos'] -= star_pos
    
    # Scale mass by R^3
    R = utils.strip_units(f['rxy'])
    m0 = f['mass'].copy()
    f['mass'] *= (R+tau)**3
    
    # Save temporary snapshot
    f_prefix = str(np.random.randint(np.iinfo(int).max))
    f_name = f_prefix + '.std'
    
    # Save temporary .param file
    if param is not None:
        
        param_name = f_prefix + '.param'
        param_dict = utils.configparser(param, 'param')
        utils.configsave(param_dict, param_name)
        
    f.write(filename=f_name, fmt=pynbody.tipsy.TipsySnap)
        
    f['pos'] += star_pos
    f['mass'] = m0
    
    command = 'totipnat < {} | skid -tau {:.2e} -d {:.2e} -m {:d} -s {:d} -o {}'\
    .format(f_name, tau, rho_min, n_smooth, n_smooth, f_prefix)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if verbose:
        
        for line in iter(p.stdout.readline, ''):
            print line,
            
    p.wait()
    
    # Load clumps
    clumps = isaac.loadhalos(f_prefix + '.grp')
    
    # Cleanup
    for name in glob.glob(f_prefix + '*'):
        
        os.remove(name)
        
    return clumps