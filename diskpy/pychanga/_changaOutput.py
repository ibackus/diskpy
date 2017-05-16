# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 23:26:04 2015

@author: ibackus
"""
import glob
import os
import re
import warnings
import subprocess
import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import datetime

from diskpy.utils import configparser

BUF_SIZE = int(1e6)

def snapshot_time(f, paramname=None):
    """
    Gets the physical time of a snapshot.  t=0 corresponds to the initial
    conditions.
    """
    if isinstance(f, str):
        
        f = pynbody.load(f, paramname=paramname)
        
    t_unit = SimArray(1, f.infer_original_units('yr'))
    # Note, for ChaNGa outputs, t0 = t_unit (ie, 1 in simulation units)
    # To correct for this, we subtract off one time unit from the 
    # snapshot's time
    t = f.properties['time'] - t_unit
    
    return t

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

def load_acc(filename, param_name = None, low_mem = False):
    """
    Loads accelerations from a ChaNGa acceleration file (.acc2), ignoring the
    star particle.  ASSUMES A SINGLE STAR PARTICLE

    IF param_name is None, a .param file is searched for, otherwise param_name
    should be a string specifying a .param file name

    IF no param_file is found, the defaults are used:
        length unit: AU
        mass unit  : Msol

    Setting low_mem=True decreases memory usage by about 2x but also increases
    readtime by about 2x
    """

    if param_name is None:

        prefix = filename.split('.')[0]

        param_list = glob.glob('*' + prefix +'*param')

        if len(param_list) > 0:

            param_name = param_list[0]

        elif len(glob.glob('*.param')) > 0:

            param_name = glob.glob('*.param')[0]

        else:

            warnings.warn('Could not find .param file.  Assuming default units')

    if param_name is not None:

        # If a param name is set or a param file has been found:
        print 'Loading param file: {}'.format(param_name)
        param = configparser(param_name, ftype='param')

    else:

        # Set the default parameters
        param = {}
        # Assume AU as length unit
        param['dKpcUnit'] = pynbody.units.au.ratio('kpc')
        # Assume mass units as Msol
        param['dMsolUnit'] = 1.0

    # Figure out units
    G = pynbody.units.G
    l_unit = param['dKpcUnit']*pynbody.units.kpc
    m_unit = param['dMsolUnit']*pynbody.units.Msol
    t_unit = ((l_unit**3) * G**-1 * m_unit**-1)**(1,2)
    a_unit = l_unit * t_unit**-2

    if low_mem:
        # Perform ASCII read
        with open(filename, 'r') as f:
            # First line in tipsy format is a header
            n_particles = int(f.readline().strip())
            # Pre-allocate for speed
            acc = SimArray(np.zeros(3*n_particles, dtype=np.float32), a_unit)
            # Buffered read
            i0 = 0
            tmp_lines = f.readlines(BUF_SIZE)
            while tmp_lines:
                acc[i0:i0 + len(tmp_lines)] = tmp_lines
                i0 += len(tmp_lines)
                tmp_lines = f.readlines(BUF_SIZE)
                
        return acc.reshape([n_particles, 3], order='F')[0:-1]

    else:

        # Load acceleration file as numpy array
        acc = np.genfromtxt(filename, skip_header=1).astype(np.float32)
        n_particles = len(acc)/3
        # Reshape and make it a SimArray with proper units
        acc = SimArray(acc.reshape([n_particles, 3], order='F'), a_unit)

        return acc[0:-1]
        
def walltime(filename, verbose=True):
    """
    Reads walltime information from a ChaNGa .log file.

    ** ARGUMENTS **

    filename : str
        Filename of the .log file to load

    ** RETURNS **

    wall_per_step : array
        Wall time per step in seconds
    """

    log_file = np.genfromtxt(filename, comments='#', delimiter=' ')
    wall_per_step = log_file[:,-1]
    walltime_total = datetime.timedelta(seconds = wall_per_step.sum())
    walltime_avg = datetime.timedelta(seconds = wall_per_step.mean())
    
    if verbose:
        print 'Total walltime: '
        print str(walltime_total)
        print 'Average walltime per step:'
        print str(walltime_avg)

    return wall_per_step
    
def read_rung_dist(filename):
    """
    Reads the rungdistribution from ChaNGa stdout that has been saved to a disk
    by grepping the file for "Rung dist" and interpreting the output.
    
    Paramters
    ---------
    
    filename : str
        Filename for the ChaNGa stdout output
    
    Returns
    -------
    
    rung_dist : array
        Rung distribution for all big steps.  An array of integers of shape 
        (num big steps, num rungs)
    """
    filename = 'snapshot.out'
    command = "grep -i -e 'rung dist' " + filename
    output = subprocess.check_output(command, shell=True)
    output = output.strip('\n').split('\n')
    
    for i, line in enumerate(output):
        # Extract all the numbers from the rung distribution string
        a = line.split('(')[1].strip().split(',)')[0].split(', ')
        
        for j, b in enumerate(a):
            
            a[j] = int(b)
            
        # Save all numbers back to output
        output[i] = a
        
    rungdist = np.array(output, dtype=int)
    
    return rungdist