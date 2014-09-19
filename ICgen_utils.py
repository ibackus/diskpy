# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 17:37:03 2014

@author: ibackus
"""
__version__ = "$Revision: 1 $"
# $Source$

# External modules
import subprocess
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import os

# ICgen modules
from ICglobal_settings import global_settings

class larray(list):
    """
    A simple subclass of list for containing listified arrays
    (see ICgen_utils.listify)
    Should be instantiated by ICgen_utils.listify
    
    USAGE:
    
    Creating an larray object:
    
        a = larray() # blank larray
        a = larray(shape) # sets a.shape = shape
    
    Then one can append to an larray just as with a list
    
        a.append(stuff)
        
    To return return a normal array:
    
        array = a.delistify()
    """
    
    def __init__(self, shape = None):
        
        list.__init__(self)
        self.shape = shape
    
    def delistify(self):
        """
        Return an array made from self.  See ICgen_utils.delistify
        """
        
        return delistify(self)

def listify(array, max_element=10**7):
    """
    Breaks up an array or SimArray into chunks and saves as an larray object
    (essentially, a list).  Useful for pickling very large arrays which otherwise
    may throw an error with pickle
    
    Whenever possible the array is NOT copied, rather a view is returned.  This
    depends on the functionality of array.ravel() (see numpy.ravel)
    
    **ARGUMENTS**
    
    array : array_like or SimArray
        Input array to listify
    max_element : int
        Maximimum number of elements per chunk (i.e., per list item)
        
    **RETURNS**
    
    list_array : larray
        A listified array.  
        
    See Also : ICgen_utils.delistify
    
    """
    
    # Initialize
    shape = array.shape
    n_elements = np.prod(shape)
    
    if n_elements <= max_element:
        
        return array
        
    else:
        
        out_list = larray(shape)
        array = array.ravel()
        
        # Number per slice
        N = int(max_element)
        
        counter = 0
        i = 0
        
        while counter < n_elements:
            
            out_list.append(array[counter:counter+N])
            i += 1
            counter += N
        
        return out_list
    
def delistify(in_list):
    """
    Reverses listify.
    Takes an larray object (a listified array) and returns the original array
    
    **ARGUMENTS**
    
    in_list : larray
        Listified array
        
    **RETURNS**
    
    array : array_like or SimArray
        array from the original data
        
    See Also : ICgen_utils.listify
    """
    
    if isinstance(in_list, larray):
        
        if isinstance(in_list[0], SimArray):
            
            array = SimArray(np.concatenate(in_list), in_list[0].units)
            
        else:
            
            array = np.concatenate(in_list)
        
        return array.reshape(in_list.shape)
        
    else:
        
        return in_list
        
def est_eps(smoothlength_file):
    """
    Estimate gravitational softening length (eps) from a ChaNGa output .smoothlength
    file.  eps is estimated as 1/2 the mean smoothing length
    
    **ARGUENTS**
    
    smoothlength_file : str
        Filename of the .smoothlength file
        
    **RETURNS**
    
    eps : float
        Estimate of the gravitational softening length in simulation units
    """
    # Open ChaNGa output file containing smoothing lengths for all particles
    f = open(smoothlength_file, 'r')
    # Total number of particles (incl. star)
    nParticles = int(f.readline().strip())
    # Allocate smoothing length array    
    smoothlength = np.zeros(nParticles, dtype=np.float32)
    # Read smoothing lengths
    for i, line in enumerate(f):
        
        smoothlength[i] = float(line.strip())
        
    # Calculate eps, ignoring star particle
    mean_smooth = (smoothlength.sum() - smoothlength[-1])/(nParticles-1)
    eps = mean_smooth/2
    
    return eps   

def changa_run(command, verbose = True, logfile_name=None, force_wait=False):
    """
    A wrapper for running ChaNGa
    
    **ARGUMENTS**
    
    command : str
        A full command line command for running ChaNGa.  Can be produced from 
        defaults using ICgen_utils.changa_command
    verbose : bool
        (optional) Flag for printing ChaNGa output to stdout.
        If True - stdout is printed.  This will effectively makes changa_run
        wait on ChaNGa completion
    logfile_name : str
        (optional) If set, saves ChaNGa output to file
    force_wait : bool
        (optional) Default = False
        If set, forces wait on ChaNGa before completion
    
    **RETURNS**
    
    p : subprocess.Popen
        A process object created by subprocess.Popen for the ChaNGa command
    """
    
    if logfile_name is not None:
        
        logfile = open(logfile_name, 'w')
        logfile.close()
        logfile = open(logfile_name, 'a')
    
    if verbose:
        
        output = subprocess.PIPE
        p = subprocess.Popen(command.split(), stderr=output, stdout=output)
        
        for line in iter(p.stdout.readline, ''):
            
            print line,
            if logfile_name is not None:
                
                logfile.write(line)
                
        p.wait()
        
    else:
        
        if logfile_name is not None:
            
            output = logfile
            
        else:
            
            output = subprocess.PIPE
            
        p = subprocess.Popen(command.split(), stderr=output, stdout=output)
        
    if force_wait:
        
        p.wait()
        
    return p
    
def changa_command(param_name, preset=None, changa_bin=None, changa_args='', runner_args=''):
    """
    A utility for created command line commands for running ChaNGa
    
    **ARGUMENTS**
    
    param_name : str
        Filename of the .param file used for ChaNGa
    preset : str
        if None, the default preset is used
        Presets are defined in global_settings
    changa_bin : str
        Path to the ChaNGa binary to use.  If None, defaults are used
        Overrides preset binary
    changa_args : str
        Additional user supplied arguments for ChaNGa
    runner_args : str
        Additional user supplied arguments for the runner (ie charmrun or mpirun)
        
    **RETURNS**
    
    command : str
        A command line command for running ChaNGa
    """
    
    # Contains all the presets
    preset_dict = global_settings['changa_presets']
    
    # Load the preset    
    if preset is None:
        
        preset = preset_dict['default']
        
    preset_list = preset_dict[preset]
    
    # Get full path to ChaNGa binary
    if changa_bin is None:
        
        changa_bin = preset_list[2]
    
    changa_bin = os.popen('which ' + changa_bin).read().strip()
    
    if '' == changa_bin:
        
        raise RuntimeError, 'Could not find ChaNGa.  Try different preset'
    
    # Merge user defined extra arguments    
    runner_args = ' '.join([preset_list[1], runner_args])
    changa_args = ' '.join([preset_list[3], changa_args])
    runner = preset_list[0]
    
    command = ' '.join([runner, runner_args, changa_bin, changa_args, param_name])
    command = ' '.join(command.split())
    
    return command            
            
def arg_cat(arg_list):
    """
    STILL ALPHA!!! 
    
    arg_str = arg_cat([args1, args2, ...])
    
    Concatenates a list of various command line arguments.  arg_list should
    be a list containing command line argument strings.
    
    Priority is given to later arguments.  So arg_list[2] overwrites arg_list[1]
    if they share any flags
    
    **EXAMPLES**
    
    args1 = '-n 20 +cd 13 ++b fire testit'
    args2 = '-n 20 +cd 300'
    print arg_cat([args1, args2])
    
        returns:
        +cd 300 -n 20 testit ++b fire
    """
    
    args_dict = {}
    
    # Loop through all sets of arguments
    for args in arg_list:
        
        # Split args if it's not already a list
        if isinstance((args), str):
            
            args = args.split()
            
        # Parse arguments
        counter = 0
        
        while counter < len(args):
            
            key = args[counter]
            
            if (key[0] == '-') or (key[0] == '+'):
                # We have a flag followed by its value
                val = args[counter+1]
                counter += 2
                
            else:
                # We just have an argument
                val = ''
                counter += 1
                
            args_dict[key] = val
                
    args_str = ''
    
    for key, val in args_dict.iteritems():
        
        if val == '':
            
            # Tack on to the end
            args_str = ' '.join([args_str, key])
            
        else:
            
            # Place at beginning
            args_str = ' '.join([key, val, args_str])
            
    args_str = ' '.join(args_str.split())
    
    return args_str