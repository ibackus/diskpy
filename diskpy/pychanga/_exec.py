# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:04:27 2016

@author: ibackus
"""
import numpy as np
import re
import os
import subprocess
from subprocess import Popen, PIPE

from diskpy import global_settings

def est_time_step(param_name, preset='default', dDelta0=100, changa_args='', runner_args=''):
    """
    A routine to automatically estimate a reasonable time-step size for ChaNGa.
    The idea is to have about half the particles fall into the lowest rung (ie 
    the big time step).  This is done by calculating the rung distribution for 
    a large time step by running ChaNGa and killing ChaNGa once it has output 
    rung distribution.
    
    NOTE: this is still fairly alpha.  A better version would probably not
    just run ChaNGa and then try to kill it.  To be safe, a local ChaNGa preset
    should probably be used.
    
    **ARGUMENTS**
    
    param_name : str
        Filename for a ChaNGa .param file which defines parameters for the
        snapshot.  The snapshot must already be saved to disk
    preset : str
        changa_runner preset to use.  See diskpy.global_settings
    dDelta0 : int or float
        Some large time step that should place all the particles at higher
        rungs.
    changa_args : str
        Additional command line arguments to pass to changa.  CANNOT include
        -n (number of time steps) or -dt (timestep size)
    runner_args : str
        Additional command line arguments to pass to the runner, ie to 
        charmrun or mpirun
        
    **RETURNS**
    
    dDelta : float
        Estimated reasonable time step that places half the particles on the
        lowest rung (ie the big time step)
    """
    
    settings = global_settings['changa_presets'][preset]
    changa_name = settings[2]
    runner_name = settings[0]
    
    changa_args += ' -n 1 -dt {0}'.format(dDelta0)
    command = changa_command(param_name, preset, changa_args=changa_args, runner_args=runner_args)
    
    rung_line = ''
    p = changa_run(command, verbose=False)
    
    for line in iter(p.stdout.readline, ''):
        
        if 'rung distribution' in line.lower():
            
            # Kill the runner
            kill_command = 'pkill -9 ' + runner_name
            pkill = subprocess.Popen(kill_command.split(), \
            stdout=subprocess.PIPE)
            pkill.wait()
            
            # Kill ChaNGa
            kill_command = 'pkill -9 ' + changa_name
            pkill = subprocess.Popen(kill_command.split(), \
            stdout=subprocess.PIPE)
            pkill.wait()
            
            rung_line = line.strip()
            break
        
    if rung_line == '':
        
        raise RuntimeError('ChaNGa failed to output rung distribution')
        
    rung_list = re.findall('\d+', rung_line)
    rung_hist = np.array(rung_list).astype(float)
    rung_edges = np.arange(len(rung_hist) + 1, dtype=float)
    
    s = np.cumsum(rung_hist)
    Ntot = s[-1]
    
    # Find first bin which gives us more than half the total number
    for i, n in enumerate(s):
        
        if n >= 0.5*Ntot:
            
            ind = i
            break
    
    # Calculate the median rung    
    rung_med = rung_edges[ind] + (0.5*Ntot - s[ind-1])/rung_hist[ind]
    
    # Now estimate a time step that will fit about half the particles on the
    # lowest rung (ie the big time step)
    
    dDelta = dDelta0 * 2.0**(-rung_med+1)
    
    return dDelta

def changa_run(command, verbose = True, force_wait=False, return_success=False):
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
    force_wait : bool
        (optional) Default = False
        If set, forces wait on ChaNGa before completion
    return_success : bool
        If set, returns the success of the run, i.e. checks to see if the 
        simulation finished properly.
    
    **RETURNS**
    
    p : subprocess.Popen
        A process object created by subprocess.Popen for the ChaNGa command
    success : bool
        IF return_success is set to True, the status of the run is returned
        instead of p.
    status : bool
    """
    output = subprocess.PIPE
    p = subprocess.Popen(command.split(), stderr=output, stdout=output)
    success = False
    
    if verbose or return_success:
        
        for line in iter(p.stdout.readline, ''):
            
            if "Done." in line:
                
                success = True
                
            if verbose:
                
                print line,
                
        p.wait()
        
    if force_wait:
        
        p.wait()
        
    
    if return_success:
        
        return success
        
    else:
        
        return p

def changa_command(param_name, preset='default', changa_bin=None, \
changa_args='', runner_args='',restart_dir=None):
    """
    A utility for created command line commands for running ChaNGa
    
    **ARGUMENTS**
    
    param_name : str
        Filename of the .param file used for ChaNGa
    preset : str
        if None or 'default', the default preset is used
        Presets are defined in global_settings
    changa_bin : str
        Path to the ChaNGa binary to use.  If None, defaults are used
        Overrides preset binary
    changa_args : str
        Additional user supplied arguments for ChaNGa
    runner_args : str
        Additional user supplied arguments for the runner (ie charmrun or mpirun)
    restart_dir : str
        (optional) If set, this will be treated as a restart.  All changa args
        will be ignored.
        
    **RETURNS**
    
    command : str
        A command line command for running ChaNGa
    """
    
    # Contains all the presets
    preset_dict = global_settings['changa_presets']
    
    # Load the preset    
    if (preset is None) or (preset is 'default'):
        
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
    
    if restart_dir is not None:
        
        changa_args = '+restart {0}'.format(restart_dir)
        
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