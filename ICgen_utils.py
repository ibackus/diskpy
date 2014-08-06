# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 17:37:03 2014

@author: ibackus
"""

import subprocess
import multiprocessing
import os

def changa_run(param_name, preset='local', changa_args=None, runner_args=None, changa_bin=None):
        
        if preset == 'local':
            
            # Set up defaults
            
            # Number of processes to use
            proc = multiprocessing.cpu_count() - 1
            
            # location of the ChaNGa binary
            if changa_bin is None:
                
                changa_bin = os.popen('which ChaNGa').read().strip()
            
            if proc < 1:
                
                proc = 1
                
            # Now parse added arguments and compare to defaults
            
            'charmrun {} +p {} ++local {} {} {}'.format()
            
def arg_cat(arg_list):
    """
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
        
        args_str = ' '.join([args_str, key, val])
            
    args_str = ' '.join(args_str.split())
    
    return args_str