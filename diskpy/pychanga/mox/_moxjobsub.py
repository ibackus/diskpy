#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Helper utilies for submitting jobs on mox (hyak 2) via the slurm workload
manager.

Created on Thu Apr 13 11:11:19 2017

@author: ibackus
"""

import os

from diskpy import global_settings
from diskpy.utils import configparser, configsave, units_from_param
from diskpy.pychanga import changa_command

def sbatch_script(workdir=None, param='snapshot.param', nodes=1, ppn=None, 
                  walltime=48, jobname='sbatch_job', email=None, 
                  changa_preset='default', partition='vsm', savepath=None,
                  changa_args='', runner_args='', **kwargs):
    """
    A not very robust function to generate SBATCH submission scripts for ChaNGa
    jobs on hyak.  Some of the requirements include:
    
    SBATCH (it's on mox akak hyak2, don't worry)
    
    By default, any directory with 'lastcheckpoint' in it will be restarted!
    If you don't want to just restart a simulation, delete lastcheckpoint from
    the simulation directory!
    
    Parameters
    ----------
    
    workdir : str
        (optional) Directory of the simulation.  Default is current working
        directory
    param : str
        Filename of the .param file (not full path)
    nodes : int
        Number of computation nodes to request
    ppn : int
        (optional) Number of cores per node
    walltime : float or int
        Walltime to request in hours
    jobname : str
        Name of the job
    email : str
        Email address for PBS to send updates to
    changa_preset : str
        ChaNGa preset to use (see diskpy.global_settings and 
        diskpy.pychanga.changa_command)
    partition : str
        mox partition (e.g. per group) to submit job to.
    savepath : str
        (optional) If supplied, save this script as an executable to this 
        path.
    changa_args, runner_args : str
        Optional extra command line arguments to pass to changa and the runner
        (i.e. charmrun, mpirun...)
    **kwargs
    flag pairs.  CAREFUL: these must not conflict with other flags.  Flag,val
    pairs (ie -flag val) are passed as: flag=val.  To add an extra command 
    line parameter to sbatch, e.g. -flag val
        
    Returns
    -------
    
    PBS_script : str
        A string for the PBS script.  Can be saved easily to file
        
    """
    if workdir is None:
        
        workdir = os.getcwd()
        
    # Setup filenames
    param_full = '$workdir/' + param
    outfile = '$workdir/' + 'changa.out'
    fprefix = os.path.splitext(param)[0]
    
    # Load the preset    
    if (changa_preset is None) or (changa_preset is 'default'):
        
        changa_preset = global_settings['changa_presets']['default']
        
    preset = global_settings['changa_presets'][changa_preset]
    
    # Set up the walltime for PBS
    hours = int(walltime)
    mins = int((walltime*60)%60)
    secs = int((walltime*3600)%60)
    walltime_str = '{0:d}:{1:02d}:{2:02d}'.format(hours,mins,secs)
    
    # Write the script!
    script = '#!/bin/bash -l\n'
    # Some sbatch flags
    script += \
        '#SBATCH -N {0}\n'.format(nodes) + \
        '#SBATCH -J {0}\n'.format(jobname) +\
        '#SBATCH -t {0}\n'.format(walltime_str) +\
        '#SBATCH -p {0}\n'.format(partition)
        
    if ppn is not None:
        
        script += '#SBATCH -n {0}\n'.format(nodes*ppn)
        
    if email is not None:        
        
        script += '#SBATCH --mail-type=ALL --mail-user={0}\n'.format(email)
    
    # Parse kwargs
    if kwargs is not None:
        
        for key, val in kwargs.iteritems():
            
            script += '#SBATCH -{0} {1}\n'.format(key, val)
            
    # Runtime initialization
    script += \
        'workdir={0}\n'.format(workdir) +\
        'outfile={0}\n'.format(outfile) +\
        'param={0}\n'.format(param_full) +\
        'cd $workdir\n' + \
        'changbin=$(which {0})\n'.format(preset[2])
    resumeCommand = changa_command('$param', preset=changa_preset,\
    restart_dir='{0}.chk$last'.format(fprefix), changa_bin='$changbin',
    changa_args=changa_args, runner_args=runner_args)
    changa_args += ' -wall {0}'.format(int(walltime*60))
    startCommand = changa_command('$param', preset=changa_preset, \
    changa_bin='$changbin', changa_args=changa_args, runner_args=runner_args)
    
    # Now assume that we want to restart if there is a checkpoint
    script += \
        'if [ -e "lastcheckpoint" ]\n' +\
        'then\n' +\
        '    echo "lastcheckpoint exists -- restarting simulation..."\n' +\
        '    last=`cat lastcheckpoint`\n' +\
        '    {0} >> $outfile 2>&1\n'.format(resumeCommand) +\
        'else\n' +\
        '    echo "lastcheckpoint doesnt exist -- starting new simulation..."\n' +\
        '    {0} &> $outfile \n'.format(startCommand) +\
        'fi\n'
    
    if savepath is not None:
        
        with open(savepath, 'w') as f:
            
            f.write(script)
            print "Submission script saved to:", savepath
            
        os.chmod(savepath, 0774)
        
    return script