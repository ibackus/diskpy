# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 23:25:31 2015

@author: ibackus
"""

import os
import glob
import re
from subprocess import Popen, PIPE
import pynbody

from diskpy import global_settings
from diskpy.utils import configparser, configsave, units_from_param

from .. import changa_command

def pbs_script(workdir=None, param='snapshot.param', nodes=1, ppn=12, walltime=48, \
jobname='PBS_job', backfill=False, email=None, changa_preset='default', **kwargs):
    """
    A not very robust function to generate PBS submission scripts for ChaNGa
    jobs on hyak.  Some of the requirements include:
    
    mpi version of ChaNGa
    gcc_4.4.7-ompi_1.6.5
    PBS (it's on hyak, don't worry)
    
    By default, any directory with 'lastcheckpoint' in it will be restarted!
    If you don't want to just restart a simulation, delete lastcheckpoint from
    the simulation directory!
    
    Parameters
    ----------
    
    *Required*
    
    workdir : str
        (optional) Directory of the simulation.  Default is current working
        directory
        
    *Optional*
    
    param : str
        Filename of the .param file (not full path)
    nodes : int
        Number of computation nodes to request
    ppn : int
        Number of cores per node
    walltime : float or int
        Walltime to request in hours
    jobname : str
    backfill : bool
        Boolean flag for whether to use the backfill (default is TRUE)
    email : str
        Email address for PBS to send updates to
    changa_preset : str
        ChaNGa preset to use (see diskpy.global_settings and 
        diskpy.ICgen.ICgen_utils.changa_command)
        
    **kwargs
    flag pairs.  CAREFUL: these must not conflict with other flags.  Flag,val
    pairs (ie -flag val) are passed as: flag=val
        
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
        
    preset = global_settings['changa_presets'][changa_preset]
    
    # Set up the walltime for PBS
    hours = int(walltime)
    mins = int((walltime*60)%60)
    secs = int((walltime*3600)%60)
    walltime_str = '{0:d}:{1:02d}:{2:02d}'.format(hours,mins,secs)
    
    # Write the script!
    
    # Start with the shebang
    script = '#!/bin/bash\n'
    
    # Some PBS flags
    script += '#PBS -N {0}\n\
#PBS -j oe\n\
#PBS -l nodes={1}:ppn={2},feature={2}core\n\
#PBS -l walltime={3}\n\
#PBS -V\n'.format(jobname, nodes, ppn, walltime_str)
    
    # Email stuff
    if email is not None:
        
        script += '#PBS -M {0}\n'.format(email)
        script += '#PBS -m be\n'
        
    # Choose whether to use the backfill
    if backfill:
        
        script += '#PBS -q bf\n'
        
    # Parse kwargs
    if kwargs is not None:
        
        for key, val in kwargs.iteritems():
            
            script += '#PBS -{0} {1}\n'.format(key, val)
            
    # Runtime initialization
    script += 'module load gcc_4.4.7-ompi_1.6.5\n\
export MX_RCACHE=0\n\
workdir={0}\n\
cd $workdir\n\
changbin=$(which {1})\n'.format(workdir, preset[2])

    resumeCommand = changa_command(param_full, preset=changa_preset,\
    changa_args='+restart {0}.chk$last -wall {1}'.format(fprefix, int(walltime*60)))
    startCommand = changa_command(param_full, preset=changa_preset, \
    changa_args='-wall {0}'.format(int(walltime*60)))
    
    # Now assume that we want to restart if there is a checkpoint
    script += 'if [ -e "lastcheckpoint" ]\n\
then\n\
    echo "lastcheckpoint exists -- restarting simulation..."\n\
    last=`cat lastcheckpoint`\n\
    {0} >> {1} 2>&1\n'.format(resumeCommand, outfile)
    script += 'else\n\
    echo "lastcheckpoint doesnt exist -- starting new simulation..."\n\
    {0} >& {1}\n\
fi\n'.format(startCommand, outfile)
    
    return script

def subber(directories, scriptname, scriptstr=None, subfile=None):
    """
    Submits scriptname, contained in all the directories, to the submission
    queue using qsub.  Optionally, the script can be provided as a string (or
    a list of strings for each simulation) in the variable scriptstr
    
    Optionally a bash script can be saved to subfile instead of submitting
    the scripts.  This can be useful for batch submission when using a
    computation node.
    
    Note: the submission scripts will be made executable for all users
    
    
    **ARGUMENTS**
    
    directories : str, list, ndarray, tuple
        The directory or directories containing the submission script
    scriptname : str
        Script to submit.  Should be present in all the directories
    scriptstr : str, list, or None
        Default = None (do nothing!)
        Either a string containing the PBS submission script (see PBS_script)
        or a list of such strings.  These will be saved to directories/scriptname
    subfile : str or None
        (optional) If a string, filename for a bash script to be saved instead 
        of executing the qsub commands
    """
    
    # Make the directories list iterable
    if isinstance(directories, str):
        
        directories = [directories]
        
    # Make the submission scripts an iterable list
    if isinstance(scriptstr, str):
        
        scriptstr = [scriptstr] * len(directories)
    
    # Get current working directory
    cwd = os.getcwd()
    
    # Change directories to full paths
    fullpaths = []    
    for directory in directories:
        
        fullpaths.append(os.path.abspath(directory))
     
    # Submission command
    command = 'qsub ' + scriptname
    
    if subfile is not None:
        
        bashscript = open(subfile, 'w')
        bashscript.write('#/bin/bash\n')
    
    # Submit all the scripts
    for i, fullpath in enumerate(fullpaths):
        
        os.chdir(fullpath)
        
        # If submission scripts have been provided as strings, write them out
        if scriptstr is not None:
            
            f = open(scriptname, 'w')
            f.write(scriptstr[i])
            f.close()
            
        # Make them executable
        os.system('chmod a+x ' + scriptname)
            
        if subfile is not None:
            
            bashscript.write(command + '\n')
            
        else:
            
            # Submit the script to PBS
            p = Popen(command.split(), stderr=PIPE, stdout=PIPE)
            for line in iter(p.stdout.readline, ''):
                
                print line,
    
            p.wait()
        
    # Finalize
    if subfile is not None:
        
        bashscript.close()
        os.system('chmod a+x ' + subfile)
        
    os.chdir(cwd)

def make_continue_sub(simdir='.', paramname='snapshot.param', \
newparam='continue.param', t=None, t_extra=None, oldsub='subber.sh', \
newsub='cont_subber.sh'):
    """
    Makes a submission script for continuing a simulation from a previous output.
    Also makes a .param file for the continued simulation.  The simulation
    will be continued in the same directory, with snapshot numbering scheme
    for outputs being the same.
    
    Parameters for the original simulation cannot be altered (except the number
    of steps you want to continue the simulation by).  PBS runtime parameters
    also cannot be changed (number of nodes, walltime, etc...)
    
    Any checkpoints will be deleted.
    
    Requires a submission script be present for the original simulation
    
    NOTE: if nSteps, nSteps_extra are not set, the total number of steps
    to simulate is not changed.
    
    
    **ARGUMENTS**
    
    simdir : str
        The simulation directory
    paramname : str
        Filename of the .param file for the simulation
    newparam : str
        filename for the .param file for the continued simulation
    t : float or SimArray
        Total simulation time to run.
        If no units are specified, it is in simulation units
    t_extra : float or SimArray
        Extra simulation time to run
        If no units are specified, it is in simulation units
        OVERIDES t!!!
    oldsub : str
        Filename for the original submission script
    newsub : str
        Filename for the new submission script
        
    **RETURNS**
    
    sub_path : str
        Full path to the PBS submission script
    
    """
    
    # Lazy man's way of dealing with files in another directory
    cwd = os.getcwd()
    
    os.chdir(simdir)
    
    # Load param file
    param = configparser(paramname, 'param')
    fprefix = param['achOutName']
    
    # Find all the outputs.  They should be of the format fprefix.000000
    search_exp = '^' + fprefix + '.(?:(?<!\d)\d{6}(?!\d))$'
    flist = []
    
    for fname in glob.glob(fprefix + '*'):
        
        if re.search(search_exp, fname) is not None:
            
            flist.append(fname)
    
    # Find the number of the last output (the last 6 chars should be an int)
    flist.sort()
    iStartStep = int(flist[-1][-6:])
    param['iStartStep'] = iStartStep
    param['achInFile'] = flist[-1]    
    dDelta = param['dDelta']
    
    # Set the number of steps to run        
    if t_extra is not None:
        
        # Convert to simulation units if needed
        if pynbody.units.has_units(t_extra):
            
            t_unit = units_from_param(param)['t_unit']
            t_extra.convert_units(t_unit)
        
        # Assign output
        param['nSteps'] = iStartStep + int(round(t_extra/dDelta))
        
    elif t is not None:
        
        # Convert to simulation units if needed
        if pynbody.units.has_units(t):
            
            t_unit = units_from_param(param)['t_unit']
            t.convert_units(t_unit)
            
        # Assign output
        param['nSteps'] = int(round(t/dDelta))
    
    # Save new param file
    configsave(param, newparam, ftype='param')
    
    # Delete old checkpoints

    for checkpoint in glob.glob(fprefix + '.chk*'):
        
        print 'removing ' + checkpoint
        os.system('rm -rf ' + checkpoint)
        
    if os.path.exists('lastcheckpoint'):
        
        print 'removing lastcheckpoint'
        os.remove('lastcheckpoint')
    
    # Create a submission script for the simulation continuation
    oldsubfile = open(oldsub, 'r')
    newsubfile = open(newsub, 'w')
    
    for line in oldsubfile:
        
        newsubfile.write(line.replace(paramname, newparam))
        
    oldsubfile.close()
    newsubfile.close()
    
    # Make the submission script executable
    os.chmod(newsub, 0777)
    
    sub_path = os.path.abspath(newsub)
    
    # Change back to original working directory
    os.chdir(cwd)
    
    return sub_path