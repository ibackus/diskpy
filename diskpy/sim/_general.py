# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:58:16 2016

@author: ibackus
"""
import os
import re
import glob
from diskpy.utils import configparser, logparser
from diskpy.pychanga import getpar, get_fnames
from diskpy.sim import FileName, FileList

def findParams(directory=None, verbose=True):
    """
    Scans a directory to find a .param file for a simulation.  The most
    recently changed param file is taken as the simulation param file
    """
    if directory is None:
        
        directory = os.getcwd()
        
    searchstr = os.path.join(directory, '*.param')
    files = [f for f in glob.glob(searchstr) if os.path.isfile(f)]
    files.sort(key=os.path.getmtime)
    
    if len(files) == 0:
        
        return None
        
    if verbose and len(files) > 1:
        
        print "Multiple .param files found in {0}.  Taking most recent file"\
        .format(directory)
        
    param = os.path.relpath(files[-1], directory)
    return param
    
def listOutputs(directory=None, params={}):
    """
    Finds the output snapshots for a simulation run with the params in a 
    given directory.  Uses diskpy.pychanga.get_fnames to search the directory
    
    This returns a list of the paths relative to the directory
    """
    if directory is None:
        
        directory = os.getcwd()
        
    fnames = []
    
    fprefix = getpar('achOutName', params)
    fprefix = os.path.join(directory, fprefix)
        
    fnames = get_fnames(fprefix, directory)
    
    for i, fname in enumerate(fnames):
        
        fnames[i] = os.path.relpath(fname, directory)
    
    return fnames
    
def inferLogPath(paramname=None, logfile=None):
    """
    Tries to get the log path (relative to the base simulation directory) from
    either the paramfile or (if supplied) the logfile
    
    paramname and logfile should be relative to directory
    
    Returns a FileName
    """
    # Cast as FileNames
    paramname = FileName(paramname)
    logfile = FileName(logfile)
    # check
    if (logfile.path is None) and (paramname.path is None):
        
        raise ValueError, "No param names or log files supplied"
    
    # try to guess the logfile name
    params = loadParam(paramname)
    
    if logfile.path is None:
        
        logfilename = getpar('achOutName', params) + '.log'
        logfile = FileName(logfilename, paramname.directory)
        
    return logfile
    
def loadParam(paramname):
    """
    Attempts to load the param file stored in directory.  If paramname is None
    an empty dict is returned
    """
    param = {}
    
    if paramname.path is not None:
        
        param.update(configparser(paramname.abspath(), 'param'))
        
    return param
    
def loadLog(logfile):
    """
    Attempts to load the logfile stored in directory.  If logfile is None, 
    an empty dict is returned
    """
    log = {}
    # Load the log file
    if logfile.path is not None:
        
        if os.path.exists(logfile.abspath()):
            
            log.update(logparser(logfile.abspath()))
            
    return log