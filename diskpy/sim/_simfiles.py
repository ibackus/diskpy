# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 16:10:47 2016

@author: ibackus
"""
import os

import diskpy
from diskpy.pychanga import getpar

from diskpy.sim import FileName, FileList
from diskpy.sim._general import loadParam, inferLogPath, \
loadLog, listOutputs, findParams



class SimFiles():
    """
    Parameters
    ----------
    paramname : str
    directory : str
    logfile : str
    """
    def __init__(self, directory, paramname=None, logfile=None, 
                 notefile='notes.txt'):
        
        self.verbose = True
        self.files = {}
        self.files['notes'] = self.asfilename(notefile)
        self.set_directory(directory)
        # load the params
        self.files['param'] = self.asfilename(paramname)
        self.files['log'] = self.asfilename(logfile)
        self.updateparam()
        # get the filenames
        self.scanfolder()
        
    def __repr__(self):
        
        string = '<SimFiles>\n'
        string += 'directory: {0}\n'.format(self.directory())
        string += 'param file: {0}\n'.format(self.files['param'])
        string += 'log file: {0}\n'.format(self.files['log'])
        string += 'Inital Conditions: {0}\n'.format(self.files['IC'])
        string += 'Time steps: {0}\n'.format(len(self))
        
        return string
        
    def __len__(self):
        
        return len(self.files['outputs'])
        
    def asfilename(self, path):
        """
        Casts a path (str, None, or FileName) as a FileName object stored in
        the simulation's directory
        """
        
        if isinstance(path, FileName):
            
            return path
            
        else:
            
            return FileName(path, self.directory)
        
    def updateparam(self, newparam=None, newlog=None):
        """
        Re-load the runtime parameters used for the simulation, possibly using
        new files
        
        If no param file exists, .param files will be searched for
        
        Parameters
        ----------
        newparam : str or FileName
            Path to the new .param file
        newlog : str or FileName
            Path to the new .log file
        """
        # Cast as filenames
        newparam = self.asfilename(newparam)
        newlog = self.asfilename(newlog)
        
        if newparam.path is not None:
            
            self.files['param'] = newparam
            
        if newlog.path is not None:
            
            self.files['log'] = newparam
            
        paramname = self.files['param']
        logfile = self.files['log']
        # Try to find param if required
        if (paramname.path is None) and (logfile.path is None):
            
            paramname = self.asfilename(findParams(self.directory(), self.verbose))
            
            if paramname.path is None:
                
                raise RuntimeError, "could not find a .param file in {0}"\
                .format(self.directory())
        
        param = loadParam(paramname)
        logfile = inferLogPath(paramname, logfile)
        log = loadLog(logfile)
        usedparam = {}
        usedparam.update(param)
        usedparam.update(log)
        
        self.params = {'param': param, 'log': log, 'usedparam': usedparam}
        self.files['log'] = logfile
        self.files['param'] = paramname
        
    def parseparam(self, paramname=None):
        
        if paramname is None:
            
            paramname = self.files['param']
            
        
    def scanfolder(self):
        
        outputs = listOutputs(self.directory(), self.params['usedparam'])
        self.files['outputs'] = FileList(outputs, self.directory)
        IC = self.asfilename(self.params['usedparam']['achInFile'])
        if not os.path.exists(IC.abspath()):
            
            IC = self.asfilename(None)
            
        self.files['IC'] = IC
        
    def printNoteFile(self):
        
        notefile = self.files['notes']
        
        if notefile.exists():
            
            print '<Notes stored in {0}>'.format(notefile.path)
            
            with open(notefile.abspath(), 'r') as f:
                
                print f.read()
        
        
    def directory(self):
        
        return self._directory
        
    def set_directory(self, directory):
        """
        Changes the directory the sim is assumed to be in.  Useful if, for
        instance, you move the simulation
        """
        directory = os.path.abspath(directory)
        self._directory = directory
        