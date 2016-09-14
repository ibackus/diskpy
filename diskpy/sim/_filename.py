# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:52:00 2016

@author: ibackus
"""
import os

class FileName():
    """
    A simple class for storing/handling filenames relative to base directory.
    
    Parameters
    ----------
    path : str or FileName
        Path to the file relative to the basedir.  path can optionally be an
        absolute path
        If a FileName, returns a shallow copy of path
    basedir : str or method
        Either the base directory the files are stored in or a method which 
        returns the based directory
    """
    def __init__(self, path, basedir=''):
        
        self.setdirectory(basedir)
            
        if isinstance(path, FileName):
            
            basedir = self.setdirectory(path._basedir)
            path = path.path
            
        elif (path is not None) and (os.path.isabs(path)):
            
            path = os.path.relpath(path, self.directory())
            
        self.path = path
        
    def __call__(self):
        
        return self.path
        
    def __repr__(self):
        
        return "<FileName: {0}>".format(self.path)
        
    def __str__(self):
        
        return str(self.path)
        
    def __lt__(self, y):
        
        if isinstance(y, FileName):
            
            y = y.abspath()
            
        return self.abspath().__lt__(y)
        
    def directory(self):
        """
        Returns the base directory that the file-path is relative to
        """
        return getdirectory(self._basedir)
        
    def setdirectory(self, basedir):
        """
        Set the base directory this file path is stored under
        
        Can be None, a str, or a method which returns a string
        """
        self._basedir = basedir
        
    def abspath(self):
        """
        Returns the absolute path of the file
        """
        directory = self.directory()
        fullpath = os.path.join(directory, self.path)
        return os.path.abspath(fullpath)
    
    def relpath(self, other):
        """
        Returns the paths of this relative to another directory
        """
        path = self.abspath()
        
        return os.path.relpath(path, other)
        
    def exists(self):
        
        if self.path is None:
            
            return False
            
        else:
            
            return os.path.exists(self.abspath())
        
    def ext(self):
        """
        Returns the file extension
        """
        return os.path.splitext(self.path)
        
def getdirectory(basedir):
    """
    Returns the directory from basedir which can be a string or a method which
    returns a string
    """
    try:
        
        return basedir()
        
    except TypeError:
        
        return str(basedir)