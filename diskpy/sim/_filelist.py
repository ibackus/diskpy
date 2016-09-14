# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:54:09 2016

@author: ibackus
"""

import os
from diskpy.utils import CustomList
from diskpy.sim import FileName

class FileList(CustomList):
    """
    FileList is a class for storing a list of filenames relative to a base
    directory.  FileList objects have access much like python lists, but with
    added functionality.  File names are supplied as either FileName objects
    or as strings and are stored as FileName objects.
    
    Parameters
    ----------
    fnames : list or iterable
        A list or iterable containing filenames as strings or as FileName 
        objects
    basedir : str or method
        Either the base directory the files are stored in or a method which 
        returns the based directory
    """
    def __init__(self, fnames=[], basedir=None):
        
        if basedir is None:
            
            basedir = os.getcwd()
            
        self._basedir = basedir
        CustomList.__init__(self, fnames)
        self._className = "FileList"
        
    def __getitem__(self, ind):
        
        if isinstance(ind, slice):
            
            indices = ind.indices(len(self))
            return FileList([self._getOneItem(i) for i in range(*indices)])
            
        else:
            
            return self._getOneItem(ind)
            
    def _setOneItem(self, ind, filename):
        
        if isinstance (filename, str):
            # Cast it as a FileName obj
            filename = FileName(filename, self.directory)
        
        CustomList._setOneItem(self, ind, filename)
    
    def setdirectory(self, basedir):
        """
        Set the base directory this file path is stored under
        
        Can be None, a str, or a method which returns a string
        """
        self._basedir = basedir
        
    def directory(self):
        """
        Returns the directory this file list is stored in
        """
        return getdirectory(self._basedir)
            
    def path(self):
        """
        Returns a list of the file paths relative to the base directory
        """
        return [fname.path for fname in self]
        
    def abspath(self):
        """
        Returns the absolute path this is stored in
        """
        return [fname.abspath() for fname in self]
        
    def relpath(self, other):
        """
        Returns the paths of this relative to another directory
        """
        return [fname.relpath(other) for fname in self]
        
    def ext(self):
        """
        Returns the file extensions
        """
        return [fname.ext() for fname in self]
        
def getdirectory(basedir):
    """
    Returns the directory from basedir which can be a string or a method which
    returns a string
    """
    try:
        
        return basedir()
        
    except TypeError:
        
        return str(basedir)