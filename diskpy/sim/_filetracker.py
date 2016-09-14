# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:22:07 2016

@author: ibackus
"""
import numpy as np
import pynbody
from pynbody.snapshot import SimSnap
from diskpy.sim import FileName, FileList
import diskpy

class FileTracker():
    
    def __init__(self):
        
        self.filenames = []
        self.isloaded = []
        self.snapshots = []
        
    def __len__(self):
        
        return len(self.filenames)
        
    def addFile(self, fname):
                
        if fname not in self.filenames:
            
            self.filenames.append(fname)
            self.isloaded.append(False)
            
    def markLoaded(self, fname):
        
        i = self.filenames.index(fname)
        self.isloaded[i] = True
    
    def numLoaded(self): 
        
        return np.sum(self.isloaded)
        
    def deleteFile(self, fname):
        
        i = self.filenames.index(fname)
        self.filenames.pop(i)


def overrideLoad(f, filetracker, filename=None):
    
    f._old_load_array = f._load_array
    if filename is None:
        
        filename = f.filename
        
    filetracker.addFile(filename)
    
    def load(*args, **kwargs):
        
        output = f._old_load_array(*args, **kwargs)
        filetracker.markLoaded(filename)
        
    f._load_array = load
    return