# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 14:22:07 2016

@author: ibackus
"""

import pynbody
from pynbody.snapshot import SimSnap

class FileTracker():
    
    def __init__(self):
        
        self.files = []
        self.isloaded = []
        
    def __len__(self):
        
        return len(self.files)
        
    def addFile(self, fname):
        
        if fname not in self.files:
            
            self.files.append(fname)
            self.isloaded.append(False)
            
    def markLoaded(self, fname):
        
        i = self.files.index(fname)
        self.isloaded[i] = True
        
    

def overrideLoad(f, filetracker, filename=None):
    
    f._old_load_array = f._load_array
    fname = f.filename
    filetracker.addFile(f.filename)
    
    def load(*args, **kwargs):
        
        output = f._old_load_array(*args, **kwargs)
        filetracker.markLoaded(fname)
        
    f._load_array = load
    return