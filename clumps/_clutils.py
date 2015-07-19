# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 16:01:48 2015

@author: ibackus
"""
import numpy as np

def loadhalos(fname=''):
    """
     Load halo (.grp) file generated from fof
     Should be an ascii list of numbers, where the first row contains the
     total number of particles (gas+star+dark) and the remaining rows define
     which halo each particle belongs to
     """
    if fname == '':
        # Empty filename
        pass
    grp = np.loadtxt(fname,dtype=np.uint16)
    grp = grp[1:]   # (ignore the number of particles)

    return grp