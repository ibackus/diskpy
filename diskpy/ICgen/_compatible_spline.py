# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:51:09 2016

@author: ibackus
"""

import scipy
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from warnings import warn
import numpy as np

from distutils.version import LooseVersion as version

# InterpolatedUnivariateSpline in older versions of scipy don't handle regions
# outside the bbox and don't have the argument ext=0
if version(scipy.__version__) < '0.16.0':
    warn('WARNING: Deprecated version of scipy.  Please upgrade to at least 0.16.0')
    
    spl = spline
    
    # OVERRIDE spline if scipy is outdated
    def spline(x, y, w=None, bbox=[None]*2, k=3, ext=0, check_finite=False):
        
        splobj = spl(x, y, w, bbox, k)
        bbox = list(splobj._data[3:5])
        
        if ext == 'zeros':
            
            ext = 0
        
        def splfcn(xpts):
            
            xpts = np.asarray(xpts)
            yinterp = splobj(xpts)
            mask = (xpts < bbox[0]) | (xpts > bbox[1])
            if np.any(mask):
                yinterp[mask] = ext
            
            return yinterp
            
        return splfcn
