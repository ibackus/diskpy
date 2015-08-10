# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:23:03 2015

@author: ibackus
"""

from diskpy.utils import get_module_names as _getnames

__all__ = _getnames(__file__)

for modulename in __all__:
    
    exec 'import {0}'.format(modulename)