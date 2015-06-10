# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 17:35:36 2014

@author: ibackus
"""

import numpy as np

fname = 'hermite_spline_coeffs.dat'
f =open(fname,'r')

a_list = []
i = 0
order_list = []

for line in f:
    
    l = line.strip().split(',')
    order_list.append(int(l[0]))
    
    for n in range(len(l)):
        
        l[n] = float(l[n].strip())
        
    a_list.append(np.array(l[1:],dtype='float'))

order = np.array(order_list)
