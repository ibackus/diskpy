# -*- coding: utf-8 -*-
"""
Calculates cubic spline interpolations for sigma(r) and probability(r)
probability = 2*pi*r*sigma

Created on Mon Jan 27 13:00:52 2014

@author: ibackus
"""
import pynbody
import numpy as np
import cPickle as pickle
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad

import isaac

def sigma(fName=None,Rd=1.0,r_in=0.05,Md=1.0,cutlength=1.0,kind='power'):
    """
    ****************
    
    By default, generates a spline interpolation of sigma vs r
    Returns sigma vs r as an object whose call method
    utilizes cubic spline interpolation (see scipy.interpolation.interp1d)
    
    fName should contain a pickled dictionary with the entries:
        'sigma': <sigma evaluated at r>
        'r':     <r for the bins>
        
    If the input sigma has units, sigma vs r will be returned in units of
    Msol/au^2
    
    *****************    
    
    If fName=None, user must define a function to calculate sigma(r)
            
    """
    if fName == None:
        # USER MUST DEFINE A FUNCTION TO CALCULATE SIGMA(R)
        def sigout(r):
            # Output is returned  
            return output
    else:
        print 'Loading {0}'.format(fName)
        inDict  = pickle.load(open(fName,'rb'))
        sigmaBinned = inDict['sigma']
        rBinned = inDict['r']
        if pynbody.units.has_units(sigmaBinned):
            sigmaBinned.convert_units('Msol au**-2')
        print 'Calculating spline interpolation (slow for many data points)'
        sigout = interp1d(rBinned,sigmaBinned,kind='cubic',fill_value=0.0,\
        bounds_error=False)
    return sigout
    
def prob(fName=None):
    """
    By default, returns un-normalized probability vs r as a function whose call method
    utilizes cubic spline interpolation (see scipy.interpolation.interp1d)
    
    fName should contain a pickled dictionary with the entries:
        'sigma': <sigma evaluated at r>
        'r':     <r for the bins>   
    
    probability is calculated as sigma*2*pi*r
    
    If fName=None, user must define a function to calculate prob(r)
    """
    if fName == None:
        def pr(r):
            # Enter fuction here.  output is returned
            return output
    else:
        inDict  = pickle.load(open(fName,'rb'))
        rBinned = inDict['r']
        prBinned = inDict['sigma']*(2*np.pi*rBinned)
        print 'Calculating spline interpolation (slow for many data points)'
        pr = interp1d(rBinned,prBinned,kind='cubic',fill_value=0.0,\
        bounds_error=False)        
        #pr = UnivariateSpline(rBinned,prBinned,k=3,s=0)
    return pr

def cdfinv_r(fName=None,pdf=None):
    """
    Calculates the inverse of the cumulative distribution function for
    probability as a function of r.
    
    *** Arguments ***
    
    * fName *   File name for retrieving sigma(r).  If None, the user must
        define a function to calculate the inverse CDF.  Otherwise,
        fName should contain a pickled dictionary with the entries:
            'sigma': <sigma evaluated at r>
            'r':     <r for the bins>   
            
    * pdf *     Function to calculate pdf as a function of radius.  If None,
        pdf is calculated using the function prob() (defined in calc_sigma.py)
        Otherwise, pdf should be callable, for instance as a spline interp.
        
    
    *** Returns ***
        
    Returns a spline interpolation of the inverse CDF.  This is normalized
    such that the spline interpolation is in [0,1].
    """
    if fName is None:
        def finv(r):
            # Define function here
            return output
    else:
        print 'calculating CDF'
        # Calculate the CDF from prob
        inDict = pickle.load(open(fName,'rb'))
        r = inDict['r']
        r[0] = 0.0
        nr = len(r)
        if pdf is None:
            pdf = prob(fName)
        f = np.zeros(nr)
        for n in range(nr):
            f[n] = quad(pdf,r[0],r[n])[0]
        f /= f.max()
        print 'calculating inverse CDF'
        # Calculate the inverse CDF.
        # Assume CDF is approximately monotonic and sort to force it to be
        ind = f.argsort()
        f = f[ind]
        r = r[ind]
        # Drop values where CDF is constant (ie, prob = 0)
        mask = np.ones(nr,dtype='bool')
        for n in range(1,nr):
            if f[n] == f[n-1]:
                mask[n] = False
        f = f[mask]
        r = r[mask]
        finv = interp1d(f,r,kind='linear')
    return finv