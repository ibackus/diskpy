# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:39:23 2015

@author: ibackus
"""
import pynbody as pb
SimArray = pb.array.SimArray
import numpy as np


G = SimArray(1.0, 'G')
kB = SimArray(1.0, 'k')

def rho0_est(h, sigma):
    """
    Estimates the midplane density vs r as 
    
    .. math::
    
        \\rho(z=0) = \\frac{\\Sigma}{h \\sqrt{2 \\pi}}
        
        
    Parameters
    ----------
    h, sigma : array(s)/SimArrays (same shape)
        Disk height and surface density (estimated)
    """
    out = sigma/h
    mask = np.isnan(out)
    out[mask] = 0
    return out
    
def h_est(r, M, m, T, gamma=1):
    """
    Estimates disk height :math:`h=c_{s}/\\Omega`
    
    Parameters
    ----------
    
    r : SimArray, array, float
        Disk radii to calculate h at
    M : SimArray, float
        Central star mass
    m : SimArray, float
        Mean molecular mass (with units of mass)
    T : SimArray, array, float
        Temperature as a function of r
    gamma : float
        (optional) Adiabatic index to be used.  If temperature is not a
        functino of z, gamma should equal 1 (isothermal)
    
    Returns
    -------
    
    h : SimArray, float, array
        Disk height vs. r (1st order estimate)
    """
    a = G*M*m/(kB*T)
    hout = np.sqrt(r*r*r/a)
    return hout.in_units(r.units)