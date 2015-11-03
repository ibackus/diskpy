# -*- coding: utf-8 -*-
"""
Implements functions for estimate spiral density wave power in PPDs

Created on Mon Nov  2 13:11:22 2015

@author: ibackus
"""

import pynbody
SimArray = pynbody.array.SimArray
import numpy as np

from diskpy.pdmath import bin2dsum, dA

def powerspectrum(f, mMax=30, rbins=50):
    """
    The density power spectrum along the angular direction, summed along the
    radial direction.
    
    Parameters
    ----------
    f : SimSnap
        Snapshot of a disk
    mMax : int
        Maximum fourier mode to calculate
    rbins : int or array
        Number of radial bins or the binedges to use
        
    Returns
    -------
    
    m : array
        Array of the fourier modes from 0 to mMax (integers)
    power : array
        Power in the fourier modes, summed along radial direction.  Power is
        take to be the square of the surface density fourier transform
    """
        
    r, m, sigtransform = sigmafft(f, rbins, 2*mMax + 1)
    m = m[0,:]
    power = (abs(sigtransform)**2).sum(0)
    
    return m, power

def sigmafft(f, rbins=50, thetabins=50):
    """
    Calculates the fourier transform of the surface density along the angular
    direction.  Works by binning sigma in r, theta (using sigmacylindrical) and
    doing a fourier transform along the theta axis.
    
    Parameters
    ----------
    f : SimSnap
        Simulation snapshot of a disk
    rbins, thetabins : int or array like
        Number of bins or binedges to use
    
    Returns
    -------
    rmesh : SimArray
        2D Meshgrid of radial binedges
    mmesh : Array
        2D meshgrid of m binedges, where m is an integer (the mth fourier mode)
    sigfft : SimArray
        2D meshgrid of surface density fourier transformed along the angular
        direction.
    
    Notes
    -----
    
    The returned arrays are indexed such that array[i,j] gives the value of
    array at radial bin i and theta bin j.
    """
    
    rmesh, thetamesh, sigma = sigmacylindrical(f, rbins, thetabins)
    sigfft = SimArray(np.fft.rfft(sigma), sigma.units)
    nm = sigfft.shape[1]
    m = np.arange(0, nm)
    rmesh2, mmesh = np.meshgrid(rmesh[:,0], m)
    
    return rmesh2.T, mmesh.T, sigfft

def sigmacylindrical(f, rbins=50, thetabins=50):
    """
    Estimates the surface density, binned in cylindrical coordinates
    
    Parameters
    ----------
    
    f : SimSnap (see pynbody)
        Snapshot of a disk
    rbins, thetabins : int or arraylike
        Number of bins or binedges
    
    Returns
    -------
    
    rmesh : SimArray
        2D Meshgrid of radial binedges
    thetamesh : Array
        2D meshgrid of angle binedges
    sigma : SimArray
        2D meshgrid of surface density.
    
    Notes
    -----
    
    The returned arrays are indexed such that array[i,j] gives the value of
    array at radial bin i and theta bin j.
    """
    r = f.g['rxy']
    theta = np.asarray(np.arctan2(f.g['y'], f.g['x']))
    theta = theta % (2*np.pi)
    # Default theta bin edges are 0 to 2 pi, evenly spaced
    if isinstance(thetabins, int):
        
        thetabins = np.linspace(0, 2*np.pi, thetabins + 1)
        
    # Now bin/sum particle masses
    msum, redges, thetaedges = bin2dsum(r, theta, f.g['mass'], rbins, thetabins)
    # sigma = mass/area
    sigma = msum / dA(redges, thetaedges)
    # Do mesh grid
    rmesh, thetamesh = np.meshgrid(redges, thetaedges)
    return rmesh.T, thetamesh.T, sigma
