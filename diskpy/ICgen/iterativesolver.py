#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Contains an iterative solver class to iteratively solve

Created on Thu Apr 20 10:36:48 2017

@author: ibackus
"""
from scipy.special import erf
from scipy.integrate import cumtrapz

import numpy as np
import pynbody
SimArray = pynbody.array.SimArray

G = SimArray(1.0,'G')
kB = SimArray(1.0, 'k')

class IterativeSolver():
    r"""
    Calculates the vertical density profile for a thin, extended gaseous disk
    orbiting a central star by solving vertical hydrostatic equilibrium 
    using an iterative method.
    
    The generator equation used is:
        
    .. math::
        \rho_{n+1}(z) &= \rho_{n+1}(0)\exp\left[ -a \left( \frac{1}{R} - 
        \frac{1}{\sqrt{R^2 + z^2}} \right) - 2 b II_{\rho}(z)\right] \\
        a             &\equiv \frac{G M m}{k_B T} \\
        b             &\equiv \frac{2\pi G m}{k_B T} \\
        II_{\rho}(z)  &\equiv \int_0^z\int_0^{z'} \rho(z'') dz'' dz'
        
    The normalization :math:`\rho_{n+1}(0)` is calculated by the constraint 
    that the integral of the density (numerically estimated) is equal to
    the surface density.
    
    This generator typically converges for reasonable disk/star parameters.
    I've tried other generators, and they all failed horribly.
    
    The initial guess is for a thin disk ignoring self gravity:
        
        .. math::
            \rho_0(z) &= \frac{\Sigma}{\sqrt{2 \pi} h} \exp(-z^2/2 h^2) \\
            h         &\equiv \sqrt{R^3/a} = c_s/\Omega
            
    The 2nd integral of this is known which provides a basis for the first
    order solution:
        
        .. math::
            II_{\rho,0}(z) &= z I_{\rho}(z) + h^2 \rho{0}(z) 
                - \frac{\Sigma h}{\sqrt{2 \pi}} \\
            I_{\rho,0}(z) &= \frac{\Sigma}{2}\mathrm{erf}(z/\sqrt{2}h)
            
    Which is then normalized numerically
    """
    def __init__(self, IC, R, z=None):
        
        sigma = IC.sigma(R)
        T = IC.T(R)
        M = IC.settings.physical.M
        m = IC.settings.physical.m
        # Initial calculations
        a = G*M*m/(kB*T)
        b = 2*np.pi*G*m/(kB*T)
        h = np.sqrt(R*R*R/a).in_units('au')
        self.h = h
        self._rhounits = (sigma/h).units
        # Setup z
        if z is None:
            zmax = IC.settings.rho_calc.zmax
            nz = IC.settings.rho_calc.nz
            if zmax is None:
                # Default value of zmax
                zmax = 6. * h
            z = np.linspace(0, zmax, nz)
            if not pynbody.units.has_units(z):
                raise ValueError, 'z has no units.  maybe zmax doesnt either'
        z = self._assume_h_units(z)
        rho0 = self.rho0(z, h, sigma)
        rzfactor = a * (1./R - 1./np.sqrt(R**2 + z**2))
        rzfactor.convert_units('1')
        # Save stuff
        self.a = a
        self.b = b
        self.R = R.in_units(h.units)
        self.T = T
        self.M = M
        self.m = m
        self.sigma = sigma
        self.rho = [rho0]
        self.z = z
        self._rzfactor = rzfactor
        
    def _assume_h_units(self, z):
        """
        For all the operations here, I will assume the units of z are the 
        units of h
        """
        if pynbody.units.has_units(z):
            return z.in_units(self.h.units) 
        else:
            z = SimArray(z, self.h.units)
            return z
        
    def iterate(self):
        """
        Perform one iteration.  Guesses are appended as SimArrays to the list
        self.rho
        """
        currentOrder = len(self.rho) - 1
        
        if currentOrder == 0:
            
            rho = self._rho1()
            
        else:
            
            rho = self._generator(self.rho[-1])
            
        self.rho.append(rho)
        
    def fit(self, maxiter=50, ftol=1e-8):
        """
        A convenience utility which iterates over self.iterate() and saves
        results to a results dict
        
        This will also delete any guess except the first one
        """
        self.rho = [self.rho[0]]
        i = 0
        converged = False
        
        while (i < maxiter) and (not converged):
            
            self.iterate()
            change = abs(similarity(self.rho[-2], self.rho[-1], True))
            converged = change < ftol
            i += 1
            
        if not converged:
            
            print 'WARNING: maxiter reached.  did not converge'
            
        self.results = {'rho': self.rho[-1], 'z': self.z}
    
    def _generator(self, rho):
        """
        Generates a next-order guess for rho
        """
        # Load        
        z = self.z
        sigma = self.sigma
        part1 = self._rzfactor
        # Double integral of rho
        IntIntRho = cumtrapz(cumtrapz(rho, z, initial=0.), z, initial=0.)
        IntIntRho = SimArray(IntIntRho, rho.units * (z.units)**2)
        # Calculate rho
        part2 = (2 * self.b * IntIntRho).in_units('1')
        rho = np.exp(-part1 - part2)
        # Normalize
        rho *= 0.5 * sigma/np.trapz(rho, z)
        return rho
        
    def rho0(self, z, h=None, sigma=None):
        """
        Initial guess (a gaussian).  this models a thin disk with no 
        self-gravity.  Generally, this will overestimate the disk scale height.
        """
        
        if h is None:
            h = self.h
        if sigma is None:
            sigma = self.sigma
        
        rho = sigma/(h * np.sqrt(2*np.pi))
        rho = rho* np.exp(-z**2/(2*h**2))
        
        return rho.in_units(self._rhounits)
    
    def _rho1(self):
        r"""
        Calculates the rho1 guess which is known semianlytically by running
        the generator on the initial guess
        """
        b = self.b
        z = self.z
        h = self.h
        sigma = self.sigma
        part1 = self._rzfactor
        
        # Calculate 2nd integral of rho0
        Irho = 0.5 * sigma * erf(z/(h * np.sqrt(2)))
        IntIrho = Irho*z + (h**2)*self.rho[0] - sigma*h/np.sqrt(2*np.pi)
        
        part2 = 2 *  b * IntIrho
        part2.convert_units('1')
        rho = np.exp(-part1 - part2)
        rho1norm = 0.5 * sigma/np.trapz(rho, z)
        
        return (rho * rho1norm).in_units(self._rhounits)
    
def similarity(p1, p2, forcePositive=False):
    """
    Attempts to quantify the similarity between 2 PDFs using a normalized and
    shifted version of the bhattacharyya distance, defined as:
        
    1 - D(p1, p2)/D(p1, p1)
    
    Where D(x, y) is the Bhattacharyya distance between x and y.
    
    If force positive is set, negative values are set to zero.
    
    A value << 1 indicates similar pdfs
    """
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    if forcePositive:
        p1[p1<0] = 0
        p2[p2<0] = 0
    elif np.any(p1 < 0) or np.any(p2 < 0):
        raise ValueError, "Negative values present.  try using force Positive"
        
    bhattacharyya = -np.log(np.sqrt(p1*p2).sum())
    norm = -np.log(p1.sum())
    return 1 - bhattacharyya/norm


def estHeight2(rho, z):
    """
    A simple way to estimate disk height from a density profile
    """
    h2 = np.trapz(rho*z**2, z) / np.trapz(rho, z)
    return np.sqrt(h2).in_units('au')

def heightProfile(IC):
    """
    Generates a disk height profile from the vertical density profile
    
    Returns
    -------
    h, r : SimArray
        Disk height and radius
    """
    r = IC.rho.r_bins
    nr = len(r)
    h = SimArray(np.zeros(nr), 'au')
    for i in range(nr):
        z = IC.rho.z_bins[:, i]
        rho = IC.rho.rho_binned[:,i]
        h[i] = np.sqrt( (rho*(z**2)).sum()/rho.sum() )
    return h, r