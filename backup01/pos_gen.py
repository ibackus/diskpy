# -*- coding: utf-8 -*-
"""
Defines a function to randomly generate particle positions according to 
the desired surface density profile (sigma vs r) and the vertical profile
(rho vs r,z).

Created on Mon Jan 27 18:48:04 2014

@author: ibackus
"""

import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import cPickle as pickle
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
import scipy.optimize as opt

# LOCAL IMPORTING
import calc_sigma
import calc_rho_zr

def make(rhoFileName,sigFileName,nParticles,rlim = None,zlim = None,\
savename = None):
    """
    Randomly generate particle positions according to the surface density
    profile contained in sigFileName and the vertical profile contained in
    rhoFileName.  PDFs are interpolated using a cubic spline, and are 
    assumed to be 0 outside the regions defined the the input files
    
    * Arguments *
    rhoFileName - file name for a file containing a pickled dictionary
        for rho calculated at points on a grid.  The dictionary should
        contain:
            dict['rho'] : 2D array, rho at all pairs of points (z,r)
            dict['z']   : a 1D array of z points
            dict['r']   : a 1D array of r points
            
    sigFileName - file name for a file containing a pickled dictionary
        for sigma (surface density) calculated at points r.  Dictionary
        should contain:
            dict['sigma']   : 1D array, rho at all points r
            dict['r']       : 1D array of r points
    
    nParticles - total number of particles to generate
    
    rlim - cylindrical radius bounds for generating particles.  if None,
        then the default rlims are the range of values contained in 
        rhoFileName.
        
    zlim - positive z-bounds for generating particles.  Particles will be 
        always generated symmetrically, within these bounds (the sign of
        z is randomly chosen).  If None, bounds are taken from rhoFileName
    
    savename - if None, no file is saved.  Otherwise, the output is pickled
        and saved to 'savename'
    
    * Output *
    Returns a dictionary containing the randomly generated positions, with
    the following keys:
        
        'x'     : x-positions of the particles
        'y'     : y-positions
        'z'     : z-positions
        'r'     : monotonically increasing R-positions
        'theta' : angular positions
    
    
    !!!Everything should be in units of Msol and au!!!
    pynbody simarrays are acceptable
    """
    # ------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------
    a = pickle.load(open(rhoFileName,'rb'))
    rho = a['rho']
    z = a['z']
    r = a['r']
    if rlim == None:
        rlim = [r.min(),r.max()]
    if zlim == None:
        zlim = [z.min(),z.max()]
    # ------------------------------------------------------------
    # Generate PDFs
    # ------------------------------------------------------------
    #print('Generating rho(z,r) spline interpolation')
    #rho_zr = RectBivariateSpline(z,r,rho)
    cdfinv_zr = calc_rho_zr.cdfinv_zr(rhoFileName)
    
    # Probability as a function of r (2*pi*r*sigma, up to a constant factor)
    # Pr is normalized, but we need the MAX value of it
    print 'Calculating probability vs r spline'    
    Pr = calc_sigma.prob(sigFileName)
    def PrNeg(x):
        return -1.0*Pr(x)
    print 'Finding maximum probability for r = {0} to {1}'.format(rlim[0],rlim[1])
    #print PrNeg(rlim[1])
    xmax = opt.fminbound(PrNeg,float(rlim[0]),float(rlim[1]))
    Prmax = Pr(xmax)
    print 'Max. prob = {0} \nat r = {1}'.format(Prmax,xmax)
    # ------------------------------------------------------------
    # Generate random R values according to PDF from sigma vs r
    # ------------------------------------------------------------
    print('Calculating R') 
    R = np.zeros(nParticles)
    cdfinv_r = calc_sigma.cdfinv_r(sigFileName,Pr)
    R = cdfinv_r(np.random.uniform(size=nParticles))
    # ------------------------------------------------------------
    # Generate random Z values for the R values found above,
    # according to rho(z,r)
    # ------------------------------------------------------------
    print('Calculating Z')
    Z = np.zeros(nParticles)
    cdfvals = np.random.uniform(size=nParticles)
    rbins = np.digitize(R,r)
    dr = r[1]-r[0]
    for i in range(nParticles):
        rbin = rbins[i]
        zlo = cdfinv_zr[rbin-1](cdfvals[i])
        zhi = cdfinv_zr[rbin](cdfvals[i])
        Z[i] = zlo + ((zhi-zlo)/dr)*(R[i] - r[rbin-1])
    # ------------------------------------------------------------
    # Format the output
    # ------------------------------------------------------------
    # Assign units
    if pynbody.units.has_units(r):
        rUnits = r.units
    else:
        rUnits = pynbody.units.au
    if pynbody.units.has_units(z):
        zUnits = z.units
    else:
        zUnits = pynbody.units.au
    if pynbody.units.has_units(Z):
        Z.convert_units(zUnits)
    else:
        Z = SimArray(Z,zUnits)
    if pynbody.units.has_units(R):
        R.convert_units(rUnits)
    else:
        R = SimArray(R,rUnits)
    Z = Z*np.sign(np.random.randn(nParticles))
    theta = 2.0*np.pi*np.random.rand(nParticles)
    x = R*np.cos(theta)
    y = R*np.sin(theta)
    
    # Output
    outDict = {'x': x,'y': y, 'z': Z, 'r': R, 'theta': theta}
    if savename is not None:
        pickle.dump(outDict,open(savename,'wb'))
    return outDict
            
    
    