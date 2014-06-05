# -*- coding: utf-8 -*-
"""
A script for generating initial conditions for an axisymmetric
protoplanetary disc, assuming an isothermal equation of state and
steady-state.

Please note, as long as the previous steps have been completed, you
can start at any step

Look carefully at the settings saved in settings_file (by default,
 ICgen_settings.py)

The temperature profile is set in calc_temp.py

The surface density profile is either (1) provided as a pickled dictionary
containing sigma (surface density) and r, where sigma is evaluated
at r or (2) defined as a function in calc_sigma.py

IMPORTANT: to prevent any units errors, everything should be in units
of Msol and 'au', except the mean molecular weight and Temperatures.

*****     STEPS     *****
1) Numerically calculate the density, rho(z,r)
2) Randomly generate particle positions according to rho(z,r)
3) Re-scale particle mass, calculate velocities, and save tipsy snapshot
    -particles are lightened to allow us to slowly increase particle mass
    over time in ChaNGa

@author: ibackus
"""
# DEBUGGING
import sys
if globals().has_key('init_modules'):
	for m in [x for x in sys.modules.keys() if x not in init_modules]:
		del(sys.modules[m]) 
else:
	init_modules = sys.modules.keys()
# END DEBUGGING
# External modules
import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import scipy.interpolate as interp
import os
import cPickle as pickle
# ICgen modules
import calc_rho_zr
#import pos_gen
import pos_gen_grid as pos_gen
import calc_temp
import calc_sigma
ICgenDir = os.path.dirname(os.path.realpath(__file__))

# File containing all the settings definitions
settings_file = 'ICgen_settings.py'

execfile(settings_file)

""" **************************************************
EXECUTION
************************************************** """
if initial_step < 2:
    print '***********************************************'
    print 'STEP 1 OF 3: Generating PDF(r,z) [rho(z,r)]'
    print '***********************************************'
    rho = calc_rho_zr.rho_zr(sigmaFileName, nr=nr, nz=nz, zmax=zmax, \
    m=m, M=M, rho_tol=rho_tol, output=rhoFileName, rmin=rmin, rmax=rmax, \
    T0=T0, r0=r0, Tpower=Tpower)
else:
    rho = pickle.load(open(rhoFileName,'rb'))
if initial_step < 3:
    print '***********************************************'
    print 'STEP 2 OF 3: Generating random positions'
    print '***********************************************'
    pos = pos_gen.make(rhoFileName, sigmaFileName, nParticles, zlim=zlim, \
    rlim=rlim, savename=posFileName)
else:
    pos = pickle.load(open(posFileName,'rb'))
if initial_step < 4:
    print '***********************************************'
    print 'STEP 3 OF 3: Generating tipsy snapshot'
    print '***********************************************'
    # -------------------------------------------------
    # Initialize/calculate disc properties
    # -------------------------------------------------
    # Constants
    G = SimArray(1.0,'G')
    kB = SimArray(1.0,'k')
    # Assign variable names
    r = pos['r']
    z = pos['z']
    # Mass of the star
    mStar = M
    # Mass of the SPH gas particles
    if Mdisc is None:
        # Estimate mDisc from sigma
        sigma = calc_sigma.sigma(sigmaFileName)
        rbins = np.linspace(0,pos['r'].max(),nr)
        Mdisc = np.trapz(2*np.pi*rbins*sigma(rbins),rbins)
        Mdisc = SimArray(float(Mdisc),'Msol')
        print 'Mdisc = {0}'.format(Mdisc)
    mGas = np.ones(nParticles)*Mdisc/nParticles
    mGas *= mScale
    print 'Scaling particle masses by {0}'.format(mScale)
    # Find the total mass interior to every particle
    N_interior = np.array(np.argsort(pos['r']).argsort())
    mInt = mGas[[0]]*N_interior + mStar
    # -------------------------------------------------
    # Calculate velocities
    # -------------------------------------------------
    print 'Estimating radial derivative of rho'
    # Calculate the radial derivative of rho(z,r) from the PDF calculated
    # above (for finding pressure gradients)
    rho_spline = interp.RectBivariateSpline(rho['z'],rho['r'],rho['rho'])
    #tck = rho_spline.tck + rho_spline.degrees
    dr = rho['r'][[1]] - rho['r'][[0]]
    drho_dr_grid = np.gradient(rho['rho'])[1]/dr
    drho_dr_spline = interp.RectBivariateSpline(rho['z'],rho['r'],drho_dr_grid)
    drho = np.zeros(nParticles)
    rho_array = np.zeros(nParticles)
    for n in range(nParticles):
#        if np.mod(n,100) == 0:
#            print 'Estimating radial derivative of rho {0} of {1}'\
#            .format(n,nParticles)
        drho[n] = drho_dr_spline(z[n],r[n])
        rho_array[n] = rho_spline(z[n],r[n])
    rho_array = SimArray(rho_array,rho['rho'].units)
    drho = SimArray(drho,drho_dr_grid.units)
    print 'Calculating circular velocity'
    # Generate temperature for every particle
    T = calc_temp.T(pos['r'], T0=T0, r0=r0, Tpower=Tpower)
    # Calculate keplerian velocity for each particle
    v2grav = G*mInt/r
    v2dens = (kB*T/m)*(r*drho/rho_array)
    #       ignore nans and infs
    v2dens[(np.isnan(v2dens)) | (np.isinf(v2dens))] = 0.0
    v2temp = (kB*T/m)*Tpower
    v = np.sqrt(v2grav + v2dens + v2temp).in_units('2.98e+01 km s**-1')
    # Sometimes, at large r, the velocities due to the pressure and temp
    # Gradients become negative.  If this is the case, set them to 0
    nanind = np.isnan(v)
    v[nanind] = 0.0
    print 'Setting {0} particle velocities to 0 (to avoid nans)'.format(nanind.sum())
    # -------------------------------------------------
    # Assign output
    # -------------------------------------------------
    vel = SimArray(np.zeros([nParticles,3]),'2.98e+01 km s**-1')
    vel[:,0] = -np.sin(pos['theta'])*v
    vel[:,1] = np.cos(pos['theta'])*v
    # Generate positions
    xyz = SimArray(np.zeros([nParticles,3]),'au')
    xyz[:,0] = pos['x']
    xyz[:,1] = pos['y']
    xyz[:,2] = pos['z']
    # Generate snapshot
    snapshot = pynbody.new(star=1,gas=nParticles)
    snapshot.gas['vel'] = vel
    snapshot.gas['pos'] = xyz
    snapshot.gas['temp'] = T
    snapshot.gas['mass'] = mGas
    snapshot.gas['metals'] = metals
    snapshot.gas['eps'] = eps
    
    snapshot.star['pos'] = SimArray([[ 0.,  0.,  0.]],'au')
    snapshot.star['vel'] = SimArray([[ 0.,  0.,  0.]], '2.98e+01 km s**-1')
    snapshot.star['mass'] = mStar
    snapshot.star['metals'] = SimArray([ 1.])
    snapshot.star['eps'] = SimArray([ 0.01], 'au')
    # -------------------------------------------------
    # Save Output
    # -------------------------------------------------
    # Save snapshot
    snapshot.write(filename=snapshotName,fmt=pynbody.tipsy.TipsySnap)
    print 'Tipsy snapshot saved to {0}'.format(snapshotName)
    # Start saving log file
    logFile = open(logFileName,'w')
    logFile.write('Generated from settings file: {0}\n'.format(os.path.abspath(settings_file)))
    m_amu = float((m*0.992733/pynbody.units.m_p).in_units('1'))
    # Save director file
    if paramName is not None:
        try:
            import isaac
            defPar = os.path.join(ICgenDir,'default.param')
            pars = isaac.configparser(defPar,ftype='param')
            pars['achInFile'] = snapshotName
            pars['achOutName'] = os.path.splitext(snapshotName)[0]
            # Convert molecular mass from mass of proton to amu
            pars['dMeanMolWeight'] = m_amu
            isaac.configsave(pars,paramName,ftype='param')
            print 'ChaNGa parameter file saved to {0}'.format(paramName)
            logFile.write('ChaNGa parameter file saved to {0}\n'.format(paramName))
        except:
            print 'Could not save .param file.  Check that the module isaac\
            was properly imported'
    #Continue writing to log file
    logFile.write('Tipsy snapshot saved to {0}\n'.format(snapshotName))
    logFile.write('Surface density profile saved to {0}\n\
rho(z,r) saved to {1}\n\
random positions saved to {2}\n'.format(sigmaFileName,rhoFileName,posFileName))
    logFile.write('Mstar = {0} {1}\n\
Mean molecular mass: {2} amu\n'.format(mStar,mStar.units,m_amu))
    logFile.close()
    # Generate example changa shell script
print '***********************************************'
print 'FINISHED!'
print '***********************************************'
