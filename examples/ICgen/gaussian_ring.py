# -*- coding: utf-8 -*-
"""
An example script for generating initial conditions.  This assumes ICgen has
already been properly setup

This generates initial conditions with a surface density profile of the form:

Sigma ~ exp(-(r-rd)^2/2a^2)

This generates a ring centered at 10 AU with a width of 1 AU.  The disk and star
masses are 0.1, 1 Msol respectively.  50k particles are used.  A powerlaw
temperature will be used

Created on Wed Sep 17 11:23:32 2014

@author: ibackus
"""
from diskpy import ICgen
import pynbody
SimArray = pynbody.array.SimArray

if __name__ == "__main__":
    
    # Initialize a blank initial conditions (IC) object:
    IC = ICgen.IC()
    
    # Echo the default settings to get an idea of what parameters can be
    # changed
    IC.settings()
    
    # Set up the surface density to be a gaussian ring
    IC.settings.sigma.kind = 'gaussring'
    
    # Now echo the sigma settings (defaults)
    IC.settings.sigma()
    
    # The important parameters to set are:
    #   Rd : Center of the disk ring, should be a pynbody SimArray with units
    #   ringwidth : width of the ring (standard deviation) SimArray with units
    #   m_disk : disk mass (SimArray with units)
    #   n_points : number of radial points to calculate sigma at
    IC.settings.sigma.n_points = 1000
    IC.settings.sigma.Rd = SimArray(10.0, 'au')
    IC.settings.sigma.ringwidth = SimArray(1.0, 'au')
    IC.settings.sigma.m_disk = SimArray(0.1, 'Msol')
    
    # Lets be careful and save what we've done.  This will save the ICs to
    # IC.p in the current directory
    IC.save()
    
    # Lets look at the settings for calculating density/vertical hydroequilibrum
    # The defaults should be fine, but if higher resolution is wanted, try
    # dropping r_bin_tol by a factor of 10
    # Echo defaults:
    IC.settings.rho_calc()
    
    # Set up the position generation for 50k particles
    IC.settings.pos_gen.nParticles = int(5e4)
    IC.settings.pos_gen.method = 'grid' # optional 'random'
    
    # Set up the temperature profile to use.
    # We'll use something of the form T = T0(r/r0)^Tpower
    IC.settings.physical.kind = 'powerlaw'
    IC.settings.physical.Tpower = -0.5  # exponent
    IC.settings.physical.T0 = SimArray(150.0, 'K')  # temperature at r0
    IC.settings.physical.Tmin = SimArray(10.0, 'K') # Minimum temperature
    IC.settings.physical.r0 = SimArray(1.0, 'au')
    
    # Let's set the star mass and gas mass assuming H2
    IC.settings.physical.M = SimArray(1.0, 'Msol') # star mass in solar masses
    IC.settings.physical.m = SimArray(2.0, 'm_p') # mass of H2
    
    # Lets pick a preset for changa to run on
    IC.settings.changa_run.preset = 'local'
    IC.settings.changa_run.preset = 'changa_uw-copy'
    
    # Save our work to IC.p
    IC.save()
    IC.settings()
    
    # We should be done, all we have to do now is tell the ICs to generate.
    # There are 2 ways to do this, and it may be fairly slow.
    
    # 1) One way is simply to call:
    IC.generate()
    # The ICs should already be saved, but just to be extra safe:
    IC.save()
    # This will run through the whole procedure and save a tipsy snapshot
    # to snapshot.std with a basic .param file saved to snapshot.param
    
    ## 2) Otherwise we can go step by step
#    IC.maker.sigma_gen() # Generate surface density profile and CDF
#    IC.maker.rho_gen() # Calculate density according to hydrodynamic equilibrium
#    IC.maker.pos_gen() # Generate particle positions
#    IC.maker.snapshot_gen() # Generate the final tipsy snapshot with velocities etc
#    IC.save()
    
    # This will run through the whole procedure and save a tipsy snapshot
    # to snapshot.std with a basic .param file saved to snapshot.param
