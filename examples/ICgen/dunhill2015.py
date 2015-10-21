"""
@author: dflemin3
Script for generating disk ICs about a stellar system
using ibackus's ICgen routines.
See https://github.com/ibackus/ICgen
"""
from diskpy.ICgen import ICgen, binary
import pynbody
SimArray = pynbody.array.SimArray
if __name__ == '__main__':
 
    # Initialize a blank initial conditions (IC) object:
    IC = ICgen.IC()
    
    IC.settings.physical.M = SimArray(3.6, 'Msol') #Total stellar mass in solar masses
    IC.settings.physical.m = SimArray(2.35, 'm_p') #mean molecular mass

    #Define masses of primary, secondary as pynbody SimArrays
    #Note, m1 + m2 == IC.settings.physical.M 
    #Only need to set if you're considernig a circumbinary system
    m1 = SimArray(2.2,'Msol')
    m2 = IC.settings.physical.M - m1
    
    #Scale the mass of the disk to be some fraction of the star mass
    IC.settings.snapshot.mScale = 0.005    

    #Define whether the star is a single star or binary
    IC.settings.physical.starMode = 'binary'
    
    #Set binary system parameters.  If single star, comment this out
    #Define list of orbital elements of the following form:
    #X = [e, a [AU], i, Omega, w, nu] where all angles are in degrees
    X = [0.6, 0.22, 0.0, 0.0, 0.0, 0.0]
    IC.settings.physical.binsys = binary.Binary(X,m1,m2,'kepler')
    
    # Lets generate a disk with powerlaw from [Rin,Rd] au followed by a cutoff
    # Set up the surface density profile settings.  Notice that right now the
    # Lets use a simple powerlaw with cut-offs at the interior and edge of the
    # disk
    # The important parameters to set are:
    #   Rd : Disk radius, should be a pynbody SimArray
    #   Qmin : Minimum estimated Toomre Q in the disk
    #   n_points : number of radial points to calculate sigma at
    #	power: sigma ~ r^(power)
    IC.settings.sigma.kind = 'powerlaw'
    IC.settings.sigma.power = -1.7
    IC.settings.sigma.Qmin = 1.5
    IC.settings.sigma.n_points = 1000
    
    IC.settings.sigma.Rd = SimArray(2.0,'au') #Outer edge of powerlaw part of disk
    IC.settings.sigma.rmax = 2.0 #Set rmax 
    IC.settings.sigma.rin = 0.25 #Set inner disk radius
    IC.settings.cutlength = 0.01 #Set exp cutoff length scale
    IC.settings.pos_gen.method = 'random' #Instead of grid sampling, use random
    
    #This will save the ICs to
    # IC.p in the current directory
    IC.save()
    
    # Change the settings used for numerically calculating the gas density    
    # Set the number of gas particles
    IC.settings.pos_gen.nParticles = 100000
    
    # Set up the temperature profile to use.  Available kinds are 'powerlaw'
    # and 'MQWS'
    # We'll use something of the form T = T0(r/r0)^Tpower
    IC.settings.physical.kind = 'powerlaw'
    IC.settings.physical.Tpower = -1  # exponent
    IC.settings.physical.T0 = SimArray(750, 'K')  # temperature at r0
    IC.settings.physical.Tmin = SimArray(150.0, 'K') # Minimum temperature
    IC.settings.physical.r0 = SimArray(1.0, 'au')
    
    # Lets have changa run on the local preset
    IC.settings.changa_run.preset = 'local'
    
    # Save our work to IC.p
    IC.save()
    IC.settings()
    
    #Actually generate snapshot
    IC.generate()
    IC.save()
    # This will run through the whole procedure and save a tipsy snapshot
    # to snapshot.std with a basic .param file saved to snapshot.param
