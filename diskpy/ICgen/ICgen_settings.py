# -*- coding: utf-8 -*-
"""
Defines the settings class for ICgen

USAGE:

Load settings from file:
    
    import ICgen_settings
    settings = ICgen_settings.settings('filename')
    
Load settings from defaults (Defined in ICgen_settings.py):

    import ICgen_settings
    settings = ICgen_settings.settings()

AFTER LOADING:    
Access settings:

    # Return number of particles:
    settings.pos_gen.nParticles
    
    # Etc...
    
Echo settings:

    # Echo ALL settings
    settings()
    
    # Echo filenames
    settings.filenames()
    
    # Etc...
    
READ/WRITE

Write (save) settings to disk:

    settings.write('filename')
    
Load (read) settings from disk:

    settings.load('filename')
"""

__version__ = "$Revision: 5 $"
# $Source$
__iversion__ = int(filter(str.isdigit,__version__))

from ICgen_utils import checkversion, pickle_import

import pynbody
SimArray = pynbody.array.SimArray
import cPickle as pickle
import numpy as np
import os

_dir = os.path.dirname(os.path.realpath(__file__))

""" **************************************************
SETTINGS

ALL of these variable names MUST be set!  
************************************************** """

class filenames:
    """
    Filenames
    
    To echo settings:
        filenames()
    To access settings:
        filenames.settingname
    """
    
    def __init__(self):
        

        # Initial conditions filename
        self.IC_file_name = 'IC.p'
        # Filename to save tipsy snapshot to
        self.snapshotName = 'snapshot.std'
        # Filename to save ChaNGa .param file to.  If None, no file saved
        # To edit settings, you can change default.param or the output .param file
        self.paramName = 'snapshot.param'
        # Default .director filename
        self.directorName = 'snapshot.director'
        
    def __call__(self):
        
        print_settings(self, 'Filenames:')
            
class physical:
    """
    Defines default physical parameters
    
    To echo settings:
        physical()
    To access settings:
        physical.settingname
    """
    
    def __init__(self, kind=None, binsys=None, starMode = 'single'):
        #editted by dflemin3 07/10/2015
        # Molecular mass of the gass.  If m = None, Assumed to be H2, m = 2.00132 m_p
        self.m = SimArray(2.00132,'m_p')
        # Mass of the star (or total mass of binary system).  If M = None, Assumed to be 0.33 Msol
        self.M = SimArray(0.33 , 'Msol')
        # CONSTANTS FOR CALCULATING TEMPERATURE.  T(r) = T0(r/r0)^Tpower.
        # See calc_temp.py.  If None, defaults to settings in calc_temp.py
        self.T0 = SimArray(332.406,'K')
        self.r0 = SimArray(0.5,'au')
        self.Tpower = -0.59
        self.Tmin = SimArray(0, 'K')    # Minimum temperature cut-off
        self.Tmax = SimArray(np.inf, 'K')
        # Equation of state parameters
        self.eos = 'isothermal' # isothermal or adiabatic
        self.gamma = 1.4
        
        # Binary parameters
        self.starMode = starMode #single star or binary?
        #Binary Orbital Parameters...store in Binary class from binary.py
        self.binsys = binsys
        
        if kind is None:
            
            kind = 'powerlaw'
            
        self.kind = kind # Type of temperature profile.  See calc_temp.py
        
    def __call__(self):
        
        print_settings(self, 'General physical parameters:')
        
class sigma:
    """
    Settings for generating a surface density profile
    
    To echo:
        sigma()
    To access settings:
        sigma.settingname
    """
    def __init__(self, kind='powerlaw'):
        
        # Kind of surface density profile to use.  Options:
        #   'powerlaw'
        #   'MQWS'
        #   'viscous'
        # Setting self.kind automatically initializes the default settings
        # for that kind
        self.kind = kind
        
    def _set_defaults(self):
        # Set default attributes:
        # First delete all attributes but kind
        kind = self.kind
        for key in self.__dict__.keys():
            
            if key != 'kind':
                
                self.__dict__.pop(key, None)
                
        # Now set attributes
        if kind == 'powerlaw':
            
            self.Rd = SimArray(1.0,'au')
            self.rin = 0.5
            self.rmax = 2.3
            self.cutlength = 0.3
            self.Qmin = 1.5
            self.n_points = 1000
            self.power = -1 #Power on powerlaw of form sigma ~ r^(power)
       
     
        if (kind == 'mqws') | (kind == 'MQWS'):
            
            self.rin = 4.0
            self.rout = 20.0
            self.rmax = None
            self.power = -1
            self.m_disk = SimArray(0.1, 'Msol')
            self.n_points = 1000
            self.Qmin = 1.5
            
        if kind == 'viscous':
            
            self.Rd = SimArray(1.0, 'au')
            self.rin = 0.1
            self.rmax = 2.0
            self.power = -1
            self.m_disk = SimArray(0.1, 'Msol')
            self.n_points = 1500
            self.gamma = 0.9
            
        if kind == 'gaussring':
            
            self.Rd = SimArray(1.0, 'au')
            self.ringwidth = SimArray(0.1, 'au')
            self.m_disk = SimArray(0.1, 'Msol')
            self.n_points = 1000
            
        # innercut and outercut determine where to apply a hard cut to the
        # surface density (if anywhere).  for example, to set sigma=0 for
        # R > 2 au, do: outercut = SimArray(2,'au')
        self.innercut = None
        self.outercut = None
            
    def __setattr__(self, attr, value):
        """
        Override the default __setattr__ so that changing self.kind causes
        the defaults to be set for that kind
        """
        self.__dict__[attr] = value
        
        if attr == 'kind':
            
            self._set_defaults()
            
        
    def __call__(self):
        
        print_settings(self, 'Sigma profile parameters:')

class rho_calc:
    """
    Settings for calculating rho(z,r)
    
    To echo settings:
        rho_calc()
    To access settings:
        rho_calc.settingname
    """
    
    def __init__(self):
        

        # The number of radial data points to calculate rho(z,r)
        # If None (default), reasonable r bins are automatically estimated
        self.nr = None
        # The number of vertical points to calculate rho(z,r) at
        self.nz = 1000
        # The maximum z (assumed to be au) to calculate rho(z,r) at.
        # If None (default, safe) zmax is automatically estimated at each rbin
        self.zmax = None
        # Minimum number of rbins to use
        self.min_r_bins = 150
        # Accuracy requirement for r bins
        self.r_bin_tol = 1e-3
        
    def __call__(self):
        
        print_settings(self, 'Rho calculation settings:')
        
class pos_gen:
    """
    Settings for generating random positions [STEP 2]
    
    To echo settings:
        pos_gen()
    To access settings:
        pos_gen.settingname
    """
    
    def __init__(self):
        
        # Number of random positions to generate:
        self.nParticles = 40411
        # The method for generating positions
        self.method = 'grid'
        
    def __call__(self):
        
        print_settings(self, 'Position generator settings:')
        
class snapshot:
    """
    Settings for generating tipsy files (includes particle mass, temp, and vel.)
    [STEP 3]
    """
    
    def __init__(self, nParticles=0):
        
        # Factor by which to scale the disc mass before time evolving it up to
        # a final mass of Mdisc (see above).  Should be between 0 and 1
        # Default: 1.0 (no scaling done).  Final particle masses are scaled
        # by this factor
        self.mScale = 1.0
        self.nSmooth = 32
        
        # Other defaults that shouldn't need to be changed
        self.metals = 1.0
        
    def __call__(self):
        
        print_settings(self, 'Tipsy snapshot generator settings:')
                
class changa_run:
    """
    Settings for running ChaNGa (needed when calculating velocities, eps, so on)
    """
    
    def __init__(self, preset = 'local'):
        
        # Run configuration to use.  See ICgen_utils.changa_command for options
        self.preset = preset
        # Additional arguments for all ChaNGa calls
        self.changa_args = ''
        # Additional arguments for all runner (mpirun, charmrun, ...) calls
        self.runner_args = ''
        # Display ChaNGa output
        self.verbose = True
        # Save to log file
        self.logfile_name = None
        
    def __call__(self):
        
        print_settings(self, 'ChaNGa run options:')
        
class settings:
    """
    settings for ICgen
    
    USAGE:
    
    Load settings from file:
        
        import ICgen_settings
        settings = ICgen_settings.settings('filename')
        
    Load settings from defaults (Defined in ICgen_settings.py):
    
        import ICgen_settings
        settings = ICgen_settings.settings()
        
    Load settings for a given surface density profile
    
        import ICgen_settings
        settings = ICgen_settings.settings('powerlaw')
        or try:
        settings = ICgen_settings.settings('MQWS')
    
    AFTER LOADING:    
    Access settings:
    
        # Return sigmaFileName:
        settings.filenames.sigmaFileName
        
        # Return number of particles:
        settings.pos_gen.nParticles
        
        # Etc...
        
    Echo settings:
    
        # Echo ALL settings
        settings()
        
        # Echo filenames
        settings.filenames()
        
        # Etc...
        
    READ/WRITE
    
    Write (save) settings to disk:
    
        settings.write('filename')
        
    Load (read) settings from disk:
    
        settings.load('filename')
    """    
    def __init__(self, settings_filename=None, kind=None):
        
        self.__version__ = __iversion__
        
        if settings_filename is None:
            # Load the defaults
            self.filenames = filenames()
            self.sigma = sigma(kind)           
            self.rho_calc = rho_calc()
            self.pos_gen = pos_gen()
            self.snapshot = snapshot(self.pos_gen.nParticles)
            self.physical = physical(kind)
            self.changa_run = changa_run()
            
        else:
            
            self.load(settings_filename)
        
    def __call__(self):
        
        self.filenames()
        print ''
        # Check if sigma is a setting (needed for compatibility)
        if hasattr(self,'sigma'):
            
            self.sigma()
            print ''
            
        self.rho_calc()
        print ''
        self.pos_gen()
        print ''
        self.snapshot()
        print ''
        self.physical()
        print ''
        self.changa_run()
        
    def save(self,settings_filename = None):
        """
        Save settings to settings_filename.  If settings_filename is None
        (or is not set), settings_filename will be taken as the filename of
        the most recent save/load
        """
        
        if settings_filename is None:
            # Try filename stored in settings object (self)            
            if not hasattr(self,'settings_filename'):
                
                raise RuntimeError,'No settings_filename defined.  Cannot save'
                
            else:
                
                settings_filename = self.settings_filename
                
        f = open(settings_filename,'wb')
        pickle.dump(self.__dict__,f,2)
        f.close()
        print 'Settings saved to {0}'.format(settings_filename)
        # Update the settings_filename
        self.settings_filename = settings_filename
        
    def load(self, settings_filename):
        """
        Load settings from settings_filename.
        """
        if isinstance(settings_filename, str):
            # Load the settings
            tmp_dict = pickle_import(settings_filename, _dir)
            
        else:
            # Assume the an already loaded settings dictionary is being passed
            tmp_dict = settings_filename
        
        # Sigma has to be handled specially because the settings are dynamic
        if 'sigma' in tmp_dict:
            
            self.sigma.kind = tmp_dict['sigma'].kind
        
        # Get the names of the settings objects
        keys = []
        
        for key in self.__dict__.keys():
            
            if key[0] != '_':
                
                keys.append(key)
                
        # Now update the dictionary with the settings object dictionaries
        for key in keys:
            
            if key in tmp_dict:
                
                old = getattr(self, key)
                new = tmp_dict[key]
                # Update the old (possibly default) settings
                old.__dict__.update(new.__dict__)
                setattr(self, key, old)
                
        
def print_settings(setting, header=None):
    
    print '------------------------'
    print header
    print ''
    
    for key,val in setting.__dict__.iteritems():
            
            if pynbody.units.has_units(val):
                print '{0} : {1} {2}'.format(key,val,val.units)
            else:
                print '{0} : {1}'.format(key,val)
