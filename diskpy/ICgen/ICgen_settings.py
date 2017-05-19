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
    
Echo setting defaults (with descriptions)
    settings.print_defaults()
    settings.sigma.print_defaults()
    
READ/WRITE

Write (save) settings to disk:

    settings.write('filename')
    
Load (read) settings from disk:

    settings.load('filename')
    
    
------------------

A note on implementation.  There are 2 base settings classes: settingsBase and
kindSettingsBase.  The first is just for generic settings, the second is
to allow users to have multiple layouts (see class sigma)

See the base classes for more info
"""

__version__ = "$Revision: 7 $"
# $Source$
__iversion__ = int(filter(str.isdigit,__version__))

from ICgen_utils import pickle_import

import pynbody
SimArray = pynbody.array.SimArray
import cPickle as pickle
import numpy as np
import os

_dir = os.path.dirname(os.path.realpath(__file__))

def default_settings():
    """
    Echos the current default settings with descriptions for ICgen settings
    """
    x = settings()
    x.print_defaults()
    

class settingsBase:
    """
    
    """
    def __init__(self, descr=None, **kwargs):
        
        self._description = descr
    
    def __call__(self):
        
        print_settings(self)
        
    def __repr__(self):
        
        return repr_settings(self)
    
    def _set_defaults(self, defaults=None, **kwargs):
        """
        Sets defaults from a defaults list.  Defaults lists are a list of 
        length 3 tuples:
        [ ('variable_name', default_val, 'Doc string for variable'),
          ('variable_name2', default_val2, 'Other doc string') ]
        
        kwargs can change available variables but not create new ones
        
        if defaults is None, defaults are taken from self._defaults
        """
        if defaults is None:
            if hasattr(self, '_defaults'):
                defaults = self._defaults
        if defaults is None:
            # No defaults found
            return
            
        # Extract defaults
        default_dict = {default[0]: default[1] for default in defaults}
        self.__dict__.update(default_dict)
        # Now read in kwargs
        goodkeys = self.__dict__.keys()
        for k in kwargs.keys():
            
            if k not in goodkeys:
                
                raise ValueError, 'Unrecognized option {}. Allowed options:{}'\
                    .format(k, goodkeys)
                
            # Ignore Nones
            if kwargs[k] is None:
                
                kwargs.pop(k, None)
                
        self.__dict__.update(kwargs)
        
    def print_defaults(self):
        """
        Prints defaults defined in self._defaults.   Defaults lists are a list of 
        length 3 tuples:
        [ ('variable_name', default_val, 'Doc string for variable'),
          ('variable_name2', default_val2, 'Other doc string') ]
        
        Also prints childrens defaults if they are settingsBase instances
        """
        # Print your own defaults
        if hasattr(self, '_defaults') and self._defaults is not None:
            
            string = repr_defaults(self._defaults)
            print string
            
        # if you contain other settings objects, print their defaults
        for k, v in self.__dict__.iteritems():
            
            if isinstance(v, settingsBase):
                
                print '\n----------------------------'
                print k
                print ''
                v.print_defaults()
        
class kindSettingsBase(settingsBase):
    """
    
    """
    _defaults = []
    _globaldefaults = []
    _defaultlayouts = {}
    def __init__(self, descr=None):
        """
        
        """
        settingsBase.__init__(self, descr)
        self._header = descr
        self._description = '{0}\nAvailable kinds:\n{1}'\
            .format(self._header, self._defaultlayouts.keys())
        
    def _set_defaults(self):
        """
        Deletes all public attributes (except kind) and updates self.__dict__
        with self._defaults[kind]
        """        
        kind = self.kind
        
        for key in self.__dict__.keys():
            
            if (key != 'kind') and (key[0] != '_') :
                
                self.__dict__.pop(key, None)
                
        self._defaults = self._defaultlayouts[kind]
        settingsBase._set_defaults(self)
        
    def _set_globaldefaults(self):
        
        settingsBase._set_defaults(self, self._globaldefaults)
    
    def __setattr__(self, attr, value):
        """
        Override the default __setattr__ so that changing self.kind causes
        the defaults to be set for that kind
        """
        self.__dict__[attr] = value
        
        if attr == 'kind':
            
            if value not in self._defaultlayouts.keys():
                
                raise ValueError, 'Unrecognized kind {0}'.format(value)
                
            self._set_defaults()
            
        elif attr == '_defaultlayouts':
            
            self._description = '{0}\nAvailable kinds:\n{1}'\
                .format(self._header, self._defaultlayouts.keys())
                
    def print_defaults(self):
        
        for layout in self._defaultlayouts.keys():
            
            line = 'kind: {}'.format(layout)
            if layout == self.kind:
                line += '   (CURRENT KIND)'
            print line
            print repr_defaults(self._defaultlayouts[layout], indent=2)
        
        print 'Settings for all kinds:'
        print repr_defaults(self._globaldefaults, indent=2)
    
class filenames(settingsBase):
    """
    Filenames
    
    To echo settings:
        filenames()
    To access settings:
        filenames.settingname
    """
    _defaults = [
    ('IC_file_name', 'IC.p', 'Initial conditions filename'), 
    ('snapshotName', 'snapshot.std', 'Filename to save tipsy snapshot to'),
    ('paramName', 'snapshot.param', 'Filename to save ChaNGa .param file to.'),
    ('directorName', 'snapshot.director', 'Default .director filename'), 
    ('settings_filename', 'IC_settings.p', 'settings file name')
    ]
    def __init__(self):
        
        settingsBase.__init__(self, descr='Filenames')
        self._set_defaults()
        
        
class physical(settingsBase):
    """
    Defines default physical parameters
    
    To echo settings:
        physical()
    To access settings:
        physical.settingname
    
    """
    _defaults = [
    ('m', SimArray(2.00132,'m_p'), 'Molecular mass of the gas'),
    ('M', SimArray(0.33 , 'Msol'), 'Mass of the star (or total mass of '\
         'binary system).'),
    ('T0', SimArray(332.406,'K'), 'Temp at R0 for powerlaw profile: T(r) = '\
         'T0(r/r0)^Tpower'),
    ('r0', SimArray(0.5,'au'), 'R0 for temp powerlaw: T(r) = T0(r/r0)^Tpower'),
    ('Tpower', -0.59, 'Powerlaw for temperature profile'),
    ('Tmin',  SimArray(0, 'K'), 'Minimum temperature'),
    ('Tmax',  SimArray(np.inf, 'K'), 'Maximum temperature'),
    ('kind', 'powerlaw', 'Temp profile kind.  See calc_temp.py'),
    ('eos',  'isothermal', 'Equation of state (isothermal or adiabatic)'),
    ('gamma',  1.4, 'Adiabatic index'),
    ('starMode',  'single', 'Star system type: single or binary' ),
    ('binsys',  None, 'Binary orbital parameters: to be set at runtime, see '\
         'binary.py'),
    ]
    def __init__(self, **kwargs):
        
        settingsBase.__init__(self, descr='General physical parameters:')
        self._set_defaults(**kwargs)
        
        
class sigma(kindSettingsBase):
    """
    Settings for generating a surface density profile
    
    To echo:
        sigma()
    To access settings:
        sigma.settingname
    """
    _defaultlayouts = {}
    _defaultlayouts['powerlaw'] = [
    ('Rd', SimArray(1.0,'au'), 'Disk radius - where exponential cutoff is '\
         'applied'),
    ('rin', 0.5, 'Fraction of disk radius where interior cutoff ends.  The '\
         'interior profile goes to 0 significantly inside of this'),
    ('rmax', 2.3, 'Maximum multiple of Rd to calculate surface density for'),
    ('cutlength', 0.3, 'Steepness of the exterior cutoff'),
    ('Qmin', 1.5, 'Approximate minimum Toomre Q: defines disk mass.  The real '\
         'value will differ'),
    ('n_points', 1000, 'Number of radial points to calculate surface density at'),
    ('power', -1., 'Powerlaw for Sigma = Sigma_0 (R/R0)^power')
    ]
    _defaultlayouts['mqws'] = [
    ('rin', 4.0, 'Dimensionless interior cutoff radius'),
    ('rout', 20.0, 'Dimensionless exterior cutoff radius'),
    ('rmax', None, 'Maximum radius'),
    ('power', -1, 'Powerlaw.  Note this multiplies a profile'),
    ('m_disk', SimArray(0.1, 'Msol'), 'Disk mass'),
    ('n_points', 1000, 'Number of radial points to calculate surface density at'),
    ('Qmin', 1.5, 'Approximate minimum Toomre Q.')
    ]
    _defaultlayouts['viscous'] = [
    ('Rd', SimArray(1.0, 'au'), 'Disk Radius'),
    ('rin', 0.1, 'Dimensionless interior cutoff.  Note the surface density '\
         'goes to zero significantly interior to this'),
    ('rmax', 2.0, 'Maximum multiple of Rd before cutoff'),
    ('power', -1, 'Powerlaw (note this also multiplies an exponential)'),
    ('m_disk', SimArray(0.1, 'Msol'), 'Disk mass'),
    ('n_points', 1500, 'Number of radial points to calculate surface density at'),
    ('gamma', 0.9, 'Power to set')
    ]
    _defaultlayouts['none'] = []
    _globaldefaults = [
    ('innercut', None, 'For a HARD (step function) interior cutoff applied'\
     'at this radius (dimensionless, multiple of Rd)'),
    ('outercut', None, 'For a HARD (step function) exterior cutoff applied'\
     'at this radius (dimensionless, multiple of Rd)')
    ]
    
    def __init__(self, kind='powerlaw'):
        if kind is None:
            
            kind = 'none'
            
        kindSettingsBase.__init__(self, 'Sigma profile parameters')
        self.kind = kind

class rho_calc(settingsBase):
    """
    Settings for calculating rho(z,r)
    
    To echo settings:
        rho_calc()
    To access settings:
        rho_calc.settingname
    """
    _defaults = [
    ('nr', None, 'The number of radial data for rho calc.  If None, '\
         'automatically determined'), 
    ('nz', 1000, 'The number of vertical points to calculate rho(z,r) at'), 
    # The maximum z (assumed to be au) to calculate rho(z,r) at.
    # If None (default, safe) zmax is automatically estimated at each rbin
    ('zmax', None, 'The maximum z to calculate rho.  If not, automatically'\
         'determined'),
    # 
    ('min_r_bins', 150, 'Minimum number of rbins to use.  More will be used '\
         'if necessary'),
    ('r_bin_tol', 1e-3, 'Accuracy requirement for r bins')
    ]
    def __init__(self):
        
        settingsBase.__init__(self, 'Rho calculation settings')
        self._set_defaults()
        
        
class pos_gen(settingsBase):
    """
    Settings for generating random positions [STEP 2]
    
    To echo settings:
        pos_gen()
    To access settings:
        pos_gen.settingname
    """
    _defaults = [
    ('nParticles', 40411, 'Number of particles'),
    ('method', 'grid', \
     'Method for generating (r, theta).  Can be "grid" or "random"')
    ]
    def __init__(self):
        
        settingsBase.__init__(self, 'Position generator settings:')
        self._set_defaults()
        
class snapshot(settingsBase):
    """
    Settings for generating tipsy files (includes particle mass, temp, and vel.)
    [STEP 3]
    """
    _defaults = [
    ('mScale', 1.0, 'For growMass runs, factor (0 to 1) by which to scale '\
         'particle masses for ICs'),
    ('nSmooth', 32, 'number of neighbors used during ChaNGa calculations'), 
    ('extraparams', {}, 'Extra runtime parameters for ChaNGa during vel calc'), 
    ('metals', 1.0, 'metals')
    ]
    def __init__(self, nParticles=0):
        
        message = 'Tipsy snapshot generator settings'
        message += '\n extraparams is a dict of extra parameters for ChaNGa'
        message += '\n e.g. extraparams = {"bDoGas": 1}'
        settingsBase.__init__(self, message)
                
class changa_run(settingsBase):
    """
    Settings for running ChaNGa (needed when calculating velocities, eps, so on)
    """
    _defaults = [
    ('preset', 'local', 'Run configuration to use.  '\
         'See ICgen_utils.changa_command for options'),
    ('changa_args', '', 'Additional arguments for all ChaNGa calls'),
    ('runner_args', '', 'Additional arguments for all runner '\
         '(mpirun, charmrun, ...) calls'),
    ('verbose', True, 'Display ChaNGa output'),
    ('logfile_name', None, 'Save ChaNGa output to log file')
    ]
    def __init__(self, **kwargs):
        
        settingsBase.__init__(self, 'ChaNGa run options')
        self._set_defaults(**kwargs)
        
class settings(settingsBase):
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
        
        settingsBase.__init__(self, descr='IC settings')
        
        if settings_filename is None:
            # Load the defaults
            self.filenames = filenames()
            self.sigma = sigma(kind)           
            self.rho_calc = rho_calc()
            self.pos_gen = pos_gen()
            self.snapshot = snapshot(self.pos_gen.nParticles)
            self.physical = physical(kind=kind)
            self.changa_run = changa_run()
            
        else:
            
            self.load(settings_filename)
        
    def save(self,settings_filename = None):
        """
        Save settings to settings_filename.  If settings_filename is None
        (or is not set), settings_filename will be taken as the filename of
        the most recent save/load
        """
        
        if settings_filename is None:
            # Try filename stored in settings object (self)            
            if not (hasattr(self, 'filenames') \
                    and hasattr(self.filenames,'settings_filename')):
                
                raise RuntimeError,'No settings_filename defined.  Cannot save'
                
            else:
                
                settings_filename = self.filenames.settings_filename
                
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
                
                
def repr_defaults(x, indent=0):
    """
    Prints a defaults list.  Defaults lists are a list of length 3 tuples:
        [ ('variable_name', default_val, 'Doc string for variable'),
          ('variable_name2', default_val2, 'Other doc string') ]
    """
    if len(x) == 0:
        return ''
    strings = []
    for key, default, doc in x:
        
        if key[0] != '_':
            
            
            if pynbody.units.has_unit(default):
                
                default_str = '{0} {1}'.format(default,default.units)
                
            else:
                
                default_str = '{0}'.format(default)
                
            line = ['{}'.format(key), default_str, doc]
            strings.append(line)
    
    # Now make them column formatted
    colwidth0 = max(len(line[0]) for line in strings)
    colwidth1 = max(len(line[1]) for line in strings)
    string = ''
    for line in strings:
        string += indent*' ' + line[0].rjust(colwidth0) + ': ' \
            + line[1].ljust(colwidth1) + ' - ' + line[2] + '\n'
    return string
        

def repr_settings(setting):
    
    header = getattr(setting, '_description', None)
    string = '------------------------\n'
    string += '{0}\n\n'.format(header)
    
    for key,val in setting.__dict__.iteritems():
            
            if key[0] == '_':
                
                pass
            
            elif isinstance(val, settingsBase):
                
                string += '{0}\n'.format(val)
                
            elif pynbody.units.has_units(val):
                
                string += '{0} : {1} {2}\n'.format(key,val,val.units)
                
            else:
                
                string += '{0} : {1}\n'.format(key,val)
    
    return string

    
def print_settings(setting):
    
    print repr_settings(setting)

def _colprint(data):
    col_width = max(len(word) for row in data for word in row) + 2  # padding
    for row in data:
        print "".join(word.ljust(col_width) for word in row)