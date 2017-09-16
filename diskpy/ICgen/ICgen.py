# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 12:48:33 2014

@author: ibackus
"""

__version__ = "$Revision: 4 $"
# $Source$

__iversion__ = int(filter(str.isdigit,__version__))

# External modules
import pynbody
SimArray = pynbody.array.SimArray
import numpy as np
import os
import cPickle as pickle
from warnings import warn

# ICgen modules
import calc_temp
import pos_class
import make_snapshot
import make_snapshotBinary
import make_snapshotSType
import make_snapshotSTypeSupplied
import ICgen_settings
import ICgen_utils
import make_sigma
import sigma_profile

from rhosolver import rhosolver, loadrho

# diskpy modules
from diskpy.utils import match_units, configsave
from diskpy import global_settings

# Initial stuff
ICgenDir = os.path.dirname(os.path.realpath(__file__))

class IC:
    """
    Defines the IC class.
    
    INITIALIZING NEW INITIAL CONDITIONS

    # Initialize a blank IC object (with no surface density profile yet made)

    IC = ICgen.IC()
    
    # Initialize an IC object and surface density profile using default settings
    
    IC = ICgen.IC(profile_kind='powerlaw')
    IC = ICgen.IC(profile_kind='MQWS')
    
    # Initialize IC object from 1-D SimArrays r, sigma (surface density)
    
    IC = ICgen.IC(r, sigma)
    
    Optionally, the CDF for the surface density profile can be supplied to
    speed up generation of sigma.  To do that:
    
    IC = ICgen.IC(r, sigma, CDF)
    
    Or, the input can be a the filename of a pickled dictionary containing
    'r', 'sigma', and optionally 'CDF'
    
    Settings can also be entered manually if needed
    
    settings = pickle.load(open('settings.p','r'))
    IC = ICgen.IC(settings = settings)
    
    GENERATING INITIAL CONDITIONS
    
    """
    
    def __init__(self, r=None, sigma=None, CDF=None, profile_kind=None, settings=None):
        
        if isinstance(r,str):
            # Load the pickled sigma dictionary
            sig_dict = pickle.load(open(r,'r'))
            sigma = sig_dict['sigma']
            r = sig_dict['r']
            
            if 'CDF' in sig_dict:
                
                CDF = sig_dict['CDF']
            
            
        # Initialize
        self.__version__ = __iversion__
        
        if settings is None:
            # Load up default settings
            self.settings = ICgen_settings.settings(kind=profile_kind)
            
        else:
            
            self.settings = settings
        
        # Try to update profile_kind (might be set in settings, or from defaults)
        try:
            
            profile_kind = self.settings.sigma.kind
            
        except:
            
            pass
        
        # Add modules/attributes
        self.T = calc_temp.T(self)
        self.maker = maker(self)
        self.add = add(self)
        
        # Generate sigma spline interpolation
        if profile_kind == 'none':
            profile_kind = None
        if (r is not None) and (profile_kind is not None):
            
            self.maker.sigma_gen(r, sigma, CDF)
        
        # Define a saving function        
        def saver(filename = None):
            """
            A wrapper for ICgen.save
            """
            save(self, filename)
            
        self.save = saver
        
    def Qest(self, r=None):
        """
        Estimate Toomre Q at r (optional) for ICs, assuming omega=epicyclic
        frequency.  Ignores disk self-gravity
        """
        return Qest(self, r)
        
    def generate(self, restart=False):
        """
        Runs through all the steps to generate a set of initial conditions
        
        IF restart=True, it picks up at the last completed step
        """
        if restart:
            
            # Find the last completed step
            if hasattr(self, 'pos'): initial_step = 3
            elif hasattr(self, 'rho'): initial_step = 2
            elif hasattr(self, 'sigma'): initial_step = 1
            else: initial_step = 0
        
        else:
            
            initial_step = 0
        
        self.save()
        
        if initial_step <= 0:
            
            # Generate sigma
            self.maker.sigma_gen()
            self.save()
        
        if initial_step <= 1:
            
            # Generate rho
            self.maker.rho_gen()
            self.save()
            
        if initial_step <= 2:
            
            # Generate positions
            self.maker.pos_gen()
            self.save()
            
        if initial_step <= 3:
            
            # Generate snapshot
            self.maker.snapshot_gen()
            self.save()
            
def Qest(ICobj, r=None):
    """
    Estimate Toomre Q at r (optional) for ICs, assuming omega=epicyclic
    frequency.  Ignores disk self-gravity
    """
    if not hasattr(ICobj, 'sigma'):
        
        raise ValueError, 'Could not find surface density profile (sigma)'
        
    G = SimArray(1.0, 'G')
    kB = SimArray(1.0, 'k')
    
    if r is None:
        
        r = ICobj.sigma.r_bins
        
    sigma = ICobj.sigma(r)
    T = ICobj.T(r)
    M = ICobj.settings.physical.M
    m = ICobj.settings.physical.m
    M = match_units(M, 'Msol')[0]
    m = match_units(m, 'm_p')[0]
    gamma = ICobj.settings.physical.gamma_cs()
    
    Q = np.sqrt(M*kB*T*gamma/(G*m*r**3))/(np.pi*sigma)
    Q.convert_units('1')
    
    return Q
        
        
def save(ICobj, filename=None):
    
    if filename is None:
        
        filename = ICobj.settings.filenames.IC_file_name
    
    save_dict = {}
    
    if hasattr(ICobj, '__version__'):
        
        save_dict['version'] = ICobj.__version__
        
    # --------------------------------------------------
    # GET SETTINGS/save a copy
    # --------------------------------------------------
    settings_name = os.path.splitext(filename)[0] + '_settings.p'
    save_dict['settings'] = settings_name
    ICobj.settings.save(settings_name)
    
    # --------------------------------------------------
    # Prepare rho, if available
    # --------------------------------------------------
    if hasattr(ICobj, 'rho'):
        
        rho = ICobj.rho
        # Generate a dictionary containing rho_binned, z_bins, r_bins
        rho_dict = {\
        'rho': rho.rho_binned,\
        'z': rho.z_bins,\
        'r': rho.r_bins}
        # Update save dictionary
        save_dict['rho'] = rho_dict
    
    # --------------------------------------------------
    # Prepare sigma, if available
    # --------------------------------------------------
    if hasattr(ICobj, 'sigma'):
        
        sigma = ICobj.sigma
        # Update save dictionary
        save_dict['sigma'] = sigma.input_dict
        save_dict['CDF'] = sigma._CDF

    # --------------------------------------------------
    # Prepare pos if possible
    # --------------------------------------------------
    if hasattr(ICobj, 'pos'):
        
        save_dict['pos'] = ICobj.pos
        
    # --------------------------------------------------
    # Prepare param/dict if possible
    # --------------------------------------------------
    if hasattr(ICobj, 'snapshot_param'):
        
        save_dict['snapshot_param'] = ICobj.snapshot_param
        param_name = ICobj.settings.filenames.paramName
        configsave(ICobj.snapshot_param, param_name)
        print 'param file saved to {0}'.format(param_name)
        
    if hasattr(ICobj, 'snapshot_director'):
        
        save_dict['snapshot_director'] = ICobj.snapshot_director
        director_name = ICobj.settings.filenames.directorName
        configsave(ICobj.snapshot_director, director_name)
        print 'director file saved to {0}'.format(director_name)
        
    # --------------------------------------------------
    # SAVE
    # --------------------------------------------------
    # Save snapshot if possible
    if hasattr(ICobj, 'snapshot'):
        
        fmt = pynbody.tipsy.TipsySnap
        fname = ICobj.settings.filenames.snapshotName
        save_dict['snapshotName'] = fname
        
        # Sometimes, saving the snapshot once raises an error.  Saving again
        # can fix this for some reason
        
        try:
            
            ICobj.snapshot.write(fmt = fmt, filename = fname)
            
        except ValueError:
            
            ICobj.snapshot.write(fmt = fmt, filename = fname)
        
    # Save the save dictionary
    pickle.dump(save_dict,open(filename,'wb'), protocol=2)
    print 'Initial conditions saved to {0}'.format(filename)        
    
def load(filename):
       
    # Load everything available from filename
    input_dict = ICgen_utils.pickle_import(filename, ICgenDir)
    
    # Get version/update IC if necessary
    if 'version' in input_dict:
        
        version = input_dict['version']
        
    else:
        
        version = 0
    
    _upgrade_version(input_dict, version)
    
    # Parse the input dictionary
    
    if 'sigma' in input_dict:
        
        # Load sigma stuff
        sigma = input_dict['sigma']['sigma']
        r = input_dict['sigma']['r']
        
        # Initialize ICobj
        if 'CDF' in input_dict:
            
            CDF = input_dict['CDF']
            ICobj = IC(r, sigma, CDF)
            
        else:
        
            ICobj = IC(r, sigma)
            
    else:
        # Initialize blank IC object, sigma has not been generated yet
        ICobj = IC()
        
    if 'settings' in input_dict:
        
        # 'settings' is the settings filename.  load that up
        ICobj.settings.load(input_dict['settings'])
        
    if 'rho' in input_dict:
        
        print 'loading rho'
        ICobj.add.rho(input_dict['rho'])
        
    if 'pos' in input_dict:
        
        print 'loading pos'
        ICobj.pos = input_dict['pos']
        ICobj.pos._parent = ICobj
        
    if 'snapshotName' in input_dict:
        
        print 'loading snapshot'
        fname = input_dict['snapshotName']
        
        try:
        
            ICobj.snapshot = pynbody.load(fname)
            if 'snapshot_param' in input_dict:
        
                print 'loading param'
                ICobj.snapshot_param = input_dict['snapshot_param']
                
            if 'snapshot_director' in input_dict:
                
                print 'loading director'
                ICobj.snapshot_director = input_dict['snapshot_director']
            
        except IOError:
            
            warn('Could not find snapshot ({0})'.format(fname))
            print 'Could not find snapshot ({0})'.format(fname)
        
    

    return ICobj
            
            
def _upgrade_version(IC_input, version):
    """
    Used for backwards compatibility.  If an initial conditions object was
    generated using an outdated version of ICgen, hopefully when it's being
    loaded it can be made up-to-date before initializing the ICobj
    
    version is the version of the IC_input
    """
    if version < 1:
        
        # Assume IC_input is a dictionary
        # Needs the changa_run settings
            
        # Initialize empty, up-to-date settings
        new_settings = ICgen_settings.settings()
        old_settings = IC_input['settings']
        
        if 'settings' in IC_input:
            
            # Load old settings, ignoring the snapshot gen settings
            for key in old_settings.__dict__.keys():
                
                if key != 'snapshot':
                    
                    new_settings.__dict__[key] = old_settings.__dict__[key]
                    
        # update the input dictionary
        IC_input['settings'] = new_settings
        
    if version < 2:
        
        # Add a .director file name to settings.
        
        warn('These ICs seem out of date.  Attempting to upgrade to current version')
        
        if not hasattr(IC_input['settings'].filenames, 'directorName'):
            
            paramName = IC_input['settings'].filenames.paramName
            directorName = os.path.splitext(paramName)[0] + '.director'
            IC_input['settings'].filenames.directorName = directorName
            
    if version < 4:
        
        IC_input['settings'] = IC_input['settings'].settings_filename
        
        
class add:
    """
    Contains modules to load data/parameters
    """
    
    def __init__(self, ICobj):
        
        self._parent = ICobj
        
    def rho(self,rho_dict):
        """
        Generates a rho object and stores it in ICobj.rho
        
        rho_dict should be a dictionary containing:
            'z':    1D array of z values
            'r':    1D array of r values
            'rho':  2D array of rho evaluated at z,r
            
        Exaple:
        
        rho_dict = pickle.load(open('rhofile.p', 'rb')) # Load up a rho dict
        ICobj.add.rho(rho_dict)     # create ICobj.rho
        """
        loadrho(self._parent, rho_dict)
                
        print 'rho stored in <IC instance>.rho'
        

class maker:
    """
    A Wrapper containing various functions for generating initial conditions.
    Outputs of the functions are saved to the IC object.  The IC object is 
    referenced as self._parent.  So to access temperature, simply call
    self._parent.T(r)
    """
    
    def __init__(self, ICobj):
        
        self._parent = ICobj
        
    def sigma_gen(self, r=None, sigma=None, CDF=None):
        """
        A Wrapper for make_sigma.sigma_gen
        See make_sigma.sigma_gen for documentation
        
        Upon executing, generates sigma, pdf, and cdf_inv and saves to ICobj
        
        USAGE:
        
        # Generate sigma object from r, sigma arrays
        sigma = ICobj.maker.sigma_gen(r, sigma)
        
        # Use pre-calculated CDF array
        sigma = ICobj.maker.sigma_gen(r, sigma, CDF)
        
        # Generate using a profile defined in sigma_profile.py
        ICobj.settings.sigma.kind = 'powerlaw'
        sigma = ICobj.maker.sigma_gen()
        
        r and sigma should be 1-D SimArrays.  sigma is the surface density
        evaluated at r
        """
        # Generate sigma
        if r is None:
            
            r, sigma = sigma_profile.make_profile(self._parent)
            
        sigma = make_sigma.sigma_gen(r, sigma, CDF)
        # Copy sigma to the parent (IC) object
        self._parent.sigma = sigma
        
        
        # and now, set up the interior temperature profile
        self._parent.T.setup_interior()
        
        print 'Sigma stored in <IC instance>.sigma'
        
    def rho_gen(self, **kwargs):
        """
        A wrapper for calc_rho_zr.
        
        Upon executing, generates rho and rho cdf inverse
        
        Optional kwargs can be supplied and passed to the rootfinder used
        (scipy.optimize.newton_krylov)
        """
        
        # Check that sigma has been generated
        if not hasattr(self._parent, 'sigma'):
            
            raise RuntimeError,'Must load/generate sigma before calculating rho'
            
        rho = rhosolver(self._parent)
        rho.solve(**kwargs)
        # Save to ICobj
        self._parent.rho = rho
        
        print 'rho stored in <IC instance>.rho'
        
    def pos_gen(self, method = None):
        """
        A wrapper for generating positions according to rho and sigma
        
        Initializes a pos object (see pos_class.py) and saves it to ICobj.pos
        
        IF called with method not set, the method used is:
            ICobj.settings.pos_gen.method
        """
        # Generate positions object
        pos = pos_class.pos(self._parent, method)
        # Save it to ICobj
        self._parent.pos = pos
        
    def snapshot_gen(self):
        """
        A wrapper for generating a tipsy snapshot from the initial conditions
        
        Uses make_snapshot.py
        """
        
        # Generate snapshot for either a single star or binary depending on 
        # IC.settings.physical.starMode
        if self._parent.settings.physical.starMode == "single":		
            
            snapshot, snapshot_param, snapshot_director = \
            make_snapshot.snapshot_gen(self._parent)
            
        elif self._parent.settings.physical.starMode == "binary":
            
            snapshot, snapshot_param, snapshot_director = \
            make_snapshot.snapshot_gen(self._parent)            
#            snapshot, snapshot_param, snapshot_director = \
#            make_snapshotBinary.snapshot_gen(self._parent)
            
        elif self._parent.settings.physical.starMode == "stype":
            
            snapshot, snapshot_param, snapshot_director = \
            make_snapshot.snapshot_gen(self._parent)
#            snapshot, snapshot_param, snapshot_director = \
#            make_snapshotSType.snapshot_gen(self._parent)
            
        elif self._parent.settings.physical.starMode == "stypeSupplied":
            
            snapshot, snapshot_param, snapshot_director = \
            make_snapshotSTypeSupplied.snapshot_gen(self._parent)
            
        else:
            
            print "Invalid starMode given in ICobj.  Assuming default single star."
            snapshot, snapshot_param, snapshot_director = \
            make_snapshot.snapshot_gen(self._parent)
            
        # Save to ICobj
        self._parent.snapshot = snapshot
        self._parent.snapshot_param = snapshot_param
        self._parent.snapshot_director = snapshot_director
        
