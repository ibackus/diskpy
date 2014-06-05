# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 12:48:33 2014

@author: ibackus
"""

# DEBUGGING
import sys

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
import calc_rho
import pos_gen_grid as pos_gen
import calc_temp
import calc_sigma
import pos_class
import make_snapshot
import isaac

import ICgen_settings
import make_sigma

ICgenDir = os.path.dirname(os.path.realpath(__file__))
param_filepath = os.path.join(ICgenDir, 'default.param')
param_default = isaac.configparser(param_filepath, ftype='param')

class IC:
    
    def __init__(self):
        
        # Initialize
        # Load up default settings
        self.settings = ICgen_settings.settings()
        # Add modules/attributes
        self.T = calc_temp.T(self)
        self.maker = maker(self)
        self.add = add(self)
        
        def saver(filename = None):
            """
            A wrapper for ICgen.save
            """
            save(self, filename)
            
        self.save = saver
        
def save(ICobj, filename=None):
    
    if filename is None:
        
        filename = ICobj.settings.filenames.IC_file_name
    
    save_dict = {}
    # --------------------------------------------------
    # GET SETTINGS
    # --------------------------------------------------
    save_dict['settings'] = ICobj.settings
    
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
        kind = sigma.kind
        save_dict['settings'].sigma.kind = kind
        
        # Change how sigma is saved depending on the method (kind) used to
        # generate sigma
        if kind == 'file':
            # Update save dictionary
            save_dict['sigma'] = sigma.input_dict

    # --------------------------------------------------
    # Prepare pos if possible
    # --------------------------------------------------
    if hasattr(ICobj, 'pos'):
        
        save_dict['pos'] = ICobj.pos
        
    # --------------------------------------------------
    # Prepare param if possible
    # --------------------------------------------------
    if hasattr(ICobj, 'snapshot_param'):
        
        save_dict['snapshot_param'] = ICobj.snapshot_param
        param_name = ICobj.settings.filenames.paramName
        isaac.configsave(ICobj.snapshot_param, param_name)
        print 'param file saved to {}'.format(param_name)
        
    # --------------------------------------------------
    # SAVE
    # --------------------------------------------------
    # Save snapshot if possible
    if hasattr(ICobj, 'snapshot'):
        
        fmt = pynbody.tipsy.TipsySnap
        fname = ICobj.settings.filenames.snapshotName
        save_dict['snapshotName'] = fname
        ICobj.snapshot.write(fmt = fmt, filename = fname)
        
    # Save the save dictionary
    pickle.dump(save_dict,open(filename,'wb'))
    print 'Initial conditions saved to {}'.format(filename)
    
def load(filename):
    
    # Initialize a blank IC object
    ICobj = IC()
    
    # Load everything available from filename
    input_dict = pickle.load(open(filename,'rb'))
    
    # Parse the input dictionary
    if 'settings' in input_dict:

        print 'loading settings'        
        ICobj.settings = input_dict['settings']
        
    if 'rho' in input_dict:
        
        print 'loading rho'
        ICobj.add.rho(input_dict['rho'])
        
    if 'sigma' in input_dict:
        
        print 'loading sigma'
        ICobj.add.sigma(input_dict['sigma'])
        
    if 'pos' in input_dict:
        
        print 'loading pos'
        ICobj.pos = input_dict['pos']
        
    if 'snapshotName' in input_dict:
        
        print 'loading snapshot'
        fname = input_dict['snapshotName']
        ICobj.snapshot = pynbody.load(fname)
        
    if 'snapshot_param' in input_dict:
        
        print 'loading param'
        ICobj.snapshot_param = input_dict['snapshot_param']

    return ICobj
            
        
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
        # Create rho object (includes a spline interpolation)
        rho_binned = rho_dict['rho']
        z_bins = rho_dict['z']
        r_bins = rho_dict['r']
        
        self._parent.rho = calc_rho_zr.rho_from_array(self._parent, rho_binned, z_bins, r_bins)
        
        print 'rho stored in <IC instance>.rho'
        
    def sigma(self, sigma_input):
        """
        Generates a sigma object and stores it in ICobj.rho
        
        IF ICobj.settings.sigma.kind = 'file':
            sigma_input should be a dictionary containing 'sigma' and 'r'
        
        Example:
        
        sigma_dict = pickle.load(open('sigmafile.p','rb'))
        IC.settings.sigma.kind = 'file'
        IC.add.sigma(sigma_dict)
        """
        sigma_out = make_sigma.sigma_gen(self._parent.settings, sigma_input)
        # Assign to parent ICobj:
        self._parent.sigma = sigma_out
        print 'Sigma stored in <IC instance>.sigma'
        

class maker:
    """
    A Wrapper containing various functions for generating initial conditions.
    Outputs of the functions are saved to the IC object.  The IC object is 
    referenced as self._parent.  So to access temperature, simply call
    self._parent.T(r)
    """
    
    def __init__(self, ICobj):
        
        self._parent = ICobj
        
    def sigma_gen(self):
        """
        A Wrapper for make_sigma.sigma_gen
        See make_sigma.sigma_gen for documentation
        
        Upon executing, generates sigma, pdf, and cdf_inv according to
        settings.sigma.kind and saves to ICobj
        
        USAGE:
        
        ICobj.maker.sigma_gen()
        """
        # Generate sigma
        sigma = make_sigma.sigma_gen(self._parent.settings)
        # Copy sigma to the parent (IC) object
        self._parent.sigma = sigma
        
        print 'Sigma stored in <IC instance>.sigma'
        
    def rho_gen(self):
        """
        A wrapper for calc_rho_zr.
        
        Upon executing, generates rho and rho cdf inverse
        """
        
        # Check that sigma has been generated
        if not hasattr(self._parent, 'sigma'):
            
            raise RuntimeError,'Must load/generate sigma before calculating rho'
            
        # Numerically calculate rho(z,r) for a given sigma.  rho(z,r)
        # obeys vertical hydrostatic equilibrium (approximately)
        rho_array, z, r = calc_rho_zr.rho_zr(self._parent)
        # Create a complete rho object.  Includes rho spline and CDF inverse
        rho = calc_rho_zr.rho_from_array(self._parent, rho_array, z, r)
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
        
        # Generate snapshot
        snapshot, snapshot_param = make_snapshot.snapshot_gen(self._parent)
        # Save to ICobj
        self._parent.snapshot = snapshot
        self._parent.snapshot_param = snapshot_param