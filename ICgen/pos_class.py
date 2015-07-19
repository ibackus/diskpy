# -*- coding: utf-8 -*-
"""
Defines a function to randomly generate particle positions according to 
the desired surface density profile (sigma vs r) and the vertical profile
(rho vs r,z).

Created on Mon Jan 27 18:48:04 2014

@author: ibackus
"""

__version__ = "$Revision: 1 $"
# $Source$

__iversion__ = int(filter(str.isdigit,__version__))

# External packages
import pynbody
SimArray = pynbody.array.SimArray
import numpy as np

# ICgen packages
import isaac
import ICgen_utils

class pos:
    """
    position class.  Generates particle positions from rho and sigma
    
    USAGE:
    # method = 'grid' or 'random'    
    pos = pos_class.pos(ICobj, method)
    
    ICobj should be an initial conditions object (ICgen.IC) with rho already
    calculated.

    """
    
    def __init__(self, ICobj, method = None, generate=True, seed=None):
        
        self._seed = seed
        # Set version
        self.__version__ = __iversion__
        # Link to parent initial conditions object
        self._parent = ICobj
        # Check that sigma and rho have been generated
        if not hasattr(ICobj, 'rho'):
            
            raise NameError,'rho could not be found in the IC object'
        
        if not hasattr(ICobj,'sigma'):
            
            raise NameError,'sigma could not be found in the IC object'
            
        if method == None:
        
            self.method = ICobj.settings.pos_gen.method
            
        else:
            
            self.method = method
            # Update settings in ICobj
            ICobj.settings.pos_gen.method = method
            
        self.nParticles = ICobj.settings.pos_gen.nParticles
        print 'Generating {0} particle positions using method: {1}'.format(\
        self.nParticles, self.method)
        
        # Generate positions
        self._generate_r()
        self.xyz = SimArray(np.zeros([self.nParticles, 3], dtype=np.float32), self.r.units)
        self._generate_z()
        self._generate_theta()
        self._cartesian_pos()
        
        # To save on memory, delete theta.  It can be re-calculated later
        # if absolutely needed
        del self.theta
        
    def __getstate__(self):
        """
        This is required to make the object pickle-able
        """
        
        # Define a dictionary containing everything needed.  
        # Ignore self.parent
        state = self.__dict__.copy()
        state.pop('_parent', None)
        
        # Now handle the possibly large arrays (too large to pickle)
        for key,val in state.iteritems():
            
            if isinstance(val, np.ndarray):
                
                state[key] = ICgen_utils.listify(val, 1001)
                
        return state
        
    def __setstate__(self, d):
        """
        This is required to make the object un-pickleable
        """
        for key, val in d.iteritems():
            
            if isinstance(val, ICgen_utils.larray):
                
                d[key] = val.delistify()
                
        self.__dict__ = d
        
    
    def _generate_r(self):
        """
        Generate radial positions
        """
        
        print 'Generating r positions'
        cdf_inv_r = self._parent.sigma.cdf_inv
        
        if self.method == 'grid':
            
            # Generate linearly increasing values of m, using 2 more than
            # necessary to avoid boundary issues
            m = np.linspace(0,1,self.nParticles + 2)
            # Calculate r from inverse CDF
            r = cdf_inv_r(m[1:-1]).astype(np.float32)
            # Assign output
            self.r = r
            
        if self.method == 'random':
            
            np.random.seed(self._seed)
            m = np.random.rand(self.nParticles)
            r = cdf_inv_r(m).astype(np.float32)
            self.r = r
            
    def _generate_z(self):
        """
        Generate z positions
        """
        
        print 'Generating z positions'
        # The inverse CDF over z as a function of r
        cdf_inv_z = self._parent.rho.cdf_inv
        # Random numbers between 0 and 1
        np.random.seed(self._seed)
        m = np.random.rand(self.nParticles)
        # Calculate z
        z = cdf_inv_z(m, self.r)
        # Randomly select sign of z
        z = z * np.random.choice(np.array([-1,1]), self.nParticles)
        # Assign output
        self.xyz[:,2] = z
        
    def _generate_theta(self):
        """
        Generate angular positions
        """
        
        nParticles = self.nParticles
        
        if self.method == 'grid':
            
            r = self.r
            
            dtheta = np.sqrt(2*np.pi*(1 - r[0:-1]/r[1:]))
            dtheta = isaac.strip_units(dtheta)
            theta = np.zeros(nParticles)
            
            for n in range(nParticles - 1):
                
                # NOTE: it's import to subtract (not add) dtheta.  The particles
                # will be moving counter-clockwise.  To prevent the particle
                # spirals from kinking, the particle spirals must go out
                # clockwise
                theta[n+1] = theta[n] - dtheta[n]
                
            self.theta = theta
            
        if self.method == 'random':
            
            np.random.seed(self._seed)
            theta = 2*np.pi*np.random.rand(nParticles)
            self.theta = theta
            
    def _cartesian_pos(self):
        """
        Generate x,y
        """
        
        r = self.r
        theta = self.theta
        self.xyz[:,0] = r*np.cos(theta)
        self.xyz[:,1] = r*np.sin(theta)
