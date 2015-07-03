"""
author: @dflemin3
Module for Binary star class.  Holds the Cartesian position/velocity coordinates and Kepler orbital elements for the binary in the 
reduced mass frame.  Can be initialized with orbital elements (preferred) or Cartesian position/velocity or a tipsy format snapshot
that is readable by pynbody.  Just need to pass a string to let it know which input type it's getting.

Rest of docstring: todo
"""
import numpy as np

#Import my binary star module that continue relevant routines
import AddBinary

class Binary(object):
	"""
	Defines the binary star class Binary.	
	Binary star class used for storing positions/velocities and Keplerian Orbital elements for the system.
	
	Input Units: Lengths = AU, Velocities = km/s (*29.8 == VEL_UNIT), angles in degrees	
	Output Units: Same as input	
	
	Initializing:
	
	With Cartesian:
		Binary(X,m1,m2,"Cartesian") where X is a list or numpy array of the form [x,y,z,vx,vy,vz] and m1, m2
		are the masses of the primary and secondary (M_sol), respectively
	With Kepler Orbital Elements:
		Binary(X,m1,m2,"Kepler") where X is e, a, i, Omega, w, nu which are the eccentricity, semimajor axis (AU),
		inclination (degrees), Longitude of Ascending Node (degrees), argument of periapsis (degrees), and
		the true anomaly (degrees).  Note the units!  m1, m2 are defined above.
	With Snapshot:
	
	
	"""	
	
	def __init__(self,X,m1,m2,state):
		
		if state == "Cartesian" or state == "cartesian":
			#Ensure input is proper
			assert (len(X) == 2), "len(Input Array) != 2. len = %d." % len(X)				
			
			self.r = X[0]
			self.v = X[1]
			self.m1 = m1
			self.m2 = m2
			self.state = state
			self.reshapeData()
			
		elif state == "Kepler" or state == "kepler":
			#Ensure input is proper
			assert (len(X) == 6), "len(Input Array) != 6. len = %d." % len(X)	
		
			self.assignOrbElems(X)
			self.m1 = m1
			self.m2 = m2
			self.state = state
		elif state == "Snapshot" or state == "snapshot":
			pass
		else:
			print "Invalid input data type state: %s.  All params set to 0 by default." % state
		
	#end function
		
	#Member Functions		
		
	def assignOrbElems(self,X):
		"""
		Given array of orbital elements, set them as class parameters.  Intended as a shorthand function.
		
		Input:
		X is e, a, i, Omega, w, nu which are the eccentricity, semimajor axis (AU),
		inclination (degrees), Longitude of Ascending Node (degrees), argument of periapsis (degrees), and
		the true anomaly (degrees).
			
		Output:
		None
		"""
		self.e = X[0]
		self.a = X[1]
		self.i = X[2]
		self.Omega = X[3]
		self.w = X[4]
		self.nu = X[5]

	#end function		
		
	def reshapeData(self):
		"""
		Ensure data, specifically r and v arrays, are of the proper shape for further manipulation.
		"""
		self.r = self.r.reshape((1,3))
		self.v = self.v.reshape((1,3))
		
	#end function

	def computeOrbElems(self):
		"""
		Compute the Kepler orbital elements.
		
		Input: Self (must have r, v set and in proper units)
		
		Output: sets and returns e,a,i,... 
		"""
		#Compute orbital elements from binary center of mass frame Cartesian coordinates
		zero = np.asarray([0,0,0])
		oe = AddBinary.calcOrbitalElements(self.r,zero,self.v,zero,self.m1,self.m2)

		#Set orbital elements, return them as well		
		self.assignOrbElems(oe)
		return oe
		
	#end function		
		
	def computeCartesian(self):
		"""
		Compute the Cartesian position and velocity in the reduced mass frame.
		"""
		assert (self.state == "Kepler" or self.state == "kepler")
		zero = np.asarray([0.,0.,0.])		
		M = AddBinary.calcMeanAnomaly(self.r,zero,self.v,zero,self.m1,self.m2)
		self.r, self.v = AddBinary.keplerToCartesian(self.a,self.e,self.i,self.Omega,self.w,M,self.m1,self.m2)

	#end function		
		
	def generateICs(self):
		"""
		From Kepler orbital elements, compute the initial position, velocities for two stars in ChaNGa-friendly units.
		"""
		#Define required array of 0s		
		zero = np.asarray([0.,0.,0.])
		
		
		if self.state == "Kepler" or self.state == "kepler":
			M = AddBinary.trueToMean(self.nu,self.e)
			x1,x2,v1,v2 = AddBinary.initializeBinary(self.a,self.e,self.i,self.Omega,self.w,M,self.m1,self.m2)
			
		elif self.state == "Cartesian" or self.state == "cartesian":
			M = AddBinary.calcMeanAnomaly(self.r,zero,self.v,zero,self.m1,self.m2)
			x1,x2,v1,v2 = AddBinary.reduceToPhysical(self.r,self.v,self.m1,self.m2)
			
		else:
			pass
		
		return x1, x2, v1, v2
		
	#end function
		
		
		
		
		
		
		