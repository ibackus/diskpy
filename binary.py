"""
author: @dflemin3 July 2015
Module for Binary star class.  Holds the Cartesian position/velocity coordinates and Kepler orbital elements for the binary in the 
reduced mass frame.  Can be initialized with orbital elements (preferred) or Cartesian position/velocity or a tipsy format snapshot
that was read in using pynbody.load("snapshot.name").  Just need to pass a string to let it know which input type it's getting.
"""
import numpy as np

#Import my binary star module that continue relevant routines
import AddBinary
import isaac

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
		X is a tipsy snapshot to be loaded.  Positions and orbital elements will be calculated from it.
	"""	
	
	def __init__(self,X,m1,m2,state):
		"""
		For a user-specified data type, initialize Cartesian positions in the center of mass frame and the
		associated Keplerian orbital elements for the binary system.
		"""
		
		if state == "Cartesian" or state == "cartesian":
			#Ensure input is proper
			assert (len(X) == 2), "Improper input. len(Input Array) != 2. len = %d.  State should be cartesian." % len(X)		
			self.state = state
			
			self.r = X[0]
			self.v = X[1]
			self.m1 = m1
			self.m2 = m2
			self.reshapeData()
			self.computeOrbElems()			
			
		elif state == "Kepler" or state == "kepler":
			#Ensure input is proper
			assert (len(X) == 6), "Improper input. len(Input Array) != 6. len = %d.  State should be kepler" % len(X)	
			self.state = state
		
			self.assignOrbElems(X)
			self.m1 = m1
			self.m2 = m2
			self.computeCartesian()
			self.reshapeData()
			
		elif state == "Snapshot" or state == "snapshot":		
			#Ensure input is indeed a tipsy snapshot of a binary system
			assert(len(X.stars) == 2), "Improper input.  Is this a tipsy snapshot of a circumbinary system?"			
			self.state = state			
			
			#Extract data from snapshot
			x1 = isaac.strip_units(X.stars[0]['pos'])
			x2 = isaac.strip_units(X.stars[1]['pos'])
			v1 = isaac.strip_units(X.stars[0]['vel'])
			v2 = isaac.strip_units(X.stars[1]['vel'])
			self.m1 = isaac.strip_units(X.stars[0]['mass'])
			self.m2 = isaac.strip_units(X.stars[1]['mass'])			
			
			#Compute position, velocity in center of mass frame then orbital elements
			self.r = x1 - x2
			self.v = v1 - v2
			self.reshapeData()
			self.computeOrbElems()
		
		else:
			print "Invalid input data type state: %s." % state
		
	#end function
		
	def __repr__(self):
		"""
		When invoked, return the orbital elements.
		"""		
		return "(%s,%s,%s,%s,%s,%s)" % (self.e,self.a,self.i,self.Omega,self.w,self.nu)
		
	#end function
		
	#Member Functions		
		
	def assignOrbElems(self,X):
		"""
		Given array of orbital elements, set them as class parameters.  Intended as a shorthand function.
		
		Input:
		X is e, a, i, Omega, w, nu which are the eccentricity, semimajor axis (AU),
		inclination (degrees), Longitude of Ascending Node (degrees), argument of periapsis (degrees), and
		the true anomaly (degrees).  Convert all values to float for consistency and as a sanity check.
			
		Output:
		None
		"""
		#Only assign data if it doesn't exist
		self.e = float(X[0])
		self.a = float(X[1])
		self.i = float(X[2])
		self.Omega = float(X[3])
		self.w = float(X[4])
		self.nu = float(X[5])

	#end function		
		
	def reshapeData(self):
		"""
		Ensure data, specifically r and v arrays, are of the proper shape for further manipulation.
		This is useful because sometimes they come in as a list, a (1,3) numpy array or a (3,) numpy array.
		Much easier to clean up upon initialization then have many checks in later functions.
		"""
		self.r = np.asarray(self.r).reshape((1,3))
		self.v = np.asarray(self.v).reshape((1,3))
		
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
		assert (self.state == "Kepler" or self.state == "kepler"), "Already have cartesian coords!"
		M = AddBinary.trueToMean(self.nu,self.e)
		self.r, self.v = AddBinary.keplerToCartesian(self.a,self.e,self.i,self.Omega,self.w,M,self.m1,self.m2)

	#end function		
		
	def generateICs(self):
		"""
		From Kepler orbital elements, compute the position, velocities for two stars in ChaNGa-friendly units.
		Called "generateICs" since I'll use this mainly to...generate initial conditions
		"""
		return AddBinary.reduceToPhysical(self.r,self.v,self.m1,self.m2)
				
	#end function