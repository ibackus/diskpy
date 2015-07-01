"""
author: @dflemin3
Module for Binary star class.  Holds the Cartesian position/velocity coordinates and Kepler orbital elements for the binary in the 
reduced mass frame.  Can be initialized with orbital elements (preferred) or Cartesian position/velocity or a tipsy format snapshot
that is readable by pynbody.  Just need to pass a string to let it know which input type it's getting.

Rest of docstring: todo
"""

#Import my binary star module that continue relevant routines
import AddBinary

class Binary:
	"""
	Defines the binary star class Binary.	
	Binary star class used for storing positions/velocities and Keplerian Orbital elements for the system.
	
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
	
	def __init__(self,X,flag):
		if flag == "Cartesian" or "cartesian":
			pass
		elif flag == "Kepler" or "kepler":
			pass
		elif flag == "Snapshot" or "snapshot":
			pass
		else:
			print "Invalid input data type flag: %s.  All params set to 0 by default." % flag
			
			
	