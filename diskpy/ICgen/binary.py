"""
author: @dflemin3 July 2015
Module for Binary star class.  Holds the Cartesian position/velocity coordinates and Kepler orbital elements for the binary in the
reduced mass frame.  Can be initialized with orbital elements (preferred) or Cartesian position/velocity or a tipsy format snapshot
that was read in using pynbody.load("snapshot.name").  Just need to pass a string to let it know which input type it's getting.
"""
import numpy as np

# Import my binary star module that continue relevant routines
import AddBinary
import pynbody
SimArray = pynbody.array.SimArray

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

    def __init__(self, X, m1, m2, state):
        """
        For a user-specified data type, initialize Cartesian positions in the center of mass frame and the
        associated Keplerian orbital elements for the binary system.
        """

        if state.lower() == "cartesian":
            # Ensure input is proper
            assert (len(
                X) == 2), "Improper input. len(Input Array) != 2. len = %d.  State should be cartesian." % len(X)
            self.state = state
            self.r = X[0]
            self.v = X[1]
            self.m1 = m1
            self.m2 = m2
            self.reshapeData()
            self.computeOrbElems()

        elif state.lower() == "kepler":
            # Ensure input is proper
            assert (len(
                X) == 6), "Improper input. len(Input Array) != 6. len = %d.  State should be kepler" % len(X)
            self.state = state
            self.assignOrbElems(X)
            self.m1 = m1
            self.m2 = m2
            self.computeCartesian()
            self.reshapeData()

        elif state.lower() == "snapshot":
            # Ensure input is indeed a tipsy snapshot of a binary system
            assert(len(
                X.stars) == 2), "Improper input.  Is this a tipsy snapshot of a circumbinary system?"
            self.state = state

            # Extract data from snapshot
            x1 = X.stars[0]['pos']
            x2 = X.stars[1]['pos']
            v1 = X.stars[0]['vel']
            v2 = X.stars[1]['vel']
            self.m1 = X.stars[0]['mass']
            self.m2 = X.stars[1]['mass']

            # Compute position, velocity in center of mass frame then orbital
            # elements
            self.r = x1 - x2
            self.v = v1 - v2
            self.reshapeData()
            self.computeOrbElems()

        else:
            print "Default Binary init"
            self.r = SimArray(np.zeros((1, 3)),'au')
            self.v = SimArray(np.zeros((1, 3)),'km s**-1')
            self.m1 = SimArray(0.0,'Msol')
            self.m2 = SimArray(0.0,'Msol')
            self.e = 0.0
            self.a = SimArray(0.0,'au')
            self.i = 0.0
            self.Omega = 0.0
            self.w = 0.0
            self.nu = 0.0
            self.state = 'Default'

    # end function

    def __repr__(self):
        """
        Return the orbital elements and the masses of the primary and secondary.
        """
        return "(%s,%s,%s,%s,%s,%s), mass: (%s,%s)" % (self.e,
                                                       self.a,
                                                       self.i,
                                                       self.Omega,
                                                       self.w,
                                                       self.nu,
                                                       self.m1,
                                                       self.m2)

    # end function

    # Member Functions

    def assignOrbElems(self, X):
        """
        Given array of orbital elements, set them as class parameters.  Intended as a shorthand function.

        Input:
        X is e, a, i, Omega, w, nu which are the eccentricity, semimajor axis (AU),
        inclination (degrees), Longitude of Ascending Node (degrees), argument of periapsis (degrees), and
        the true anomaly (degrees).  Convert all values to float for consistency and as a sanity check
        except for semimajor axis since that gets units we're interested in, au.

        Output:
        None
        """
        
        self.e = float(X[0])
        if type(X[1]) != pynbody.array.SimArray:
            self.a = SimArray(X[1],'au') #Assumes it's in au
        else:
            self.a = X[1]
        self.i = float(X[2])
        self.Omega = float(X[3])
        self.w = float(X[4])
        self.nu = float(X[5])

    # end function

    def reshapeData(self):
        """
        Ensure data, specifically r and v arrays, are of the proper shape for further manipulation.  
        Also ensure data has proper units
        This is useful because sometimes they come in as a list, a (1,3) numpy array or a (3,) numpy array.
        Much easier to clean up upon initialization then have many checks in later functions.
        """
        self.r = self.r.reshape((1, 3))
        self.v = self.v.reshape((1, 3))
        
        #If masses aren't SimArrays in units M_sol, make that so
        if type(self.m1) != pynbody.array.SimArray:
            self.m1 = SimArray(self.m1,'Msol')
        if type(self.m2) != pynbody.array.SimArray:
            self.m2 = SimArray(self.m2,'Msol')

    # end function

    def computeOrbElems(self):
        """
        Compute the Kepler orbital elements.

        Input: Self (must have r, v set and in proper units)

        Output: sets and returns e,a,i,...
        """
        # Compute orbital elements from binary center of mass frame Cartesian
        # coordinates
        assert (self.state != 'Kepler' or self.state != 'kepler'), "Already have orbital elements."        
        
        zeroR = SimArray([[0.0, 0.0, 0.0]],'cm')
        zeroV = SimArray([[0.0, 0.0, 0.0]],'cm s**-1')
        oe = AddBinary.calcOrbitalElements(
            self.r,
            zeroR,
            self.v,
            zeroV,
            self.m1,
            self.m2)

        # Set orbital elements, return them as well
        self.assignOrbElems(oe)
        return oe

    # end function

    def computeCartesian(self):
        """
        Compute the Cartesian position and velocity in the reduced mass frame.
        """
        assert (self.state == "Kepler" or self.state ==
                "kepler"), "Already have cartesian coordinates."
        M = AddBinary.trueToMean(self.nu, self.e)
        self.r, self.v = AddBinary.keplerToCartesian(
            self.a, self.e, self.i, self.Omega, self.w, M, self.m1, self.m2)

    # end function

    def generateICs(self):
        """
        From Kepler orbital elements, compute the position, velocities for two stars in ChaNGa-friendly units.
        Called "generateICs" since I'll use this mainly to...generate initial conditions
        """
        x1, x2, v1, v2 = AddBinary.reduceToPhysical(
            self.r, self.v, self.m1, self.m2)
        x1 = np.asarray(x1).reshape((1, 3))
        x2 = np.asarray(x2).reshape((1, 3))
        v1 = np.asarray(v1).reshape((1, 3))
        v2 = np.asarray(v2).reshape((1, 3))
        return x1, x2, v1, v2

    # end function
