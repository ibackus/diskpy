"""
David Fleming


Note:  All functions in the module only interact with the binary system itself.  For functions that work with
the disk or binary + disk, consult binaryUtils.py

Has all functions and utilities required to initialize/analyze binary star system

Initial input:
-Snapshot with central star of mass M (in Msol) located at (0,0,0)
in tipsy format (works with pynbody simArray data structure)
-Period (days), eccentricity of binary system

What This does:
Converts central star into 2 stars 1,2 under following conditions:
    m1 + m2 = M (divide mass into 2 stars)
    Center of Mass of stars remains at (0,0,0)

Used to initialize velocites of binary system given Keplerian relations
Also has additional functions to initialize and analyze binary system

#!!! Note: v_unit_vel will always be 29.785598165 km/s when m_unit = Msol and r_unit = 1 AU in kpc!!!
"""
#Constants and includes
import numpy as np
import math
from scipy import optimize
import sys
sys.path.append('/astro/users/dflemin3/Desktop/ICgen')
import isaac
import pynbody
SimArray = pynbody.array.SimArray

# Units/Constants
Msol = 1.98855e33  # g/Solar mass
BigG = 6.67259e-8  # in cgs
G = SimArray(BigG,'cm**3 g**-1 s**-2')
YEARSEC = 3.15569e7  # seconds per year
DAYSEC = 86400  # seconds per day
AUCM = 1.49597571e13  # cm/au
RAD2DEG = 180.0 / np.pi
SMALL = 1.0e-10  # less than this is zero enough

# ICgen-Specific constants (Shouldn't need to use these!)
VEL_UNIT = 29.785598165  # 29.785598165 km/s
POS_UNIT = 4.84813680873e-9  # in kpc == 1 au

# Function prototypes

# Binary Star Initilization Functions


def pToA(period=1, M=1):
    """
    Converts period (in days) into semimajor axis (in au) given Kepler law

    Parameters
    ----------
    Period : float
        [days]
    M : float
        COM mass [Msol]

    Returns
    -------
    Semimajor axis a (au)
    """
    conv = (DAYSEC * DAYSEC * Msol) / (AUCM * AUCM * AUCM)
    a = conv * period * period * BigG * M / (4.0 * np.pi * np.pi)
    return pow(a, 1.0 / 3.0)

# end function


def aToP(a=1, M=1):
    """
    Given a semimajor axis (au), convert into period (in days) given Kepler law

    Parameters
    ----------
    a : float
        semimajor axis [au]
    M : float
        mass of system [Msol]

    Returns
    -------
    Period : float
        [days]
    """
    conv = (AUCM * AUCM * AUCM) / (DAYSEC * DAYSEC * Msol)
    P = 4.0 * conv * np.pi * np.pi * a * a * a / (BigG * M)
    return np.sqrt(P)

# end function

##########################################################################
#                                                                                                        #
#	Functions for naively initializing binary stars at perihelion with arg peri = 0, LoAN = 0, inc = 0   #
#                                                                                                        #
##########################################################################


def calcPositions(M=1, a=1, e=0, p=0.5):
    """
    Given total mass of system M, semimajor axis a, and percentage of total mass contained in primary star p,
    calculate positions of binary components keep COM at origin
    (ex: If 1Msol system and MPri = Msec = 0.5 -> M =1 -> p = 0.5)
    Assume stars start a perihelion

    Parameters
    ----------
    M: float
        Total mass of system (Msol)
    a: float
        Semimajor axis (au)
    e: float
        eccentricity
    p: float
        % of total mass contained in primary (m1)

    Returns
    -------
    x1, x2: float
        Semimajor axes of binary stars assuming that they start at perihelion.
    """

    # Compute masses
    m1 = p * M
    m2 = (M - m1)

    # Calculate each star's semi-major axis a1,a2
    a1 = (m2 / M) * a
    a2 = (m1 / M) * a

    # Assume both stars star at perihelion (true anomaly = 0)
    x1 = a1 * (1 - e * e) / (1 + e)
    x2 = -1 * a2 * (1 - e * e) / (1 + e)

    return x1, x2

# end function


def calcV(m1=0.5, m2=0.5, a=1, e=0):
    """
    Given total mass M, postions of stars at perihelion x1, x2, and eccentricity e, calculate the velocities of the stars
    assuming that they are located at the perihelion and rotate in same direction as disk (CCW)

    Parameters
    ----------
    m1, m2: floats
        masses of primary and secondary (Msol)
    x1, x2: floats
        are semimajor axes of primary and secondary (au)
    e: float
        eccentricity

    Returns
    -------
    v1, v2: floats
        velocities of m1, m2 in km/s oriented for CCW rotation (in xy plane)
    """
    # Correct units and conversion factors, sqrt of positive numbers
    M = m1 + m2
    econv = (Msol) / (AUCM * 100 * 100 * 1000 * 1000)

    # Elliptical Orbit: initialize CCW velocity given orbit starts at
    # perihelion
    eps = (1 + e) / (1 - e)
    mu = (m1 * m2) / M
    vp = math.sqrt(econv * (BigG * M * eps) / (a))
    v1 = (mu / m1) * vp  # positive y direction for primary
    v2 = (-mu / m2) * vp  # negative y direction for secondary

    return v1, v2

# end function


def calcCriticalRadius(a=1, e=0, m1=0.5, m2=0.5):
    """
    Calculates the approximate bounds for where we would expect a planet to form
    given the conditions of a circumbinary disk around a short-period binary system.
    Calculates based on best fit from Holman&Wiegert+1999 (outer/P-type region)
    Assumes m1 ~ m2 and NOT m1 >> m2 or m2 >> m1

    Parameters
    ----------
    a: float
        Semimajor axis a of the binary system (au)
    e: float
        Eccentricity e of binary system
    m1, m2: floats
        Masses of binary components m1, m2 (Msol)

    Returns
    -------
    ac, pmac: floats
        Lower bounds for circumbinary planet's distace from binary COM and error terms of the following form:
        ac, pmac (error bounds are symmetric) both in au
    """

    # Compute mass ratio
    mu = m2 / (m1 + m2)

    # Compute critical following Holman&Wiegert+1999 with symmetric error for
    # the outer region
    ac = 1.60 + (5.10 * e) + (-2.22 * e * e) + (4.12 * mu) + \
        (-4.27 * e * mu) + (-5.09 * mu * mu) + (4.61 * e * e * mu * mu)
    ac *= a
    pmac = 0.04 + (0.05 * e) - (0.11 * e * e) + (0.09 * mu) - \
        (.17 * e * mu) - (0.11 * mu * mu) + (0.36 * e * e * mu * mu)
    pmac *= a

    return float(ac), float(pmac)

# end function

##########################################################################
#                                                                                                        #
#	Functions for computing Keplerian Orbital elements from Cartesian coordinates.  Most useful for      #
#	computing for binaries, but general enough to compute for gas particles or whatever.                 #
#                                                                                                        #
##########################################################################


def calcOrbitalElements(x1, x2, v1, v2, m1, m2):
    """
    Given as pynbody SimArrays the cental mass(es), the coodinate(s) and velocity(ies) of a CCW orbiting object,
    return the following orbital elements: eccentricity, semimajor axis, inclination, longitude of ascending node,
    argument of periapsis, and true anomaly.  This function is designed to work for a binary star system but is
    general enough to also work for a ~massless gas particle orbiting around a central quasi-Keplerian mass.

    Parameters
    ----------
    All input assumed to be in simulation units and are converted to desired units internally.
    x1,x2: position arrays in AU (x2 = 0 for gas particle case)
    v1,v2: velocity arrays in km/s (v2 = 0 for gas particle case)
    m1,m2: Central masses in Msol

    Returns
    -------
    e: float 
        Eccentricity (unitless)
    a: float
        Semimajor Axis in Au
    i: float
        Inclination in degrees
    Omega: float
        Longitude of Ascending node in degrees
    w: float
        Argument of Periapsis in degrees
    nu: float 
        True Anomaly in degrees
    """
    # Compute elements.  All unit conversion/processing done in sub functions
    e = calcEcc(x1, x2, v1, v2, m1, m2)
    a = calcSemi(x1, x2, v1, v2, m1, m2)
    i = calcInc(x1, x2, v1, v2)
    Omega = calcLongOfAscNode(x1, x2, v1, v2)
    w = calcArgPeri(x1, x2, v1, v2, m1, m2)
    nu = calcTrueAnomaly(x1, x2, v1, v2, m1, m2)

    return e, a, i, Omega, w, nu

# end function


def calcEcc(x1, x2, v1, v2, m1, m2, flag=True):
    """
    Given as pynbody arrays the masses of the binary system, arrays of the positions and velocities, compute
    its orbital eccentricity.

    Calculates e using following: e = sqrt(1 + (2*e*h^2)/(mu^2)
    for h = r x v, mu = G*(m1+m2), e = (v^2)/2 - (mu/|r|)

    Parameters
    ----------
    All inputs expected to be pynbody simArrays!!!
    
    x1, x2: SimArrays
        Position arrays of primary and secondary x1, x2 (in AU)
    v1, v2: SimArrays
        Velocity arrays of primary and secondary v1, v2 (in km/s)
    m1, m2: SimArrays
        Masses of primary and secondary m1, m2 (in Msol)
    Flag: bool 
        Whether or not to internally convert to cgs units

    Returns
    -------
    e: float
        Scalar eccentricity of binary system.
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        #m1 = np.asarray(isaac.strip_units(m1)) * Msol
        #m2 = np.asarray(isaac.strip_units(m2)) * Msol

    length, ax = computeLenAx(x1)

    # Relative position vector in cgs
    r = (x1 - x2)
    magR = SimArray(np.linalg.norm(r, axis=ax),'cm')
    
    # Compute standard gravitational parameter in cgs
    mu = (G * (m1 + m2)).in_units('cm**3 s**-2')

    # Compute relative velocity vector in cgs with appropriate scale
    v = (v1 - v2)
    magV = SimArray(np.linalg.norm(v, axis=ax),'cm s**-1')
    
    # Compute specific orbital energy
    eps = ((magV * magV / 2.0) - (mu / magR)).in_units('cm**2 s**-2')

    # Compute specific angular momentum vector
    h = SimArray(np.cross(r, v, axis=ax),'cm**2 s**-1')
    magH = SimArray(np.linalg.norm(h, axis=ax),'cm**2 s**-1')

    # Compute, return eccentricity
    return np.sqrt(1.0 + ((2.0 * eps * magH * magH) / (mu * mu)))

# end function


def calcSemi(x1, x2, v1, v2, m1, m2, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the binary's semimajor axis.

    Calculates a using the following: a = -mu/(2*e)
    where mu = G*(m1+m2) and e = (v^2)/2 - (mu/|r|)

    Parameters
    ----------
    (as pynbody SimArrays!)
    
    x1,x2: SimArrays
        Primary and secondary position arrays [AU]
    v1,v2: SimArrays
        Primary and secondary velocity arrays [km/s]
    m1,m2: SimArrays
        Primary and secondary masses [Msol]
    Flag: bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    a: float
        semimajor axis of binary orbit in AU
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        # Remove units since input is pynbody SimArray
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        #m1 = np.asarray(isaac.strip_units(m1)) * Msol
        #m2 = np.asarray(isaac.strip_units(m2)) * Msol

    length, ax = computeLenAx(x1)

    # Relative position vector in cgs
    r = (x1 - x2)
    magR = SimArray(np.linalg.norm(r, axis=ax),'cm')

    # Compute standard gravitational parameter in cgs
    mu = (G * (m1 + m2))

    # Compute relative velocity vector in cgs with appropriate scale
    v = (v1 - v2)
    magV = SimArray(np.linalg.norm(v, axis=ax),'cm s**-1')

    # Compute specific orbital energy
    eps = (magV * magV / 2.0) - (mu / magR)

    # Compute, return semimajor axis in AU (convert from cgs->AU)
    return (-mu / (2.0 * eps)).in_units('au') #/ (AUCM)

# end function


def calcInc(x1=1, x2=0, v1=1, v2=0, flag=True):
    """
    Given pynbody arrays for positions and velocities of primary and secondaries bodies and masses in
    a binary orbit, calculate's the orbit's inclination. (Given: Orbit starts in xy plane)

    i = arccos(h_z/|h|)

    Parameters
    ----------
    as pynbody SimArrays [preferred units]
    
    x1,x2: SimArrays
        Primary and secondary position arrays [AU]
    v1, v2: SimArrays
        Primary and secondary velocity arrays [km/s]
    Flag: bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    i: float
        Inclination in degrees
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        # Strip units from all inputs
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT

    # Compute length of array we're dealing with
    length, ax = computeLenAx(x1)

    # Relative position vector in cgs
    r = (x1 - x2)

    # Compute relative velocity vector in cgs with appropriate scale
    v = (v1 - v2)

    # Compute specific angular momentum vector
    h = np.cross(r, v, axis=ax)
    magH = SimArray(np.linalg.norm(h, axis=ax),'cm**2 s**-1')

    if(length > 1):
        h_z = h[:, 2]
    else:
        h_z = h[0, 2]

    # Orbit is CCW (h_z < 0) so take fabs to have i >= 0
    h_z = np.fabs(h_z)

    # Compute i, convert to degrees
    i = np.arccos(h_z / magH)

    return i * RAD2DEG  # return in degrees

# end function


def calcLongOfAscNode(x1=1, x2=0, v1=1, v2=0, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the binary's longitude of the ascending node Omega.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    Calculates Omega using the following: Omega = arccos(n_x/|n|) n_y > 0
    Omega = 2*pi - arccos(n_x/|n|) n_y < 0
    where n = (0,0,1) x h for h = r x v

    Parameters
    ----------    
    Assumed as pynbody SimArrays in simulation units (AU, scaled velocity, etc)
    
    x1,x2: SimArrays
        Primary and secondary position arrays [AU]
    v1,v2: SimArrays
        Primary and secondary velocity arrays [km/s]

    Returns
    -------
    Omega: float
        longitude of the ascending node in degrees
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        # Strip units from all inputs and convert to cgs
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT

    # Define unit vectors pointing along z, x and y axes respectively
    # Also ensure function can handle any number of values
    length, ax = computeLenAx(x1)

    k = np.zeros((length, 3))
    i = np.zeros((length, 3))
    j = np.zeros((length, 3))

    if(length > 1):
        j[:, 1] = 1
        i[:, 0] = 1
        k[:, 2] = 1
    else:
        j[0, 1] = 1
        i[0, 0] = 1
        k[0, 2] = 1

    # Relative position vector in cgs
    r = (x1 - x2)

    # Compute relative velocity vector in cgs with appropriate scale
    v = (v1 - v2)

    # Compute specific angular momentum vector
    h = np.cross(r, v, axis=ax)

    # Compute vector pointing to ascending node
    n = np.cross(k, h, axis=ax)
    magN = SimArray(np.linalg.norm(n, axis=ax),'cm**2 s**-1')

    # Ensure no divide by zero errors?
    #magN[magN < SMALL] = 1.0

    # Compute LoAN
    inc = calcInc(x1, x2, v1, v2)/RAD2DEG
    Omega = np.arccos(dotProduct(i, n) / magN)

    # If inclination is ~0, define LoAN as 0
    Omega[inc < SMALL] = 0.0

    # Fix phase due to arccos return range
    Omega[dotProduct(n, j) < 0] = 2.0 * np.pi - Omega[dotProduct(n, j) < 0]

    # Convert to degrees, return
    return Omega * RAD2DEG

# end function


def calcEccVector(x1=1, x2=0, v1=1, v2=0, m1=1, m2=1, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the eccentricity vector in the reduced two body system.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    Parameters
    ----------
    x1,x2: SimArrays
        Primary and secondary position arrays [length]
    v1, v2: SimArrays
        Primary and secondary velocity arrays [velocity]
    m1,m2: SimArrays
        Primary and secondary masses [mass]
    Flag: bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    Ecc: array
        Eccentricity vector in cgs
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        # Remove units in case input is pynbody SimArray
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        #m1 = np.asarray(isaac.strip_units(m1)) * Msol
        #m2 = np.asarray(isaac.strip_units(m2)) * Msol

    # Determine length of arrays
    length, ax = computeLenAx(x1)

    # Relative position vector in cgs
    r = (x1 - x2)
    magR = SimArray(np.linalg.norm(r, axis=ax).reshape(len(r),1),'cm')

    # Compute standard gravitational parameter in cgs
    mu = (G * (m1 + m2)).reshape(len(r),1)

    # Compute relative velocity vector in cgs with appropriate scale
    v = (v1 - v2)

    # Compute specific angular momentum vector
    h = SimArray(np.cross(r, v, axis=ax),'cm**2 s**-1')

    # Compute, return eccentricity vector
    return (SimArray(np.cross(v, h, axis=ax),'cm**3 s**-2') / mu) - (r / magR)

# end function


def calcArgPeri(x1=1, x2=0, v1=1, v2=0, m1=1, m2=1, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the eccentricity vector in the reduced two body system.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    Parameters
    ----------
    Assumed as pynbody SimArrays [preferred units]
    
    x1,x2 : SimArrays
        Primary and secondary position arrays [AU]
    v1,v2 : SimArrays
        Primary and secondary velocity arrays [km/s]
    m1, m2 : SimArrays
        Primary and secondary masses [Msol]
    Flag : bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    w: float
        Argument of pericenter in degrees
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        # Remove units since input is pynbody SimArray
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        #m1 = np.asarray(isaac.strip_units(m1)) * Msol
        #m2 = np.asarray(isaac.strip_units(m2)) * Msol

    # Compute eccentricity vector
    e = calcEccVector(x1, x2, v1, v2, m1, m2, flag=False)
    magE = SimArray(np.linalg.norm(e, axis=1),'1')

    length, ax = computeLenAx(x1)

    # Define unit vector pointing along z axis
    k = np.zeros((length, 3))
    if length > 1:
        k[:, 2] = 1.0
    else:
        k[0, 2] = 1.0

    # Define specific angular momentum vector
    # Relative position vector in cgs
    r = (x1 - x2)

    # Compute relative velocity vector in cgs with appropriate scale
    v = (v1 - v2)

    # Compute specific angular momentum vector
    h = SimArray(np.cross(r, v, axis=ax),'cm**2 s**-1')

    # Compute vector pointing to ascending node
    n = np.cross(k, h, axis=ax)
    magN = SimArray(np.linalg.norm(n, axis=ax),'cm**3 s**-1')

    # Ensure no divide by zero errors?
    magN[magN < SMALL] = 1.0

    # Compute argument of periapsis
    inc = calcInc(x1, x2, v1, v2)/RAD2DEG
    arg = dotProduct(n, e) / (magN * magE)
    
    # Bounds check arg    
    w = np.arccos(arg)
    #w[arg < -1.0] = np.pi
    #w[arg > 1.0] = 0.0
    w[dotProduct(e, k) < 0] = 2.0 * np.pi - w[dotProduct(e, k) < 0.0]
    w[inc < SMALL] = 0.0  # For orbit in a plane

    #Return in degrees
    return w * RAD2DEG

# end function


def calcTrueAnomaly(x1=1, x2=0, v1=1, v2=0, m1=1, m2=1, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the true anomaly in the reduced two body system.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    Parameters
    ----------
    Assumed as pynbody SimArrays [preferred units]
    
    x1,x2: SimArrays
        Primary and secondary position arrays [AU]
    v1,v2: SimArrays
        Primary and secondary velocity arrays [km/s]
    m1, m2: SimArrays
        Primary and secondary masses [Msol]
    Flag: bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    nu: float
        True anomaly in degrees
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        # Remove units since input is pynbody SimArray
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        #m1 = np.asarray(isaac.strip_units(m1)) * Msol
        #m2 = np.asarray(isaac.strip_units(m2)) * Msol

    # Compute length, correct axis
    length, ax = computeLenAx(x1)

    # Compute eccentricity vector
    e = calcEccVector(x1, x2, v1, v2, m1, m2, flag=False)
    magE = SimArray(np.linalg.norm(e, axis=ax),'1')

    # Compute radius vector
    r = (x1 - x2)
    v = (v1 - v2)
    magR = SimArray(np.linalg.norm(r, axis=ax),'cm')

    # Compute true anomaly making sure I can handle single numbers or arrays
    nu = np.arccos(dotProduct(e, r) / (magE * magR))
    if isinstance(nu, np.float64):
        if dotProduct(r, v) < 0.0:
            nu = 2.0 * np.pi - nu
    else:
        nu[dotProduct(r, v) < 0.0] = 2.0 * np.pi - nu[dotProduct(r, v) < 0.0]

    # Convert to degrees, return
    return nu * RAD2DEG

# end function


def calcEccentricAnomaly(x1=1, x2=0, v1=1, v2=0, m1=1, m2=1, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the eccentric anomaly in the reduced two body system.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    Parameters
    ----------
    Assumed as pynbody SimArrays [preferred units]
    
    x1,x2: SimArrays
        Primary and secondary position arrays [AU]
    v1,v2: SimArrays
        Primary and secondary velocity arrays [km/s]
    m1, m2: SimArrays
        Primary and secondary masses [Msol]
    Flag: bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    E: float
        Eccentric anomaly in degrees
    """
    #This if/else is stupid and I should change it eventually--dflemin3
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        # Remove units since input is pynbody SimArray
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        #m1 = np.asarray(isaac.strip_units(m1)) * Msol
        #m2 = np.asarray(isaac.strip_units(m2)) * Msol
        e = calcEcc(x1, x2, v1, v2, m1, m2, flag=False)
        nu = calcTrueAnomaly(x1, x2, v1, v2, m1, m2, flag=False)
    else:
        e = calcEcc(x1, x2, v1, v2, m1, m2, flag=False)
        nu = calcTrueAnomaly(x1, x2, v1, v2, m1, m2, flag=False)

    # Calc E
    nu = nu * (np.pi / 180.0)  # convert to radians for numpy functions
    E = np.arccos((e + np.cos(nu)) / (1.0 + e * np.cos(nu)))

    # Make sure this can handle single numbers or arrays
    if isinstance(E, np.float64):
        if nu > np.pi and nu < 2.0 * np.pi:
            E = 2.0 * np.pi - E
    else:
        E[np.logical_and(nu > np.pi, nu < 2.0 * np.pi)] = 2.0 * np.pi - E

    # Return E in degrees
    return E * RAD2DEG


def calcMeanAnomaly(x1=1, x2=0, v1=1, v2=0, m1=1, m2=1, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the Mean anomaly in the reduced two body system.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    Parameters
    ----------
    Assumed as pynbody SimArrays [preferred units]
    
    x1,x2: SimArrays
        Primary and secondary position arrays [AU]
    v1,v2: SimArrays
        Primary and secondary velocity arrays [km/s]
    m1, m2: SimArrays
        Primary and secondary masses [Msol]
    Flag: bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    M: float
        Mean anomaly in degrees
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        # Remove units since input is pynbody SimArray
        #x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        #x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        #v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        #v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        #m1 = np.asarray(isaac.strip_units(m1)) * Msol
        #m2 = np.asarray(isaac.strip_units(m2)) * Msol
        e = calcEcc(x1, x2, v1, v2, m1, m2, flag=False)
        E = calcEccentricAnomaly(x1, x2, v1, v2, m1, m2, flag=False)
    else:
        e = calcEcc(x1, x2, v1, v2, m1, m2, flag=False)
        E = calcEccentricAnomaly(x1, x2, v1, v2, m1, m2, flag=False)

    # Calculate Mean Anomaly
    E = E * (np.pi / 180.0)  # Conver E to radians for numpy
    M = E - e * np.sin(E)

    # Return M in degrees
    return (M * RAD2DEG)

# end function


def trueToMean(nu, e):
    """
    Given the true anomaly nu in degrees and the eccentricity e, compute the mean anomaly M in degrees.

    Parameters
    ----------
    nu: float
        True anomaly (degrees)
    e: float
        eccentricity

    Returns
    -------
    M: float
        mean anomaly (degrees)
    """
    # Compute eccentric anomaly E
    nu = nu * (np.pi / 180.0)  # convert to radians for numpy functions
    E = np.arccos((e + np.cos(nu)) / (1.0 + e * np.cos(nu)))

    # Make sure this can handle single numbers or arrays
    if isinstance(E, np.float64):
        if nu > np.pi and nu < 2.0 * np.pi:
            E = 2.0 * np.pi - E
    else:
        E[np.logical_and(nu > np.pi, nu < 2.0 * np.pi)] = 2.0 * np.pi - E

    # Compute, return M
    M = E - e * np.sin(E)
    return (M * RAD2DEG)

#end function

def calcMeanMotion(x1, x2, v1, v2, m1, m2, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the Mean motion in the reduced two body system.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    Parameters
    ----------
    Assumed as pynbody SimArrays [preferred units]
    
    x1,x2: SimArrays
        Primary and secondary position arrays [AU]
    v1,v2: SimArrays
        Primary and secondary velocity arrays [km/s]
    m1, m2: SimArrays
        Primary and secondary masses [Msol]
    Flag: bool
        Whether or not to internally convert to cgs units

    Returns
    -------
    n: SimArray
        Mean motion in 1/s
    """
    if flag:
        #Ensure units are in cgs
        x1 = x1.in_units('cm')
        x2 = x2.in_units('cm')
        v1 = v1.in_units('cm s**-1')
        v2 = v2.in_units('cm s**-1')
        m1 = m1.in_units('g')
        m2 = m2.in_units('g')
        
    #Put binary in reduced mass frame, compute, return
    a = SimArray(calcSemi(x1,x2,v1,v2,m1,m2,flag=False),'au').in_units('cm')
    mu = G * (m1 + m2)
    return np.sqrt(mu/np.power(a,3))
    
#end function

##########################################################################
#                                                                                               #
#	Functions for computing Cartesian coordinates from Keplerian orbital elements.              #
#                                                                                               #
##########################################################################


def keplerToCartesian(
        a,
        e,
        i,
        Omega,
        w,
        M,
        m1,
        m2,
        angleFlag=True,
        scaleFlag=True):
    """
    Given the Keplerian orbital elements, compute the cartesian coordinates of the object orbiting
    in the reduced mass frame.  Note: Requires all angles in degrees unless noted.

    Note: A little redudant that I compute M when I typically already know the true anomaly nu, but most
    other schemes know M initially instead of nu so I'll keep it for compatibility's sake.

    Parameters
    ----------
    a: float
        Semimajor axis (AU)
    e: float
        Eccentricity
    i: float
        inclination (degrees)
    Omega: float
        Longitude of Ascending Node (degrees)
    w: float
        Argument of Pericenter (degrees)
    M: float
        Mean Anomaly (degrees)
    m1, m2: float
        Masses of central object(s) (Msol)
    angleFlag, scaleFlag: bools
        Tells code to convert degrees->rad and/or scale velocity to sim units, respectively

    Returns
    -------
    x: numpy array
        Position array of the object in reduced mass frame (AU)
    v: numpy array
        velocity array (km/s) with VEL_UNIT scaling factor optional
    """
    # Compute length to allow for multiple objects
    length, ax = computeLenAx(np.asarray(e))

    # Convert everything to radians!
    if angleFlag:
        i /= RAD2DEG
        Omega /= RAD2DEG
        w /= RAD2DEG
        M /= RAD2DEG

    # Convert M->E by solving M = E-esinE for E
    # Rearrange to E - esinE - M = 0 to solve for root
    def F(E, m, ecc):
        return E - ecc * np.sin(E) - m

    # Compute eccentric anomaly in radians
    E = optimize.newton(F, M, args=(M, e))

    # Compute unit vectors P, Q along axes of PQW frame
    # These vectors will transform to proper barycentric frame given i, Omega,
    # w orbital params
    P = np.zeros((length, 3))
    Q = np.zeros((length, 3))
    if(length > 1):
        P[:, 0] = np.cos(w[:]) * np.cos(Omega[:]) - \
            np.sin(w[:]) * np.cos(i[:]) * np.sin(Omega[:])
        P[:, 1] = np.cos(w[:]) * np.sin(Omega[:]) + \
            np.sin(w[:]) * np.cos(i[:]) * np.cos(Omega[:])
        P[:, 2] = np.sin(w[:]) * np.sin(i[:])
        Q[:, 0] = -np.sin(w[:]) * np.cos(Omega[:]) - \
            np.cos(w[:]) * np.cos(i[:]) * np.sin(Omega[:])
        Q[:, 1] = -np.sin(w[:]) * np.sin(Omega[:]) + \
            np.cos(w[:]) * np.cos(i[:]) * np.cos(Omega[:])
        Q[:, 2] = np.sin(i[:]) * np.cos(w[:])

    else:
        P[0, 0] = np.cos(w) * np.cos(Omega) - np.sin(w) * \
            np.cos(i) * np.sin(Omega)
        P[0, 1] = np.cos(w) * np.sin(Omega) + np.sin(w) * \
            np.cos(i) * np.cos(Omega)
        P[0, 2] = np.sin(w) * np.sin(i)
        Q[0, 0] = -np.sin(w) * np.cos(Omega) - np.cos(w) * \
            np.cos(i) * np.sin(Omega)
        Q[0, 1] = -np.sin(w) * np.sin(Omega) + np.cos(w) * \
            np.cos(i) * np.cos(Omega)
        Q[0, 2] = np.sin(i) * np.cos(w)

    # Compute Standard Gravitational Parameter in cgs assuming masses in Msol
    mu = BigG * (m1 + m2) * Msol

    # Compute radius vector in AU...assumes a already in AU!
    r = a * (np.cos(E) - e) * P + a * np.sqrt(1.0 - e * e) * np.sin(E) * Q

    # Compute velocity vector
    tmp = np.sqrt(mu) / (np.power(a, 1.5) * (1.0 - e * np.cos(E)))
    v = -a * e * np.sin(E) * tmp * P
    v += a * np.sqrt(1.0 - e * e) * np.cos(E) * tmp * Q

    # Convert v to km/s and in sim units if flag says yes
    conv = 1.0 / (np.sqrt(AUCM) * 100.0 * 1000.0)

    if scaleFlag:
        v *= conv / VEL_UNIT
    else:
        v *= conv

    return r, v

# end function


def reduceToPhysical(r, v, m1, m2):
    """
    Function converts from reduced mass coordinates to physical, origin-centered coords
    Works for arbitrary units, number of objects with coordinates as long as m1+m2 = central mass.

    Parameters
    ----------
    r: array
        radius vector
    v: array
        velocity vector
    m1, m2: floats
        mass of central object

    Returns
    --------
    x1, x2: numpy arrays
        position vectors of 2 mutually orbiting objects
    v1, v2:numpy arrays
        velocity vectors of 2 mutually orbiting objects
    """
    # Compute reduced mass
    M = m1 + m2
    mu = (m1 * m2) / M

    # Compute positions of the form:  x2----origin--x1 where m1 >= m2
    x1 = (mu / m1) * r
    x2 = -(mu / m2) * r

    # Compute velocities for CCW orbit such that m2 v_y < 0, m1 v_y > 0 at
    # pericenter
    v1 = (mu / m1) * v
    v2 = -(mu / m2) * v

    return x1, x2, v1, v2

# end function


def initializeBinary(
        a,
        e,
        i,
        Omega,
        w,
        M,
        m1,
        m2,
        angleFlag=True,
        scaleFlag=True):
    """
    Given the initial Kepler orbital parameters, compute the Cartesian positions and velocities
    for 2 mutually orbiting stars (binaries!).

    Parameters
    ----------
    Keplerian orbital elements: floats
        in Au, degrees (see above for more details)
    m1,m2: floats
        Masses of central objects (Msol)
    angleFlag: bool
        whether or not to convert from degrees->radians (True = convert)
    scaleFlag: bool
        whether or not to put v in sim units (True = do it/default option)

    Returns
    -------
    x1, x2: array
        Positions of 2 objects about origin (Au)
    v1, v2: array
        Velocities of 2 objects about origin (km/s in Sim Units)

    """
    r, v = keplerToCartesian(
        a, e, i, Omega, w, M, m1, m2, angleFlag, scaleFlag)

    return reduceToPhysical(r, v, m1, m2)

# end function

##########################################################################
#                                                                                               #
#	Misc functions I found useful to impliment for the above even more useful functions.       #
#                                                                                               #
##########################################################################


def accretionEDot(Binary, M_dot, dt):
    """
    Given a Binary object and an accretion rate (assumed to be in M_sol/yr), compute the theoretical rate of change of the
    binary's eccentricity in 1/second.
    Assumptions:
    -radius, velocity of binary nearly constant over accretion (found to more or less apply via empirical measurements)

    Parameters
    ----------
    Binary: binary object class
    Mdot: float
        accretion rate in M_sol/yr
    dt: float
        time interval

    Returns
    -------    
    de/dt: float
        change in eccentricity in 1/second
    """
    # Convert relevant quantities into cgs
    a = Binary.a * AUCM
    e = Binary.e
    M = (Binary.m1 + Binary.m2) * Msol
    mu = BigG * M
    Mdot = M_dot*Msol
    #dm = Mdot*dt
    r = a * (1.0 + e)
    eps = -mu / (2.0 * a)
    v = np.sqrt((BigG * M / a) * ((1.0 - e) / (1.0 + e)))
    h = r * v
    #Hdot = BigG*dm*dt/r
    Hdot = Mdot * r * v / M

    edot = 2.0 * eps * h * Hdot / (mu * mu)
    edot /= np.sqrt(1.0 + (2.0 * eps * h * h) / (mu * mu))

    return edot

# end function

def binaryPrecession(s,r_in,r_out):
    """
    Computes the period of the binary precession due to an axisymmetric disk (!!!)
    given by Rafikov 2013.

    Parameters
    ----------
    s : Tipsy snapshot
    r_in, r_out : float
        inner and outer radii of the circumbinary disk [AU]
        
    Returns
    -------
    T : SimArray
        Period of binary argument of periastron precession in yr
    """
    x1 = s.stars[0]['pos'].in_units('cm')
    x2 = s.stars[1]['pos'].in_units('cm')
    v1 = s.stars[0]['vel'].in_units('cm s**-1')
    v2 = s.stars[1]['vel'].in_units('cm s**-1')
    m1 = s.stars[0]['mass'].in_units('g')
    m2 = s.stars[1]['mass'].in_units('g')    
    
    # Define required parameters in cgs
    M_bin = m1 + m2
    M_disk = np.sum(s.gas['mass']).in_units('g')    
    a = SimArray(calcSemi(x1,x2,v1,v2,m1,m2,flag=False),'au').in_units('cm') #semimajor axis in cm
    n = calcMeanMotion(x1, x2, v1, v2, m1, m2, flag=False) #mean motion in 1/s    
    r_in = SimArray(r_in,'au').in_units('cm')
    r_out = SimArray(r_out,'au').in_units('cm')    
    
    T = 8.0*np.pi*(M_bin/M_disk)*(np.power(r_out,0.5)*np.power(r_in,2.5)/(np.power(a,3)*0.5*n))
    return isaac.strip_units(T)/YEARSEC
    
#end function


def calcCircularFrequency(x1, x2, v1, v2, m1, m2, flag=True):
    """
    Given pynbody arrays for positions and velocity of primary and secondary bodies
    and masses, calculates the circular frequency in the reduced two body system.
    Usage note: Intended for binary system, but pass x2 = v2 = 0 to use with any
    location in the disk.

    omega = (L_z)/(R^2) which assumes spherical symmetry (ok assumption here)
    L = sqrt(G*M*a*(1-e^2) for ~Keplerian

    Parameters
    ----------
    Assumed as pynbody SimArrays in simulation units (AU, scaled velocity, etc)
    
    x1,x2: arrays
        Primary and secondary position arrays (in AU)
    v1, v2: arrays
        Primary and secondary velocity arrays (km/s)
    m1, m2: floats
        Primary and secondary masses (Msol)
    Flag: Whether or not to internally convert to cgs units

    Returns
    -------
    omega: float
        Circular frequency in 1/days
    """
    e = calcEcc(x1, x2, v1, v2, m1, m2, flag=True)
    a = calcSemi(x1, x2, v1, v2, m1, m2, flag=True) * AUCM
    
    if flag:
        # Remove units since input is pynbody SimArray
        x1 = np.asarray(isaac.strip_units(x1)) * AUCM
        x2 = np.asarray(isaac.strip_units(x2)) * AUCM
        v1 = np.asarray(isaac.strip_units(v1)) * 1000 * 100 * VEL_UNIT
        v2 = np.asarray(isaac.strip_units(v2)) * 1000 * 100 * VEL_UNIT
        m1 = np.asarray(isaac.strip_units(m1)) * Msol
        m2 = np.asarray(isaac.strip_units(m2)) * Msol

    length, ax = computeLenAx(x1)

    # Calculate angular momentum assuming all arrays are nx3
    r = x1 - x2
    rMag = np.sqrt(dotProduct(r, r))
    #e = ...    
    #a = calcSemi(x1, x2, v1, v2, m1, m2, flag=False) * AUCM
    L = np.sqrt(BigG * (m1 + m2) * a * (1 - e * e))

    # Convert from 1/s to 1/day, return
    return L * DAYSEC / (rMag * rMag)

# end function


def calcCOM(m1=0.5, m2=0.5, x1=1, x2=1):
    """
    Given pynbody arrays for the mass and position of the two binary stars,
    function calculates the CoM of the two particles in order to check that
    it's still roughly 0.

    Parameters
    ----------
    
    (as pynbody SimArrays)
    m1, m2: SimArrays
        Primary and secondary mass arrays (in Msol)
    x1, x2: SimArrays
        Primary and secondary position arrays (in AU)

    Returns
    -------
    center of mass: array
        Center of mass position vector (numpy array in AU)
    """

    # Strip units from inputs
    x1 = np.asarray(isaac.strip_units(x1))
    x2 = np.asarray(isaac.strip_units(x2))
    m1 = np.asarray(isaac.strip_units(m1))
    m2 = np.asarray(isaac.strip_units(m2))

    # Compute,return CoM
    return (1.0 / (m1 + m2)) * ((m1 * x1) + (m2 * x2))

# end function


def calcRocheLobe(q, a):
    """
    Given the mass ratio q = m1/m2 and the semimajor axis in AU, compute the radius of the Roche lobe
    around m1 given the Eggleton approximation http://en.wikipedia.org/wiki/Roche_lobe

    Parameters
    ----------
    q: float
        Binary mass ratio q = m1/m2
    a: float
        Binary semimajor axis a (in AU)

    Returns
    -------
    r: float
        Radius of Roche lobe around m1 (in AU)
    """
    num = 0.49 * np.power(q, 2. / 3.)
    denom = 0.6 * np.power(q, 2. / 3.) + np.log(1.0 + np.power(q, 1. / 3.))

    return (num / denom) * a

# end function


def dotProduct(a, b):
    """
    Given n x m numpy arrays, compute the dot product row-wise.

    Parameters
    ----------
    a,b: nxm numpy array

    Returns
    -------
    dot product: numpy array 
        n rows containing dot products
    """
    # Don't trust user, make sure a, b are numpy arrays of same shape
    a = np.asarray(a)
    b = np.asarray(b)
    assert(a.shape == b.shape), "Input arrays must be same shape."

    # Determine correct axis to perform operation on
    l, ax = computeLenAx(a)

    # Numpy multiplies element wise!
    return np.sum(a * b, axis=ax)

# end function


def computeLenAx(a):
    """
    Given an numpy array, check to see if it's length is 1 (aka it's a float, or int or whatever) or otherwise and return
    it's length and the axis over which calculations like normalization and cross product is to be taken.  This function
    is useful when you expect data arrays of the form x = [n,3] with weird combinations of lists/numpy arrays
    """
    if(a.shape == () or a.shape == (len(a),)):  # Just a float
        length = 1
        ax = 0
    elif (a.shape == (1, len(a))):
        length = 1
        ax = 0
    else:
        length = len(a)
        ax = 1

    return length, ax

# end function
