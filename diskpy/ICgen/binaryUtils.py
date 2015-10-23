"""
David Fleming
Utilities to process/interact with binary star system in ChaNGa sims

Note on inputs: Most (if not all?) functions are designed to be used with SimArrays as inputs.  
"""

# Imports
import numpy as np
import re
import AddBinary
import pynbody
from scipy import interpolate
from scipy.optimize import fsolve
SimArray = pynbody.array.SimArray

from diskpy.utils import strip_units


def changaFloatSearch(name,simUnits=False):
	"""
	Given the name of a file containing line dumps for ChaNGa and outputs numpy arrays containing changa dumps line-by-line.
	
	Genfromtxt will convert output into either an numpy array of strings or a numpy array of lists of strings.
	Determine what format is, use regex library to find all floats and output accordingly.

	Default usage is searching for linear momentum dumps of the form mg,vx,vy,vz for gas
	Assume appropriate flag was used to grep data into input file
	
	Input:
	name: Name of input file (something.txt)
	simUnits: whether or not to use simUnits (False -> convert to cgs)
			  Only for use for linear momentum values!

	Output:
	Numpy array containing all floats from changa output
	"""
	
	#Read in data
	data = np.genfromtxt(name,dtype="str")

	#Determine size of intermediate output numpy array
	rows = len(data)

	result = np.zeros(rows,dtype=list)

	#Loop over inputs, store all floats in numpy array to return
	#Case: Each element of data is a string
	if isinstance(data[0], basestring):
		for i in range(rows):
			tmpList = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",data[i])
			if tmpList:
				result[i] = np.asarray(tmpList).astype(float)
	#Case: Each element of data is a list of strings
	else:
		for i in range(rows):
			for item in data[i]:
				#Check if list is empty or contains floats.  If not empty, convert to numpy array of floats, store in return array
				tmpList = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",item)
				if tmpList:
					result[i] = np.asarray(tmpList).astype(float)

	#Ensure final result is some 2D numpy array so I can access it like final[i,j]
	cols = len(result[0])
	final = np.zeros(shape=(rows,cols),dtype=float)
	for i in range(rows):
		temp = result[i]
		for j in range(cols):
			final[i,j] = float(temp[j])

	return final

#end function

def linearMomentumEffects(x1, x2, v1, v2, m1, m2, accretion):
	"""
	Given initial binary system parameters and an array tracking the accretion events, calculate the effects of accretion
	on the semimajor axis and eccentricity of the binary system.

	Inputs: Assume all input arrays are in simulation units
	Masses of primary and secondary (m1, m2 in Msol)
	Position arrays of primary and secondary x1, x2 (in AU)
	Velocity arrays of primary and secondary v1, v2 (in km/s)
	Numpy array of accretion events of the form [m vx vy vz ...] for each accreted gas particle at time of accretion

	Output:
	Semimajor axis, eccentricity of binary system after accretion events
	"""
	#Extract masses and velocities of accreted gas particles from array of known format
	m_g = np.zeros(len(accretion))	
	v = np.zeros(shape=(len(accretion),3))

	for i in range(len(accretion)):
		m_g[i] = accretion[i,0]	
		for j in range(1,4):
			v[i,j-1] = accretion[i,j]

	#Strip units from all inputs, convert all into CGS
	r1 = np.asarray(strip_units(x1))*AddBinary.AUCM
	r2 = np.asarray(strip_units(x2))*AddBinary.AUCM
	v1 = np.asarray(strip_units(v1))*AddBinary.VEL_UNIT*100*1000
	v2 = np.asarray(strip_units(v2))*AddBinary.VEL_UNIT*100*1000
	m1 = np.asarray(strip_units(m1))*AddBinary.Msol
	m2 = np.asarray(strip_units(m2))*AddBinary.Msol
	m_g = m_g*AddBinary.Msol
	v = v * AddBinary.VEL_UNIT*100*1000

	#Compute relative binary system quantities
	vBin = v1 - v2
	rBin = r1 - r2
	mBin = m1 + m2

	#Loop over accretion events, apply conservation of linear momentum at each step
	for i in range(len(accretion)):
		vBin = (1.0/(mBin+m_g[i]))*(mBin*vBin + m_g[i]*v[i])
		mBin = mBin + m_g[i]	

	#Compute final semimajor axis, eccentricity
	
	#Compute r, v, standard gravitational parameter
	magR = np.linalg.norm(rBin)
	mu = AddBinary.BigG*(mBin)
	magV = np.linalg.norm(vBin)

	#Compute specific orbital energy, angular momentum
	eps = (magV*magV/2.0) - (mu/magR)
	h = np.cross(rBin,vBin)
	magH = np.linalg.norm(h)
	
	#Compute semimajor axis	in AU
	a = -mu/(2.0*eps)/(AddBinary.AUCM)

	#Compute eccentricity
	e = np.sqrt(1 + ((2*eps*magH*magH)/(mu*mu)))

	return a,e

#end function

def find_crit_radius(r,array,toFind,num=1000):
    """
    Given an array as a function of radius,	the array value to search for, the radial range to search over and the search
    resolution, find the corresponding radius.
    Assumed that array is calculated as a function of r.
	
    Parameters
    ----------
    r: array
        array of radial points
    array: array
        array of quantity of interest (could be surface density?) as a function of r
    toFind: float
        array value you're looking for
    num: int
        resolution of search

    Returns
    -------
    critial_radius: float
        radius at which array(critical_radius) = toFind (approximately)
    """
    #Require len(r) == len(array) for interpolation to work
    assert len(r) == len(array)

    #Estimate surface density as function of radius.  s=0 -> Interpolate through all data points assuming smooth curve
    array_f = interpolate.UnivariateSpline(r,array,s=0)
	
    #Compute radius array to search over
    radius = np.linspace(r.min(),r.max(),num)

    #Find, return critical radius
    return radius[np.fabs(array_f(radius)-toFind).argmin()]

#end function

def computeCOM(stars,gas,cutoff=None,starFlag=True,gasFlag=True):
    """
    Given pynbody star and gas arrays, compute the center of mass for the entire specified system.

    Parameters
    ----------    
    stars: s.stars pynbody object
    gas: s.gas pynbody object
    cutoff: float
        radius at which you only consider objects interior to it
    starFlag: bool
        whether or not to consider stars

    Return
    -------
    com: SimArray
        Center of mass (in AU for each x,y,z component)
    """
    com = pynbody.array.SimArray(np.zeros(3),'au')	

    assert starFlag == True or gasFlag == True, "At least one flag must be true."    
    
    #If there's a cutoff, select gas particles with cylindrical radius less than the cutoff
    if cutoff != None:
        mask = gas['rxy'] < cutoff
        gas = gas[mask]
	
    if starFlag: #Include stars
        #Ensure binary
        assert len(stars) == 2
	
        #Compute stellar mass, mass-weighted position
        starMass = stars[0]['mass'] + stars[1]['mass']
        starPos = stars[0]['pos'].in_units('au')*stars[0]['mass']
        starPos += stars[1]['pos'].in_units('au')*stars[1]['mass']
	
        if gasFlag:
            #Compute, return total center of mass
            com = starPos + np.sum(gas['pos'].in_units('au')*gas[0]['mass'],axis=0)
            com /= (starMass + np.sum(gas['mass']))
            return com
        else:
            com = (starPos/starMass)
            return com
          
    else: #No stars, just gas
         com = np.sum(gas['pos'].in_units("au")*gas[0]['mass'],axis=0)
         com /= np.sum(gas['mass'])
         return com

#end function 
 

def computeVelocityCOM(s,cutoff=None,starFlag=True,gasFlag=True):
    """
    Given pynbody star and gas arrays, compute the center of mass velocity
    for the entire specified system.

    Parameters
    ----------
    s : pynbody snapshot
    cutoff : float
        radius at which you only consider objects interior to it [AU]
    starFlag : bool
        whether or not to consider stars
    gasFlag : bool
        whether or not to consider gas

    Returns
    -------
    Center of mass velocity: SimArry
        in AU for each vx,vy,vz component

    Note: a lot of "strip_units" commands included to prevent throwing weird value errors.  As long as all masses
    are in solar masses and positions in AU before this is run, you won't have any problems.	
    """
    stars = s.stars
    gas = s.gas
    
    com = pynbody.array.SimArray(np.zeros(3),'km s**-1')
    
    assert starFlag == True or gasFlag == True, "At least one flag must be true."    
    
    #If there's a cutoff, select gas particles with cylindrical radius less than the cutoff
    if cutoff != None:
        mask = gas['rxy'] < cutoff
        gas = gas[mask]
	
    if starFlag: #Include stars
        #Ensure binary
        assert len(stars) == 2
	
        #Compute stellar mass, mass-weighted position
        starMass = stars[0]['mass'] + stars[1]['mass']
        starPos = stars[0]['vel'].in_units('km s**-1')*stars[0]['mass']
        starPos += stars[1]['vel'].in_units('km s**-1')*stars[1]['mass']
	
        if gasFlag:
            #Compute, return total center of mass
            com = starPos + np.sum(gas['vel'].in_units('km s**-1')*gas[0]['mass'],axis=0)
            com /= (starMass + np.sum(gas['mass']))
            return com
        else:
            com = (starPos/starMass)
            return com
          
    else: #No stars, just gas
         com = np.sum(gas['vel'].in_units("km s**-1")*gas[0]['mass'],axis=0)
         com /= np.sum(gas['mass'])
         return com

#end function

def calcDiskRadialBins(s,r_in=0,r_out=0,bins=50):
    """
    Cleanly partitions disk into radial bins and returns the bin edges and central bin values.  Note, default
    ndim = 2 so all bins are in 2D plane (i.e. radius r is polar/cylindrical radius in xy plane which makes 
    sense for thin disks)

    Parameters
    ----------
    s: Pynbody snapshot
    r_in: float 
        Inner disk radius you'll consider (AU)
    r_out: float 
        Outer disk radius you'll consider (AU)
    bins: int
        # of bins 

    Returns
    --------
    r: numpy array 
        central radial bin values (AU)
    rBinEdges: numpy array 
        edges of radial bins
    """
    #Load data, compute semimajor axis and strip units
    x1 = s.stars[0]['pos']
    x2 = s.stars[1]['pos']
    v1 = s.stars[0]['vel']
    v2 = s.stars[1]['vel']
    m1 = s.stars[0]['mass']
    m2 = s.stars[1]['mass']
    s_a = AddBinary.calcSemi(x1, x2, v1, v2, m1, m2) #Units = au

    #Default r_in, r_out if none given
    if r_in == 0 and r_out == 0:
        r_in = 1.0*s_a
        r_out = 4.0*s_a

    #Bin gas particles by radius
    pg = pynbody.analysis.profile.Profile(s.gas,max=r_out,nbins=bins)
    r = pg['rbins'].in_units('au')    
    
    mask = (r > r_in) & (r < r_out) #Ensure you're not right on binary or too far out.
    r = r[mask]

    #Make nice, evenly spaced radial bins vector
    rBinEdges = np.linspace(np.min(r),np.max(r),bins+1)
 
    #Create compute center of radial bins
    r = 0.5 * (rBinEdges[1:] + rBinEdges[:-1])

    return r, rBinEdges

#end function

def calcNetTorque(stars,gas):
    """
    Given pynbody snapshot (Tipsy format)arrays of the stars and gas of a binary surrounded by a CB disk, 
    compute the net torque on the binary due to the CB disk.  	
    This function can be used to compute the net torque/mass due to any collection of gas (total disk, an annulus, etc) on 
    the stars.

    Parameters
    ----------
    stars, gas: pynbody-readable Tipsy snapshot arrays 
        of binary + CB disk.  Assumes units are in standard sim units (Msol,au...)
   
    Returns
    -------
    net torque: numpy array 
        Net torque/mass vector (3D) acting on binary system in cgs.
    """
    #Ensure system is binary
    assert len(stars) == 2, "Only use for binary system."

    #Compute center of mass of entire binary-disk system
    com = computeCOM(stars,gas).in_units('cm')

    #Compute net force on primary star (index = 0)

    #Compute M*m/|x'-x|^3 in cgs 
    #grav = AddBinary.BigG*(stars[0]['mass']*gas['mass']/np.power(np.linalg.norm(stars[0]['pos']-gas['pos'],1),3))
    grav = AddBinary.G*(stars[0]['mass'].in_units('g')*gas['mass'].in_units('g'))
    grav /= np.power(np.linalg.norm(stars[0]['pos'].in_units('cm')-gas['pos'].in_units('cm'),1),3)

    #Scale that value by (x'-x) to make it a vector pointing to gas particles
    F1 = -1*stars[0]['pos'].in_units('cm') + gas['pos'].in_units('cm')
    #conv = (AddBinary.Msol*AddBinary.Msol)/(AddBinary.AUCM*AddBinary.AUCM)
    F1[:,0] *= grav#*conv)
    F1[:,1] *= grav#*conv)
    F1[:,2] *= grav#*conv)

    #Compute net force on stars due to gas (3 components)
    F1 = np.sum(F1,axis=0)
    
    #Compute net force on secondary star (index = 1)

    #Compute M*m/|x'-x|^3 in cgs 
    #grav = AddBinary.BigG*(stars[0]['mass']*gas['mass']/np.power(np.linalg.norm(stars[0]['pos']-gas['pos'],1),3))
    grav = AddBinary.G*(stars[1]['mass'].in_units('g')*gas['mass'].in_units('g'))
    grav /= np.power(np.linalg.norm(stars[1]['pos'].in_units('cm')-gas['pos'].in_units('cm'),1),3)

    #Scale that value by (x'-x) to make it a vector pointing to gas particles
    F2 = -1*stars[1]['pos'].in_units('cm') + gas['pos'].in_units('cm')
    #conv = (AddBinary.Msol*AddBinary.Msol)/(AddBinary.AUCM*AddBinary.AUCM)
    F2[:,0] *= grav#*conv)
    F2[:,1] *= grav#*conv)
    F2[:,2] *= grav#*conv)

    #Compute net force on stars due to gas (3 components)
    F2 = np.sum(F2,axis=0)    
    
    """
    #Compute net force on the secondary star (index = 1) in cgs
    grav = AddBinary.BigG*(stars[1]['mass']*gas['mass']/np.power(np.linalg.norm(stars[1]['pos']-gas['pos'],1),3))

    #Scale that value by (x'-x) to make it a vector in cgs units
    conv = (AddBinary.Msol*AddBinary.Msol)/(AddBinary.AUCM*AddBinary.AUCM)
    F2 = -1*stars[1]['pos'] + gas['pos']
    F2[:,0] *= (grav*conv)
    F2[:,1] *= (grav*conv)
    F2[:,2] *= (grav*conv)

    #Compute net force
    F2 = np.sum(F2,axis=0)
    """

    #Compute the center of mass distances (au->cm)
    r1 = (stars[0]['pos'].in_units('cm') - com)#*AddBinary.AUCM
    r2 = (stars[1]['pos'].in_units('cm') - com)#*AddBinary.AUCM

    #Compute torque per unit mass in cgs
    tau1 = np.cross(r1,F1)
    tau2 = np.cross(r2,F2)
    #netTau = tau1/(stars[0]['mass']*AddBinary.Msol) + tau2/(stars[1]['mass']*AddBinary.Msol)
    netTau = tau1/(stars[0]['mass'].in_units('g')) + tau2/(stars[1]['mass'].in_units('g'))

    return np.asarray(netTau)
    
#end function

def torqueVsRadius(s,rBinEdges):
    """
    Takes in pynbody snapshot s for a binary system with a CB disk 
    returns torque per unit mass vs radius and the approximate radius of the bin where
    that torque was calculated.  Note, I only care about the z component of the torque
    since this is a circumbinary disk system 
    Note: This function is best for returning proper disk radial bins.
	
    Parameters
    ----------
    s: pynbody snapshot of binary + CB disk
    Bins: int
        number of radial bins to make the calculation
    r_in, r_out: floats
        Inner and outer radii over which torque is calculated (au)

    Returns
    -------
    tau: numpy array
        Torque per unit mass as function of radius (cgs vs au)
    """
    #Create array to hold torque as function of radius
    tau = np.zeros((len(rBinEdges)-1,3))

    #For a given radius, put gas particles in bins where s.gas['r'] is in au and pynbody units are stripped
    for i in range(0,len(rBinEdges)-1):
        rMask = np.logical_and(s.gas['rxy'].in_units('au') > rBinEdges[i], s.gas['rxy'].in_units('au') < rBinEdges[i+1])
        tau[i] = np.asarray(calcNetTorque(s.stars,s.gas[rMask])) 
        
    return tau

#end function

def calcDeDt(stars,tau):
    """
    Calculates the change in binary orbital eccentricity over time at each radial bin due to
    the torque/mass from the surrounding CB disk.

    Parameters
    ----------
    stars: pynbody stars object	(sim units)
    tau: torque/mass on binary due to CB disk during a given snapshot (cgs)
    tau is of the form tau[radii,(x,y,z) components]

    Output:
    (de/dt) at each radius (unitless/second)
    """
    #Strip units from all inputs and convert to cgs
    #Ensure units are in cgs
    x1 = stars[0]['pos'].in_units('cm')
    x2 = stars[1]['pos'].in_units('cm')
    v1 = stars[0]['vel'].in_units('cm s**-1')
    v2 = stars[1]['vel'].in_units('cm s**-1')
    m1 = stars[0]['mass'].in_units('g')
    m2 = stars[1]['mass'].in_units('g')

    #Relative position vector in cgs
    r = x1 - x2
    magR = SimArray(np.linalg.norm(r),'cm')

    #Compute standard gravitational parameter in cgs
    mu = AddBinary.G*(m1+m2)

    #Relative velocity vector in cgs
    v = v1 - v2
    magV = SimArray(np.linalg.norm(v),'cm s**-1')

    #Compute specific orbital energy
    eps = (magV*magV/2.0) - (mu/magR)

    #Compute specific angular momentum vector
    h = np.cross(r,v)
    magH = SimArray(np.linalg.norm(h),'cm**2 s**-1')

    #Calculate change in e vs time due to z component of torque/mass (dh/dt for specific angular momentum h)
    dedt = np.asarray((2.0*eps*magH/(mu*mu))*tau[:,2])
    dedt *= 1.0/np.sqrt(1.0 + (2.0*eps*magH*magH/(mu*mu)))

    return dedt #(unitless/second)

#end function

def findCBResonances(s,r,r_min,r_max,m_max=4,l_max=4,bins=50):
    """
    Given Tipsy snapshot, computes the resonances of disk on binary as a function of orbital angular frequency omega.
    Disk radius, in au, is convered to angular frequency which will then be used to compute corotation 
    and inner/outer Lindblad resonances.
   
   Note: r given MUST correspond to r over which de/dt was calculated.  Otherwise, scale gets all messed up
   
   !!! NOTE: This function is awful and deprecated --- do NOT use it.  Instead, use calc_LB_resonance !!!
 
     Parameters
     ----------
     s: Tipsy-format snapshot
     r: array
         radius array over which de/dt was calculated
     r_min,r_max: floats
         min/maximum disk radius for calculations (au)
     bins: int
         number of radial bins to calculate over
     m_max,l_max: ints
         maximum orders of (m,l) LR

    Returns
    -------
    Orbital frequency: numpy array
        for corotation and inner/outer resonances and radii as float and numpy arrays
    """
    stars = s.stars
    gas = s.gas

    m_min = 1 #m >=1 for LRs, CRs
    l_min = 1 #l >=1 for LRs, CRs

    #Compute binary angular frequency
    x1 = stars[0]['pos']
    x2 = stars[1]['pos']
    v1 = stars[0]['vel']
    v2 = stars[1]['vel']
    m1 = stars[0]['mass']
    m2 = stars[1]['mass']
     
    a = strip_units(AddBinary.calcSemi(x1, x2, v1, v2, m1, m2))
    omega_b = 2.0*np.pi/AddBinary.aToP(a,m1+m2) #In units 1/day

    #Make r steps smaller for higher accuracy
    r_arr = np.linspace(r.min(),r.max(),len(r)*10)

    #Compute mass of disk interior to given r
    mask = np.zeros((len(gas),len(r_arr)),dtype=bool)
    m_disk = np.zeros(len(r_arr))
    for i in range(0,len(r_arr)):
        mask[:,i] = gas['rxy'] < r_arr[i]
        m_disk[i] = np.sum(gas['mass'][mask[:,i]])

    #Compute omega_disk in units 1/day (like omega_binary)
    omega_d = 2.0*np.pi/AddBinary.aToP(r_arr,m1+m2+m_disk)
        
    #Compute kappa (radial epicycle frequency = sqrt(r * d(omega^2)/dr + 4*(omega^2))
    o2 = omega_d*omega_d
    dr = r_arr[1] - r_arr[0]
    #dr = (r.max()-r.min())/float(bins) #Assuming r has evenly spaced bins!
    drdo2 = np.gradient(o2,dr) #I mean d/dr(omega^2)
    kappa = np.sqrt(r_arr*drdo2 + 4.0*o2)
   
    #Allocate arrays for output 
    omega_Lo = np.zeros((m_max,l_max))
    omega_Li = np.zeros((m_max,l_max))
    o_c = np.zeros(l_max)   
 
    #Find resonance angular frequency
    for m in range(m_min,m_max+1):
        for l in range(l_min,l_max+1):
            outer = omega_d + (float(l)/m)*kappa
            inner = omega_d - (float(l)/m)*kappa
            omega_Lo[m-m_min,l-l_min] = omega_d[np.argmin(np.fabs(omega_b-outer))]
            omega_Li[m-m_min,l-l_min] = omega_d[np.argmin(np.fabs(omega_b-inner))]

            #Find corotation resonance where omega_d ~ omega_b
            o_c[l-l_min] = omega_d[np.argmin(np.fabs(omega_d-omega_b/float(l)))]

    #Rescale omega_d, kappa to be of length bins again
    omega_d = np.linspace(omega_d.min(),omega_d.max(),bins)
    kappa = np.linspace(kappa.min(),kappa.max(),bins)
    return omega_Li, omega_Lo, o_c, omega_d, kappa

#end function

def calc_LB_resonance(s,m_min=1,m_max=3,l_min=1,l_max=3):
    """
    Computes the locations of various Lindblad Resonances in the disk as a 
    function of binary pattern speed.
    
     Parameters
     ----------
     s : Tipsy-format snapshot
     m_min, l_min : ints
         minimum orders of (m,l) LR
     m_max,l_max : ints
         maximum orders of (m,l) LR

    Returns
    -------
    OLR, ILR, CR: numpy arrays
        location in AU of (m,l)th order Lindblad resonances
    """

    #Compute binary angular frequency in 1/day
    x1 = s.stars[0]['pos']
    x2 = s.stars[1]['pos']
    v1 = s.stars[0]['vel']
    v2 = s.stars[1]['vel']
    m1 = s.stars[0]['mass']
    m2 = s.stars[1]['mass']
    omega_b = 2.0*np.pi/AddBinary.aToP(AddBinary.calcSemi(x1,x2,v1,v2,m1,m2),m1+m2)
    guess = 0.05 #fsolve initial guess parameter

    #Allocate space for arrays
    OLR = np.zeros((m_max,l_max))
    ILR = np.zeros((m_max,l_max))
    CR = np.zeros(l_max)

    #Define resonance functions
    def OLR_func(omega_d, *args):
        m = args[0]
        l = args[1]
        omega_b = args[2]
        
        return omega_d*(1.0 + float(l)/m) - omega_b
        
    #end function        
        
    def ILR_func(omega_d, *args):
        m = args[0]
        l = args[1]
        omega_b = args[2]
        
        return omega_d*(1.0 - float(l)/m) - omega_b        
        
    #end function

    def CR_func(omega_d, *args):
        l = args[0]
        omega_b = args[1]
        
        return omega_d - omega_b/float(l)

    #end function

    for m in range(m_min,m_max+1):
        for l in range(l_min,l_max+1):
            OLR[m-m_min,l-l_min] = fsolve(OLR_func,guess,args=(m,l,omega_b)) 
            ILR[m-m_min,l-l_min] = fsolve(ILR_func,guess,args=(m,l,omega_b))
            CR[l-l_min] = fsolve(CR_func,guess,args=(l,omega_b))
            
    #Convert from 1/day -> au
    OLR = AddBinary.pToA(2.0*np.pi/OLR,m1+m2)
    ILR = AddBinary.pToA(2.0*np.pi/ILR,m1+m2)
    CR = AddBinary.pToA(2.0*np.pi/CR,m1+m2)        
    
    return OLR, ILR, CR
    
#end function

def calcCoMVsRadius(s,rBinEdges,starFlag=False):
    """
    Calculates the system's center of mass as a function of radius.  At a given radius r, use the total enclosed
    mass of the star(s) and gas to compute the center of mass (CoM).
	
    Parameters
    ----------
    s: Tipsy-format snapshot readable by pynbody
    rBinEdges: array
        edges of the array of radii in xy plane
    starFlag: bool
        whether or not to consider stars in center of mass calculation

    Returns
    -------
    com: array
        Numpy array of len(r) * 3 containing location of CoM in Cartesian coordinates.
    """
    stars = s.stars
    gas = s.gas
    com = np.zeros((len(rBinEdges)-1,3))
	
    if starFlag: #Include stars in center of mass calculation
        #Loop through radial points, select gas within that r, calc CoM
        for i in range(0,len(rBinEdges)-1):
            com[i,:] = computeCOM(stars,gas,cutoff=rBinEdges[i],starFlag=starFlag)
    else: #Gas disk only
        for i in range(0,len(rBinEdges)-1):
            com[i,:] = computeCOM(stars,gas,cutoff=rBinEdges[i],starFlag=starFlag)
    
    return com

#end function

def calcPoissonVsRadius(s,rBinEdges):
	"""
	Given a tipsy snapshot and radial bins, compute the Poisson noise, r/sqrt(N_particles), in each radial bin.
	"""	
	gas = s.gas
	poisson = np.zeros(len(rBinEdges)-1)	
	
	for i in range(0,len(rBinEdges)-1):
		rMask = np.logical_and(gas['rxy'].in_units('au') > rBinEdges[i], gas['rxy'].in_units('au') < rBinEdges[i+1])
		N = len(gas[rMask])
		r = (rBinEdges[i] + rBinEdges[i+1])/2.0
		poisson[i] = r/np.sqrt(N)
	
	return poisson
	
#end function
	
def calcQ(cs,kappa,sigma):
    """
    Compute the Toomre Q parameter for a gaseous disk.  Implimented here since pynbody calculates it for a
    stellar disk.
    Q = (c_s * kappa)/(pi*G*Sigma) > 1 -> axisymmetric stability
    
    Parameters
    ----------    
    c_s: float
        sound speed (cm/s)
    Kappa: float
        radially epicycle frequency (1/s).  Can also be Omega for a quasi-Keplerian disk
    Sigma: float
        surface density at some radius (g/cm^2)	
	
    Output:
    Q: unitless
        Toomre parameter
    """
	
    return (cs*kappa)/(AddBinary.BigG*np.pi*sigma)
	
#end function
	
def calcQVsRadius(s,a_in,a_out,bins):
    """
    Given a tispy snapshot, compute the Toomre Q parameter at each radial point.
    
    Parameters
    ----------
    s: Tipsy snapshot
    a_in, a_out: floats
        Minimum and maximum radial points on which to calculate the profile [AU]	
    bins: int
        number of radial bins	
	
    Returns
    -------
    r: array, AU
    Q: array, unitless
    """
    #Derive quantities in correct units
    p = pynbody.analysis.profile.Profile(s.gas,max=a_out,min=a_in,nbins=bins)
    sigma = p['density'].in_units('g cm**-2')
    kappa = p['omega'].in_units('s**-1')
    cs = p['cs'].in_units('cm s**-1')
    r = p['rbins']
	
    return r, calcQ(cs,kappa,sigma)
	
#end function
	
def calcStableSigma(r,rd,Mstar,Mdisk,Q):
    """
    Compute the surfance density sigma_0 such that if the quasi-keplerian disk
    had a surface density of sigma > sigma_0 at an OLR, the disk would be 
    unstable to a m=1 mode.  Condition comes from eqn. 110 in Shu 1990.
    Note: Assumes (b_n - c)^2 ~ 1 as authors did.    
    
    Parameters
    -----------
    r: array
        radii of OLR [AU]
    rd: float
        Maximum radially extent of the disk [AU]
    Mstar, Mdisk: floats
        Masses of the central star(s) and disk, respectively [Msol]
    Q: float
        Toomre Q stability parameter evalutated at rd
    
    Returns
    -------
    sigma_0: array
        critical surfance density [Msol/AU^2]
    """    
    sigma_0 = 3.0*(Mstar + Mdisk)/(8.0*np.pi*np.pi*r*r)
    sigma_0 *= np.power(r/rd,3.0)
    sigma_0 *= np.sqrt(1.0 + 4.0*Q*Q*(np.power(rd/r,3.0) - np.power(rd/r,1.5)))
    
    return sigma_0
    
def orbElemsVsRadius(s,rBinEdges,average=False):
    """
    Computes the orbital elements for disk particles about a binary system in given radial bins.
    Assumes center of mass has v ~ 0

    Parameters
    ----------

    s: Tipsy snapshot
    rBinEdges: numpy array
        Radial bin edges [AU] preferably calculated using binaryUtils.calcDiskRadialBins
    average: bool
        True -> average over all particles in bin, false -> randomly select 1 particle in bin
        
    Returns
    -------
    orbElems: numpy array
        6 x len(rBinEdges) - 1 containing orbital elements at each radial bin
        as e, a, i, Omega, w, nu
    """
    
    #Read snapshot and pull out values of interest
    stars = s.stars
    gas = s.gas    
    M = np.sum(stars['mass'])
    zero = SimArray(np.zeros(3).reshape((1, 3)),'cm s**-1') 
    orbElems = np.zeros((6,len(rBinEdges)-1))    
    
    #Gas orbiting about system center of mass
    com = computeCOM(stars,gas)
   
    #Loop over radial bins calculating orbital elements
    for i in range(0,len(rBinEdges)-1):
        if average: #Average over all gas particles in subsection
            rMask = np.logical_and(gas['rxy'].in_units('au') > rBinEdges[i], gas['rxy'].in_units('au') < rBinEdges[i+1])
            if i > 0:
                #Include mass of disk interior to given radius
                mass = M + np.sum(gas[gas['rxy'] < rBinEdges[i]]['mass'])
            else:
                mass = M
            N = len(gas[rMask])
            g = gas[rMask]
            if N > 0:
                orbElems[:,i] = np.sum(AddBinary.calcOrbitalElements(g['pos'],com,g['vel'],zero,mass,g['mass']),axis=-1)/N
            else: #If there are no particles in the bin, set it as a negative number to mask out later
                orbElems[:,i] = -1.0
        else: #Randomly select 1 particle in subsection for calculations
            rMask = np.logical_and(gas['rxy'].in_units('au') > rBinEdges[i], gas['rxy'].in_units('au') < rBinEdges[i+1])
            if i > 0:            
                mass = M + np.sum(gas[gas['rxy'] < rBinEdges[i]]['mass'])
            else:
                mass = M
            g = gas[rMask]            
            index = np.random.randint(0,len(g))
            particle = g[index]
            orbElems[:,i] = AddBinary.calcOrbitalElements(com,particle['pos'],zero,particle['vel'],mass,particle['mass'])
            
    return orbElems
    
#end function    
    
def diskPrecession(s,radius):
    """
    Computes the precession of the disk due to the binary quadrupole moment.
    The precessions considered are kappa_r and kappa_z corresponding to the
    precession of the argument of periapsis and longitude of th ascending node,
    respectively.
    
    Precssion frequency: Omega_p = Omega - Kappa
    Omega = sqrt((G*mu/r^3)*(1 + 3*alpha/r^2)) == orbital frequency
    
    Parameters
    ----------
    
    s: Tipsy snapshot
    r: numpy array
        array of radial bins centers [AU]
        
    Returns:
    -------
    
    Kappa array: numpy array
        2 x len(rBinEdges)-1 array containing precession at each radial point in 1/s
    """
    #Compute relevant frequencies
    alpha = SimArray(0.25,'au**2')
    r = SimArray(radius,'au')
    grav = SimArray(4.0*np.pi**2,'au**3 yr**-2 Msol**-1')
    mu = (s.stars[0]['mass'] * s.stars[1]['mass'])/np.sum(s.stars['mass'])
    omega = np.sqrt((grav*mu/np.power(r,3)) * (1.0 + 3.0*alpha/np.power(r,2)))
    kappa_r = np.sqrt((grav*mu/np.power(r,3)) * (1.0 - 3.0*alpha/np.power(r,2)))
    kappa_z = np.sqrt((grav*mu/np.power(r,3)) * (1.0 + 9.0*alpha/np.power(r,2)))
    
    #Compute precession. > 0 -> preccesion, < 0 -> recession
    omega_p = SimArray(np.zeros((2,len(r))),'yr**-1')
    omega_p[0,:] = omega - kappa_r
    omega_p[1,:] = omega - kappa_z
    
    return omega_p
    
#end function
    
def diskAverage(s,r_out,bins=50,avgFlag=True):
    """
    Computes the accretion disk mass-averaged for x via the following equation:
    integral of 2*pisigma*x*r*dr / integral of 2*pi*sigma*r*dr.
    Sigma, e,a... calculated on the fly to ensure that they are all evaluated at
    the same location.
    
    Parameters
    ----------
    s : tipsy snapshot
    r_out : float
        outer radii for averaging region. If none, use entire disk
    bins : int
        how many radial bins to calculate quantities over
    avgFlag : bool
        whether or not to average over all particles in a radial bin
        for orbital element calculation
        
    Returns
    -------
    y : list
        disk-averaged Keplerian orbital elements [e,a,i,Omega,w,nu] in AU, degrees (depending on unit)
    """    
    
    #Generate radial surface density profile
    #Begin by subtracting off the center of mass position
    #cm = computeCOM(s.stars,s.gas,cutoff=r_out,starFlag=True)
    #s['pos'] -= cm
    r = s.gas['rxy'].in_units('au')

    #Particle mass
    m_gas = s.gas['mass'][[0]]
    
    N, rBinEdges = np.histogram(r, bins=bins,range=(r.min(),r_out))
    rBinEdges = SimArray(rBinEdges,'au')
    r = (rBinEdges[1:] + rBinEdges[:-1])/2.0
    dr = rBinEdges[[1]] - rBinEdges[[0]]
    
    #Compute quantities to integrate
    sig = N*m_gas/(2.0*np.pi*r*dr)
    #s['pos'] += cm
    x = orbElemsVsRadius(s,rBinEdges,average=avgFlag)
    
    #Take correct cuts of data
    mask = r < r_out
    r = r[mask]
    sig = sig[mask]
    x = x[:,mask]
    
    #Compute total mass in region    
    denom = np.trapz(sig*r,r)
    
    #Compute mass-averaged x in region
    num = np.trapz(sig*r*x[:],r)

    return num/denom
    
#end function

def forcedEccentricity(binary_sys,r):
    """
    Given a binary class object and an array of radial points in the disk, 
    compute the forced eccentricity defined by Moriwaki et al. 2004 
    eqn 9 to first order.  Extra factor of 2 to give e_pumped instead
    of e_forced.  Note: This only applies when e_binary != 0 and when
    m2/(m1 + m2) != 0.5 (i.e. only applies for eccentric, non-equal mass
    binary)
    
    Parameters
    ----------
    binary_sys : binary.Binary class object
    r : array
        array of radii in AU
        
    Returns
    -------
        e_forced : array
            array of len(r)
    """
    mu = binary_sys.m2/(binary_sys.m1 + binary_sys.m2)
    return (5./2.)*(1.0 - 2.0*mu)*binary_sys.e*binary_sys.a/r
    
#end function
