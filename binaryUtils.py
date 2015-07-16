"""
David Fleming
Utilities to process/interact with binary star system in ChaNGa sims
"""

# Imports
import numpy as np
import re
import AddBinary
import isaac
import pynbody
from scipy import interpolate

def angMomSearch(name,simUnits=False,flag="AngularMomentum="):
	"""
	Given the name of a file containing line dumps from ChaNGa, searches for appropriate angular momentum flag (see default)
	and tallys up total change in angular momentum

	Input:
	Name of input file (something.txt)
	Flag (String to search for; defaults to AngularMomentum=)
	simUnits (whether or not to put results in terms of sim units; defaults to False so cgs units used)

	Output:
	Total angular momentum accreted by all sink particles in user-specified units
		Typically used for binary star sim so I only care about those 2 sinks
	"""

	#Read in data
	data = np.genfromtxt(name,dtype="str")

	#Find lines with angular momentum in them with correct format
	mask = (np.core.defchararray.find(data,flag) > -1)
	L = data[mask]
 
	#Loop over values in array containing angular momentum lines
	res = 0
	for i in range(len(L)):
		res = res + np.asarray(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",L[i])).astype(float)

	#Convert res to cgs if user so wishes (i.e. simUnits == False)
	#Is read in in units: Msol*(28.something km/s)*AU

	if(simUnits == False):
		return res*(AddBinary.Msol*AddBinary.AUCM*1000*100*AddBinary.VEL_UNIT)
	else:
		return res

#end function

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
	r1 = np.asarray(isaac.strip_units(x1))*AddBinary.AUCM
	r2 = np.asarray(isaac.strip_units(x2))*AddBinary.AUCM
	v1 = np.asarray(isaac.strip_units(v1))*AddBinary.VEL_UNIT*100*1000
	v2 = np.asarray(isaac.strip_units(v2))*AddBinary.VEL_UNIT*100*1000
	m1 = np.asarray(isaac.strip_units(m1))*AddBinary.Msol
	m2 = np.asarray(isaac.strip_units(m2))*AddBinary.Msol
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
	
	Inputs:
	r: array of radial points
	array: array of quantity of interest (could be surface density?) as a function of r
	toFind: array value you're looking for
	num: resolution of search

	Output:
	critial_radius: radius at which array(critical_radius) = toFind (approximately)
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

def computeCOM(stars,gas,cutoff=None,starFlag=True):
	"""
	Given pynbody star and gas arrays, compute the center of mass for the entire specified system.

	Inputs: (pynbody objects!)
	stars: s.stars pynbody object
	gas: s.gas pynbody object
	cutoff: radius at which you only consider objects interior to it
	starFlag: whether or not to consider stars

	Output:
	Center of mass (in AU for each x,y,z component) as numpy array

	Note: a lot of "strip_units" commands included to prevent throwing weird value errors.  As long as all masses
	are in solar masses and positions in AU before this is run, you won't have any problems.	
	"""
	#If there's a cutoff, select gas particles with cylindrical radius less than the cutoff
	if cutoff != None:
		mask = gas['rxy'] < cutoff
		gas = gas[mask]
	
	if starFlag: #Include stars
		#Ensure binary
		assert len(stars) == 2
	
		#Compute stellar mass, mass-weighted position
		starMass = np.sum(stars['mass'])
		starPos = (stars[0]['pos']*isaac.strip_units(stars[0]['mass']) + stars[1]['pos']*isaac.strip_units(stars[1]['mass']))
	
		#Compute, return total center of mass
		return np.asarray((starPos + np.sum(gas['pos']*isaac.strip_units(np.mean(gas['mass']))))/np.sum(starMass+np.sum(gas['mass'])))
	else: #No stars, just gas
		return np.sum(gas['pos']*isaac.strip_units(np.mean(gas['mass'])))/np.sum(gas['mass'])

#end function

def calcDiskRadialBins(s,r_in=0,r_out=0,bins=50):
	"""
	Cleanly partitions disk into radial bins and returns the bin edges and central bin values.  Note, default
	ndim = 2 so all bins are in 2D plane (i.e. radius r is polar/cylindrical radius in xy plane which makes 
	sense for thin disks)

	Inputs:
	s: Pynbody snapshot
	r_in: Inner disk radius you'll consider (AU)
	r_out: Outer disk radius you'll consider (AU)
	bins: # of bins 

	Outputs:
	r: central radial bin values (AU)
	rBinEdges: edges of radial bins
	"""
	#Load data, compute semimajor axis and strip units
	x1 = isaac.strip_units(s.stars[0]['pos'])
	x2 = isaac.strip_units(s.stars[1]['pos'])
	v1 = isaac.strip_units(s.stars[0]['vel'])
	v2 = isaac.strip_units(s.stars[1]['vel'])
	m1 = isaac.strip_units(s.stars[0]['mass'])
	m2 = isaac.strip_units(s.stars[1]['mass'])
	s_a = AddBinary.calcSemi(x1, x2, v1, v2, m1, m2) #Units = au

	#Default r_in, r_out if none given
	if r_in == 0 and r_out == 0:
		r_in = 1.0*s_a
		r_out = 4.0*s_a

	#Bin gas particles by radius
	pg = pynbody.analysis.profile.Profile(s.gas,max=r_out,nbins=bins)
	r = isaac.strip_units(pg['rbins']) #Radius from origin in xy plane
	mask = (r > r.min()) & (r < r_out) #Ensure you're not right on binary or too far out.  Redundant, but whatever
	r = r[mask]

	#Make nice, evenly spaced radial bins vector
	rBinEdges = np.linspace(np.min(r),np.max(r),bins+1)
 
	#Create compute center of radial bins
	r = 0.5* (rBinEdges[1:] + rBinEdges[:-1])

	return r, rBinEdges

#end function

def calcNetTorque(stars,gas):
	"""
	Given pynbody snapshot (Tipsy format)arrays of the stars and gas of a binary surrounded by a CB disk, 
	compute the net torque on the binary due to the CB disk.  	
	This function can be used to compute the net torque/mass due to any collection of gas (total disk, an annulus, etc) on 
	the stars.

	Input:
	stars, gas: pynbody-readable Tipsy snapshot arrays of binary + CB disk.  Assumes units are in standard sim units (Msol,au...)
   
	Output:
	numpy array of Net torque/mass vector (3D) acting on binary system in cgs.
   
	"""
	#Ensure system is binary
	assert len(stars) == 2

	#Compute center of mass of entire binary-disk system
	com = computeCOM(stars,gas)


	#Compute net force on primary star (index = 0)

	#Compute M*m/|x'-x|^3 in cgs 
	grav = AddBinary.BigG*(stars[0]['mass']*gas['mass']/np.power(np.linalg.norm(stars[0]['pos']-gas['pos'],1),3))

	#Scale that value by (x'-x) to make it a vector pointing to gas particles
	F1 = -1*stars[0]['pos'] + gas['pos']
	conv = (AddBinary.Msol*AddBinary.Msol)/(AddBinary.AUCM*AddBinary.AUCM)
	F1[:,0] *= (grav*conv)
	F1[:,1] *= (grav*conv)
	F1[:,2] *= (grav*conv)

	#Compute net force on stars due to gas (3 components)
	F1 = np.sum(F1,axis=0)

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

	#Compute the center of mass distances (au->cm)
	r1 = (stars[0]['pos'] - com)*AddBinary.AUCM
	r2 = (stars[1]['pos'] - com)*AddBinary.AUCM

	#Compute torque per unit mass in cgs
	tau1 = np.cross(r1,F1)
	tau2 = np.cross(r2,F2)
	netTau = tau1/(stars[0]['mass']*AddBinary.Msol) + tau2/(stars[1]['mass']*AddBinary.Msol)

	return np.asarray(netTau)
    
#end function

def torqueVsRadius(s,rBinEdges):
	"""
	Takes in pynbody snapshot s for a binary system with a CB disk 
	returns torque per unit mass vs radius and the approximate radius of the bin where
	that torque was calculated.  Note, I only care about the z component of the torque
	since this is a circumbinary disk system 
	Note: This function is best for returning proper disk radial bins.
	
	Inputs:
	s: pynbody snapshot of binary + CB disk
	Bins: number of radial bins to make the calculation
	r_in, r_out: Inner and outer radii over which torque is calculated (au)

	Outputs:
	tau: Torque per unit mass as function of radius (cgs vs au)
	"""
	#Create array to hold torque as function of radius
	tau = np.zeros((len(rBinEdges)-1,3))

	#For a given radius, put gas particles in bins where s.gas['r'] is in au and pynbody units are stripped
	for i in range(0,len(rBinEdges)-1):
		rMask = np.logical_and(isaac.strip_units(s.gas['rxy']) > rBinEdges[i], isaac.strip_units(s.gas['rxy']) < rBinEdges[i+1])
		tau[i] = np.asarray(calcNetTorque(s.stars,s.gas[rMask])) 
        
	return tau

#end function

def calcDeDt(stars,tau):
	"""
	Calculates the change in binary orbital eccentricity over time at each radial bin due to
	the torque/mass from the surrounding CB disk.

	Inputs:
	stars: pynbody stars object	(sim units)
	tau: torque/mass on binary due to CB disk during a given snapshot (cgs)
		 tau is of the form tau[radii,(x,y,z) components]

	Output:
	(de/dt) at each radius (unitless/second)
	"""
	#Strip units from all inputs and convert to cgs
	x1 = np.asarray(isaac.strip_units(stars[0]['pos']))*AddBinary.AUCM
	x2 = np.asarray(isaac.strip_units(stars[1]['pos']))*AddBinary.AUCM
	v1 = np.asarray(isaac.strip_units(stars[0]['vel']))*1000*100*AddBinary.VEL_UNIT
	v2 = np.asarray(isaac.strip_units(stars[1]['vel']))*1000*100*AddBinary.VEL_UNIT
	m1 = np.asarray(isaac.strip_units(stars[0]['mass']))*AddBinary.Msol
	m2 = np.asarray(isaac.strip_units(stars[1]['mass']))*AddBinary.Msol

	#Relative position vector in cgs
	r = x1 - x2
	magR = np.linalg.norm(r)

	#Compute standard gravitational parameter in cgs
	mu = AddBinary.BigG*(m1+m2)

	#Relative velocity vector in cgs
	v = v1 - v2
	magV = np.linalg.norm(v)

	#Compute specific orbital energy
	eps = (magV*magV/2.0) - (mu/magR)

	#Compute specific angular momentum vector
	h = np.cross(r,v)
	magH = np.linalg.norm(h)

	#Calculate change in e vs time due to z component of torque/mass (dh/dt for specific angular momentum h)
	dedt = np.asarray((2.0*eps*magH/(mu*mu))*tau[:,2])
	dedt *= 1.0/np.sqrt(1.0 + (2.0*eps*magH*magH/(mu*mu)))

	return dedt #(unitless/second)

#end function

def estimateCBResonances(s,r_max,m_max=5,l_max=5,bins=2500):
	"""
	Given pynbody snapshot star and gas SimArrays, computes the resonances of disk on binary as a function of period.
	Disk radius, in au, is convered to angular frequency which will then be used to compute corotation and inner/outer Lindblad resonances.
	Assumption: Assumes m_disk << m_bin which holds in general for simulations considered
	For reference: Kappa, omega computed in ~ 1/day intermediate units.
	Uses approximations from Artymowicz 1994

	Inputs:
	stars,gas: Pynbody snapshot .star and .gas SimArrays (in au, Msol, etc)
	r_max: maximum disk radius for calculations (au)
	bins: number of radial bins to calculate over

	Output:
	Orbital frequency for corotation and inner/outer resonances as float and 2 arrays
	"""
	stars = s.stars
	#gas = s.gas

	#Compute binary angular frequency
	#Strip units from all inputs
	x1 = np.asarray(isaac.strip_units(stars[0]['pos']))
	x2 = np.asarray(isaac.strip_units(stars[1]['pos']))
	v1 = np.asarray(isaac.strip_units(stars[0]['vel']))
	v2 = np.asarray(isaac.strip_units(stars[1]['vel']))
	m1 = np.asarray(isaac.strip_units(stars[0]['mass']))
	m2 = np.asarray(isaac.strip_units(stars[1]['mass']))
	a = AddBinary.calcSemi(x1, x2, v1, v2, m1, m2)
	#omega_b = 2.0*np.pi/AddBinary.aToP(a,m1+m2

	#Find corotation resonance where omega_d ~ omega_b
	r_c = a #m=1 case
	o_c = 2.0*np.pi/AddBinary.aToP(r_c,m1+m2)
        
	#Find inner lindblad resonances for m = [m_min,m_max]
	#Lindblad resonance: omega = omega_pattern +/- kappa/m for int m > 1
	m_min = 1
	l_min = 1    

	omega_Lo = np.zeros((m_max-m_min,l_max-l_min))
	omega_Li = np.zeros((m_max-m_min,l_max-l_min))
   
	#Find resonance radii, convert to angular frequency
	for m in range(m_min,m_max):
		for l in range(l_min,l_max):
		#oTmp = find_crit_radius(r,omega_d-(kappa/(float(m))),omega_b,bins) #outer LR
			oTmp = np.power(float(m+1)/l,2./3.)*a
			omega_Lo[m-m_min,l-l_min] = 2.0*np.pi/AddBinary.aToP(oTmp,m1+m2)
		
		#iTmp = find_crit_radius(r,omega_d+(kappa/(float(m))),omega_b,bins) #inner LR
			iTmp = np.power(float(m-1)/l,2./3.)*a
			omega_Li[m-m_min,l-l_min] = 2.0*np.pi/AddBinary.aToP(iTmp,m1+m2)

	return omega_Li, omega_Lo, o_c #return inner, outer, co angular frequencies

#end function

def findCBResonances(s,r,r_min,r_max,m_max=4,l_max=4,bins=50):
	"""
	Given Tipsy snapshot, computes the resonances of disk on binary as a function of orbital angular frequency omega.
	Disk radius, in au, is convered to angular frequency which will then be used to compute corotation 
	and inner/outer Lindblad resonances.
   
	Note: r given MUST correspond to r over which de/dt was calculated.  Otherwise, scale gets all messed up
 
	Inputs:
	s: Tipsy-format snapshot
	r: radius array over which de/dt was calculated
	r_min,r_max: ,min/maximum disk radius for calculations (au)
	bins: number of radial bins to calculate over
	m_max,l_max: maximum orders of (m,l) LR

	Output:
	Orbital frequency for corotation and inner/outer resonances and radii as float and numpy arrays
	"""
	stars = s.stars
	#gas = s.gas

	m_min = 1 #m >=1 for LRs, CRs
	l_min = 1 #l >=1 for LRs, CRs

	#Compute binary angular frequency
	#Strip units from all inputs
	x1 = np.asarray(isaac.strip_units(stars[0]['pos']))
	x2 = np.asarray(isaac.strip_units(stars[1]['pos']))
	v1 = np.asarray(isaac.strip_units(stars[0]['vel']))
	v2 = np.asarray(isaac.strip_units(stars[1]['vel']))
	m1 = np.asarray(isaac.strip_units(stars[0]['mass']))
	m2 = np.asarray(isaac.strip_units(stars[1]['mass']))
	a = AddBinary.calcSemi(x1, x2, v1, v2, m1, m2)
	omega_b = 2.0*np.pi/AddBinary.aToP(a,m1+m2) #In units 1/day

	#Compute omega_disk in units 1/day (like omega_binary)
	omega_d = 2.0*np.pi/AddBinary.aToP(r,m1+m2)
        
	#Compute kappa (radial epicycle frequency = sqrt(r * d(omega^2)/dr + 4*(omega^2))
	o2 = omega_d*omega_d
	dr = (r.max()-r.min())/float(bins) #Assuming r has evenly spaced bins!
	drdo2 = np.gradient(o2,dr) #I mean d/dr(omega^2)
	kappa = np.sqrt(r*drdo2 + 4.0*o2)
   
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

	return omega_Li, omega_Lo, o_c, omega_d, kappa

#end function

def calcEccVsRadius(s,rBinEdges):
	"""
	Calculates the average circumbinary disk eccentricity in radial bins.
	Also general enough to work for a circumstellar disk.	

	Inputs:
	snap: Pynbody snapshot (to get star, gas components)
	r: array of radial bins to calculate on.
	bins: Number of bins to calculate on (length of r...included so I don't confuse myself)

	Outputs:
	ecc: Vector of len = len(r) containing disk eccentricity.
	"""
	ecc = np.zeros(len(rBinEdges)-1)

	for i in range(0,len(rBinEdges)-1):
		rMask = np.logical_and(isaac.strip_units(s.gas['rxy']) > rBinEdges[i], isaac.strip_units(s.gas['rxy']) < rBinEdges[i+1])
		
		#For non-zero bins
		if(np.sum(rMask) != 0):
			x1 = s.gas[rMask]['pos']
			zero = np.zeros(3)
			v1 = s.gas[rMask]['vel']
			m1 = s.stars[0]['mass']
			
			#For binary case
			if(len(s.stars) == 2):
				m2 = s.stars[1]['mass']
			else:
				m2 = 0
		
			#Calculate average e
			ecc[i] = np.sum(AddBinary.calcEcc(x1,zero,v1,zero,m1,m2))/np.sum(rMask)
		
		#No gas in bin -> no eccentricity
		else:
			ecc[i] = 0

	return ecc

#end function

def calcCoMVsRadius(s,rBinEdges,starFlag=False):
	"""
	Calculates the system's center of mass as a function of radius.  At a given radius r, use the total enclosed
	mass of the star(s) and gas to compute the center of mass (CoM).  Ideally, I'd like to see the CoM be at
	[0,0,0] every time (or within a really small number of that).
	
	Inputs: 
	s: Tipsy-format snapshot readable by pynbody
	rBinEdges: edges of the array of radii in xy plane
	starFlag: bool for whether or not to consider stars in center of mass calculation

	Output:
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
	Given a tipsy snapshot and radial bins, compute the Poisson noise, N_particles*sqrt(radius), in each radial bin.
	Expect a powerlaw trend since N_particles ~ Surface density profile.
	"""	
	gas = s.gas
	poisson = np.zeros(len(rBinEdges)-1)	
	
	for i in range(0,len(rBinEdges)-1):
		rMask = np.logical_and(isaac.strip_units(gas['rxy']) > rBinEdges[i], isaac.strip_units(gas['rxy']) < rBinEdges[i+1])
		N = len(gas[rMask])
		r = (rBinEdges[i] + rBinEdges[i+1])/2.0
		poisson[i] = np.sqrt(N)*r
	
	return poisson
	
#end function
	
def calcQ(cs,kappa,sigma):
	"""
	Compute the Toomre Q parameter for a gaseous disk.  Implimented here since pynbody calculates it for a
	stellar disk.
	Q = (c_s * kappa)/(pi*G*Sigma) > 1 -> axisymmetric stability
	Input: (all cgs)	
	c_s: sound speed (cm/s)
	Kappa: radially epicycle frequency (1/s).  Can also be Omega for a quasi-Keplerian disk
	Sigma: surface density at some radius (g/cm^2)	
	
	Output:
	Q: unitless
	"""
	
	return (cs*kappa)/(AddBinary.BigG*np.pi*sigma)
	
#end function
	
def calcQVsRadius(s,a_in,a_out,bins):
	"""
	Given a tispy snapshot, compute the Toomre Q parameter at each radial point.
	Input:
	s: Tipsy snapshot
	a_in, a_out: Minimum and maximum radial points on which to calculate the profile [AU]	
	bins = number of radial bins	
	
	Output:
	Q and radially profile (AU) it was calculated on.
	"""
	#Derive quantities in correct units
	p = pynbody.analysis.profile.Profile(s.gas,max=a_out,min=a_in,nbins=bins)
	sigma = p['density'].in_units('g cm**-2');
	kappa = p['kappa'].in_units('s**-1');
	cs = p['cs'].in_units('cm s**-1');
	r = p['rbins'];
	
	return r, calcQ(cs,kappa,sigma)
	
#end function
	
def calcStableSigma(r,rd,Mstar,Mdisk,Q):
    """
    Compute the surfance density sigma_0 such that if the quasi-keplerian disk
    had a surface density of sigma > sigma_0 at an OLR, the disk would be 
    unstable to a m=1 mode.  Condition comes from eqn. 110 in Shu 1990.
    Note: Assumes (b_n - c)^2 ~ 1 as authors did.    
    
    Input:
    r: radii of OLR [AU]
    rd: Maximum radially extent of the disk [AU]
    Mstar, Mdisk: Masses of the central star(s) and disk, respectively [Msol]
    Q: Toomre Q stability parameter evalutated at rd
    
    Output:
    sigma_0: critical surfance density [Msol/AU^2]
    """    
    
    sigma_0 = 3.0*(Mstar + Mdisk)/(8.0*np.pi*np.pi*r*r)
    sigma_0 *= np.power(r/rd,3.0)
    sigma_0 *= np.sqrt(1.0 + 4.0*Q*Q*(np.power(rd/r,3.0) - np.power(rd/r,1.5)))
    
    return sigma_0