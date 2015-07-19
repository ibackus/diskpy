"""
 -----------------------------------------------
 Some simple python code to be easily imported from python
 -----------------------------------------------
 """
import pynbody
SimArray = pynbody.array.SimArray
pb = pynbody

import copy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import warnings
import glob
import os
import datetime
import fnmatch
import logging

self_dir = os.path.dirname(os.path.realpath(__file__))
print os.path.realpath(__file__)

def snapshot_defaults(snapshot):
    """
    Applies various defaults to tipsy snapshots of protoplanetary disk
    simulations. These include:
    
        -Sets nice units
        -Calculates particle smoothing lengths using mass and rho (if available)
    
    Changes the snapshot in place
    """
    
    # Setup units
    snapshot['pos'].convert_units('au')
    snapshot['mass'].convert_units('Msol')
    snapshot['vel'].convert_units('km s**-1')
    snapshot['eps'].convert_units('au')
    snapshot.g['temp'].convert_units('K')
    
    # Calculate smoothing lengths
    if ('rho' in snapshot.g.loadable_keys()) or ('rho' in snapshot.g.keys()):
        
        snapshot.g['rho'].convert_units('Msol au**-3')
        
        if ~(np.any(snapshot.g['rho'] == 0)):
            
            snapshot.g['smooth'] = (snapshot.g['mass']/snapshot.g['rho'])**(1,3)
    
    return

def units_from_param(param):
    """
    Figures out the simulation units from a .param file
    
    **ARGUMENTS**
    
    param : str or param dict (see configparser)
        Simulation .param file or param dict loaded by configparser
        Can also be a list or numpy array of these in which case a list
        of units dicts is returned
    
    **RETURNS**
    
    units : dict
        A dictionary of the units used in the simulation, returned as
        pynbody units
    """
    
    # Define function to load the units from a given param
    def _load_units(param):
    
        # Load param if necessary
        if isinstance(param, str):
            
            param = configparser(param, 'param')
            
        # Universal G
        G = pynbody.units.G
        
        # Load units
        dKpcUnit = param['dKpcUnit']
        dMsolUnit = param['dMsolUnit']
        
        # Set up pynbody units
        m_unit = pynbody.units.Unit('{0} Msol'.format(dMsolUnit))
        l_unit = pynbody.units.Unit('{0} kpc'.format(dKpcUnit))
        t_unit = (l_unit**3/(G*m_unit))**(1,2)
        
        # Convert the time unit to something sensible
        years = t_unit.in_units('yr')
        t_unit = pynbody.units.Unit('{0} yr'.format(years))
        
        # Return
        outdict = {'l_unit':l_unit, 'm_unit':m_unit, 't_unit':t_unit}
        return outdict
        
    # Iterate over param if necessary
    if isinstance(param, (list, np.ndarray)):
        
        outlist = []
        
        for par in param:
            
            outlist.append(_load_units(par))
            
        return outlist
        
    else:
        
        # Not iterable
        return _load_units(param)

def kepler_pos(pos, vel, t, Mstar, order=10):
    """
    Estimate position at future time t assuming an elliptical keplerian orbit
    """
    
    G = SimArray(1.0, 'G')
    mu = G*Mstar
    r = np.sqrt((pos**2).sum())
    v = np.sqrt((vel**2).sum())
    # Calculate semi-major axis
    a = mu*r/(2*mu - v*v*r)
    a.convert_units(r.units)
    # Calculate eccentricity vector
    ecc = (v*v)*pos/mu - ((pos*vel).sum())*vel/mu - pos/r
    ecc.convert_units('1')
    # Calculate eccentricity
    e = float(np.sqrt((ecc**2).sum()))
    
    # Calculate initial eccentric anomaly
    # x1 = a*e^2 + r.e
    x1 = a*e**2 + (pos*ecc).sum()
    # y1 = |r x e| * sign(r.v)
    y1 = np.sqrt((np.cross(pos, ecc)**2).sum())
    y1 *= (pos*vel).sum()/abs((pos*vel).sum())
    E0 = np.arctan2(y1,x1)
    
    # Calculate mean anomaly
    M0 = E0 - e*np.sin(E0)
    a3 = np.power(a,3)
    M = (np.sqrt(mu/a3)*t).in_units('1') + M0
    
    # Calculate eccentric anomaly
    E = E0
    for i in range(order):
        
        E = M + e*np.sin(E)
        
    # Calculate (x1, y1) (relative to center of ellipse, not focus)
    x1 = (2*a - r) * np.cos(E)
    y1 = (2*a - r) * np.sin(E)
    
    # Transform to original coordinates
    
    x1hat = ecc/np.sqrt((ecc**2).sum())
    y1hat = np.cross(np.cross(pos, vel), ecc)
    y1hat /= np.sqrt((y1hat**2).sum())
    
    pos_f = (x1 - a*e)*x1hat + y1*y1hat
    return pos_f

def findfiles(filefilter='*', basedir='.'):
    """
    Recursively find files according to filefilter
    
    ** ARGUMENTS **
    
    filefilter : str
        Filter for finding files.  ie, '*.jpg' or 'file.txt'
    
    basedir : str
        Base directory to search.  Default is the current directory
        
    ** RETURNS **
    
    files : list
        A list of the full path to all files matching filefilter
    
    """
    
    matches = []
    
    for root, dirnames, filenames in os.walk(basedir):
        
        for filename in fnmatch.filter(filenames, filefilter):
            fname = os.path.join(root, filename)
            fname = os.path.realpath(fname)
            
            matches.append(fname)
            
    return matches


def walltime(filename):
    """
    Reads walltime information from a ChaNGa .log file.
    
    ** ARGUMENTS **
    
    filename : str
        Filename of the .log file to load
        
    ** RETURNS **
    
    wall_per_step : array
        Wall time per step in seconds
    """
    
    log_file = np.genfromtxt(filename, comments='#', delimiter=' ')
    wall_per_step = log_file[:,-1]
    walltime_total = datetime.timedelta(seconds = wall_per_step.sum())
    walltime_avg = datetime.timedelta(seconds = wall_per_step.mean())

    print 'Total walltime: '
    print str(walltime_total)
    print 'Average walltime per step:'
    print str(walltime_avg)
    
    return wall_per_step
    

def load_acc(filename, param_name = None, low_mem = True):
    """
    Loads accelerations from a ChaNGa acceleration file (.acc2), ignoring the
    star particle.
    
    IF param_name is None, a .param file is searched for, otherwise param_name
    should be a string specifying a .param file name
    
    IF no param_file is found, the defaults are used:
        length unit: AU
        mass unit  : Msol
        
    Setting low_mem=True decreases memory usage by about 2x but also increases
    readtime by about 2x
    """
    
    if param_name is None:
        
        prefix = filename.split('.')[0]
        
        param_list = glob.glob('*' + prefix +'*param')
        
        if len(param_list) > 0:
            
            param_name = param_list[0]
            
        elif len(glob.glob('*.param')) > 0:
            
            param_name = glob.glob('*.param')[0]
            
        else:
            
            warnings.warn('Could not find .param file.  Assuming default units')
            
    if param_name is not None:
        
        # If a param name is set or a param file has been found:
        print 'Loading param file: {}'.format(param_name)
        param = configparser(param_name, ftype='param')
        
    else:
        
        # Set the default parameters
        param = {}
        # Assume AU as length unit
        param['dKpcUnit'] = pynbody.units.au.ratio('kpc')
        # Assume mass units as Msol
        param['dMsolUnit'] = 1.0
        
    # Figure out units
    G = pynbody.units.G
    l_unit = param['dKpcUnit']*pynbody.units.kpc
    m_unit = param['dMsolUnit']*pynbody.units.Msol
    t_unit = ((l_unit**3) * G**-1 * m_unit**-1)**(1,2)
    a_unit = l_unit * t_unit**-2
    
    if low_mem:
    
        acc_file = open(filename, 'r')
        n_particles = int(acc_file.readline().strip())
        acc = SimArray(np.zeros(3*n_particles, dtype=np.float32), a_unit)
        
        for i, line in enumerate(acc_file):
            
            acc[i] = np.float32(line.strip())
            
        acc_file.close()
            
        return acc.reshape([n_particles, 3], order='F')[0:-1]
        
    else:
        
        # Load acceleration file as numpy array
        acc = np.genfromtxt(filename, skip_header=1, dtype=np.float32)
        n_particles = len(acc)/3
        
        # Reshape and make it a SimArray with proper units
        acc = SimArray(acc.reshape([n_particles, 3], order='F'), a_unit)
        
        return acc
    
def height(snapshot, bins=100, center_on_star=True):
    """
    Calculates the characteristic height (h) of a flared disk as a function
    of cylindrical radius (r).
    
    ** ARGUMENTS **
    
    snapshot : TipsySnap
        Simulation snapshot for a flared disk
    bins : int or array_like
        Specifies the bins to use.  If int, specifies the number of bins.  If 
        array_like, specifies the bin edges
    center_on_star : bool
        If true (DEFAULT), cylindrical r is calculated relative to the star
        
    ** RETURNS **
    
    r_edges : SimArray
        Radial bin edges used for calculating h.  Length N+1
    h : SimArray
        Height as a function of r, calculated as the RMS of z over a bin.
        Length N
    """
    # Center on star
    if center_on_star:
        
        star_pos = snapshot.s['pos'].copy()
        snapshot['pos'] -= star_pos
        
    else:
        
        star_pos = 0.0*snapshot.s['pos']
        
    # Calculate height
    r = snapshot.g['rxy']
    z2 = snapshot.g['z']**2
    r_edges, z2_mean, err = binned_mean(r, z2, bins=bins, ret_bin_edges=True)
    h = np.sqrt(z2_mean)
    
    # Add star_pos back to snapshot
    snapshot['pos'] += star_pos
    
    return r_edges, h
        
def sigma(snapshot, bins=100):
    """
    Calculates surface density vs r (relative to the center of mass)
    
    ** ARGUMENTS **
    
    snapshot : tipsy snapshot
    bins : int, list, array...
        Either the number of bins to use or the binedges to use
        
    ** RETURNS **
    
    sigma : SimArray
        Surface density as a function of r
    r_bins : SimArray
        Radial bin edges
    """
    
    # Begin by subtracting off the center of mass position
    cm = (snapshot['mass'][:,None] * snapshot['pos']).sum()/(snapshot['mass'].sum())
    snapshot['pos'] -= cm
    r = snapshot.g['rxy']
    # particle mass
    m_gas = snapshot.gas['mass'][[0]]
    
    N, r_bins = np.histogram(r, bins=bins)
    r_bins = match_units(r_bins, r.units)[0]
    r_center = (r_bins[1:] + r_bins[0:-1])/2
    dr = r_bins[[1]] - r_bins[[0]]
    
    sig = N*m_gas/(2*np.pi*r_center*dr)
    
    # Add star position back to positions
    snapshot['pos'] += cm
    
    return sig, r_bins
    
    
def Q2(snapshot, molecular_mass = 2.0, bins=100, max_height=None):
    
    # Physical constants
    kB = SimArray([1.0],'k')
    G = SimArray([1.0],'G')
    # Load stuff froms snapshot
    v = snapshot.g['vt']
    r = snapshot.g['rxy']
    z = snapshot.g['z']
    T = snapshot.g['temp']
    # Calculate sound speed for all particles
    m = match_units(molecular_mass,'m_p')[0]
    cs = np.sqrt(kB*T/m)
    # Calculate surface density
    sig_binned, r_edges = sigma(snapshot, bins)
    r_cent = (r_edges[1:]+r_edges[0:-1])/2
    sig_spl = extrap1d(r_cent, sig_binned)
    sig = SimArray(sig_spl(r), sig_binned.units)
    # Calculate omega (as a proxy for kappa)
    omega = v/r
    kappa = omega
    
    #Calculate Q for all particles
    print 'kappa',kappa.units
    print 'cs',cs.units
    print 'sigma', sig.units
    Q_all = (kappa*cs/(np.pi*G*sig)).in_units('1')
    
    # Use particles close to midplane
    if max_height is not None:
        
        dummy, h = height(snapshot, bins=r_edges)
        ind = np.digitize(r, r_edges) - 1
        ind[ind<0] = 0
        ind[ind >= (len(r_edges)-1)] = len(r_edges)-2
        mask = abs(z) < (max_height*h[ind])
        Q_all = Q_all[mask]
        r = r[mask]
        
    dummy, Q_binned, dummy2 = binned_mean(r, Q_all, binedges=r_edges)
    
    return r_edges, Q_binned
    
def kappa(f, bins=100):
    """
    Estimate the epicyclic frequency from velocity
    
    **ARGUMENTS**
    
    f : TipsySnap
        Simulation snapshot
    bins : int or array-like
        Either the number of bins to use or the bin edges
        
    **RETURNS**
    
    kappa : SimArray
        epicyclic frequency
    r_edges : SimArray
        binedges used
    """    
    # Require regular spacing of bins
    if not isinstance(bins, int):
        
        dr = bins[[1]] - bins[[0]]
        eps = np.finfo(bins.dtype).eps
        
        if not np.all(bins[1:] - bins[0:-1] <= dr + 1000*eps):
            
            raise ValueError, 'Bins not uniformly spaced'
            
    r = f.g['rxy']
    v = f.g['vt']
    
    r_edges, v_mean, dummy = binned_mean(r, v, bins=bins, ret_bin_edges=True)
    dummy, rv_mean, dummy2 = binned_mean(r, r*v, bins=r_edges)
    r_cent = (r_edges[1:] + r_edges[0:-1])/2
    dr = r_edges[[1]] - r_edges[[0]]
    drv_dr = np.gradient(rv_mean, dr)
        
    kappa = np.sqrt(2*v_mean*drv_dr)/r_cent
        
    return kappa, r_edges
    
def Q(snapshot, molecular_mass = 2.0, bins=100, max_height=None, \
use_velocity=False, use_omega=True):
    """
    Calculates the Toomre Q as a function of r, assuming radial temperature
    profile and kappa ~= omega
    
    ** ARGUMENTS **
    
    snapshot : tipsy snapshot
    molecular_mass : float
        Mean molecular mass (for sound speed).  Default = 2.0
    bins : int or array
        Either the number of bins or the bin edges
    use_velocity : Bool
        Determines whether to use the particles' velocities to calculate orbital
        velocity.  Useful if the circular orbital velocities are set in the
        snapshot.
    use_omega : Bool
        Default=True.  Use omega as a proxy for kappa to reduce noise
        
    ** RETURNS **
    
    Q : array
        Toomre Q as a function of r
    r_edges : array
        Radial bin edges
    """
    
    # Physical constants
    kB = SimArray([1.0],'k')
    G = SimArray([1.0],'G')
    # Calculate surface density
    sig, r_edges = sigma(snapshot, bins)
    # Calculate sound speed
    m = match_units(molecular_mass,'m_p')[0]
    c_s_all = np.sqrt(kB*snapshot.g['temp']/m)
    # Bin/average sound speed
    dummy, c_s, dummy2 = binned_mean(snapshot.g['rxy'], c_s_all, binedges=r_edges)
    
    if use_omega:
        # Calculate keplerian angular velocity (as a proxy for the epicyclic
        # frequency, which is a noisy calculation)
        if use_velocity:
            # Calculate directly from particle's velocity
            dummy, omega, dummy2 = binned_mean(snapshot.g['rxy'], \
            snapshot.g['vt']/snapshot.g['rxy'], binedges=r_edges)
        
        else:
            # Estimate, from forces, using pynbody
            p = pynbody.analysis.profile.Profile(snapshot, bins=r_edges)    
            omega = p['omega']
        
        kappa_calc = omega
        
    else:
        
        if use_velocity:
            # Calculate directly from particle's velocities
            kappa_calc, dummy = kappa(snapshot, r_edges)
            
        else:
            # Estimate, from forces, using pynbody
            p = pynbody.analysis.profile.Profile(snapshot, bins=r_edges)    
            kappa_calc = p['kappa']
            
    return (kappa_calc*c_s/(np.pi*G*sig)).in_units('1'), r_edges
    
def Q_eff(snapshot, molecular_mass=2.0, bins=100):
    """
    Calculates the effective Toomre Q as a function of r, assuming radial temp
    profile and kappa ~= omega and scaleheight << wavelength.  This assumption
    simplifies the calculation of Q_eff (where wavelength is the wavelength of
    the disturbances of interest)
    
    ** ARGUMENTS **
    
    snapshot : tipsy snapshot
    molecular_mass : float
        Mean molecular mass (for sound speed).  Default = 2.0
    bins : int or array
        Either the number of bins or the bin edges
        
    ** RETURNS **
    
    Qeff : array
        Effective Toomre Q as a function of r for scale height << wavelength
    r_edges : array
        Radial bin edges
    """
    # Physical constants
    kB = SimArray([1.0],'k')
    G = SimArray([1.0],'G')
    # Calculate surface density
    sig, r_edges = sigma(snapshot, bins)
    # Calculate keplerian angular velocity (as a proxy for the epicyclic
    # frequency, which is a noisy calculation)
    p = pynbody.analysis.profile.Profile(snapshot, bins=r_edges)    
    omega = p['omega']
    # Calculate sound speed
    m = match_units(molecular_mass,'m_p')[0]
    c_s_all = np.sqrt(kB*snapshot.g['temp']/m)
    # Bin/average sound speed
    dummy, c_s, dummy2 = binned_mean(snapshot.g['rxy'], c_s_all, binedges=r_edges)    
    # Calculate scale height
    dummy, h = height(snapshot, bins=r_edges, center_on_star=False)
    
    a = np.pi*G*sig
    b = (2*a*h/c_s**2).in_units('1')
    Q0 = (omega*c_s/a).in_units('1')
    
    return Q0 * np.sqrt(1 + b), r_edges
    
    return ((omega*c_s/a) * np.sqrt(1 + 2*a*h/c_s**2)).in_units('1'), r_edges
    
    
    
def strip_units(x):
    """
    Removes the units from a SimArray and returns as a numpy array.  Note
    that x is copied so that it is not destroyed
    
    x can be a single SimArray or a tuple or list of SimArrays
    
    If any of the inputs are single number, they are returned as a number
    
    USAGE:
    
    array = strip_units(SimArray)
    
    array1, array2, ... = strip_units([SimArray1, SimArray2, ...])
    """
    if isinstance(x, (tuple,list)):
        
        # loop through and assign output
        x_out = []
        
        for x_i in x:
            
            if np.prod(x_i.shape) == 1:
                # There is only one element in x_i.  Make sure to return it as
                # a number  (not an array)
                if np.sum(x_i.shape) == 0:
                    # This is a zero dimensional SimArray
                    x_out.append(x_i.tolist())
                else:
                    # This is 1 dimensional SimArray
                    x_out.append(x_i[0])

            else:
                
                #This is a multi-element SimArray
                x_out.append(np.asarray(x_i.tolist()))
        
    else:
        
        if np.prod(x.shape) == 1:
            # There is only one element in x_i.  Return as a number
            if np.sum(x.shape) == 0:
                # This is a 0 dimensional SimArray
                x_out = x.tolist()
            else:
                # This a 1 dimensional SimArray
                x_out = x[0]
                
        else:
            
            x_out = np.asarray(x.tolist())
        
    return x_out
    
def set_units(x, units):
    """
    Sets the units of x to units.  If x has units, they are ignored.
    Does not destroy/alter x
    
    USAGE:
    
    SimArray = set_units(x, units)
    
    SimArray1, SimArray2, ... = set_units([x1, x2, ...], units)
    
    SimArray1, SimArray2, ... = set_units([x1, x2, ...], [units1, units2, ...])
    """
    if isinstance(x, (tuple,list)):
        
        x_out = []
        
        if not isinstance(units, (tuple, list)):
            
            units = [units]*len(x)
        
        for i in range(len(x)):
            
            x_i = x[i]
            
            if pynbody.units.has_units(x_i):
                
                x_i_array = strip_units(x_i)
                x_out.append(SimArray(x_i_array, units[i]))
                
            else:
                
                x_out.append(SimArray(x_i, units[i]))
            
    else:
        
        if pynbody.units.has_units(x):
            
            x_array = strip_units(x)
            x_out = SimArray(x_array, units)
            
        else:
            
            x_out = SimArray(x, units)
        
    return x_out
    
def setup_param(param, snapshot=None, r_orb=1.0, n_orb=10.0, n_image=None, n_snap=100, \
n_check=None):
    """
    Sets up the following for a .param file:
        
        nSteps
        dDumpFrameStep
        iOutInterval
        iCheckInterval
        
    **ARGUMENTS**
    
    param : str or param_dict (see isaac.configparser, configsave)
        parameter file for the simulation, must already have dDelta and units
        set properly
        IF a str, assumed to be a filename
    snapshot : str or TipsySnap(see pynbody) or None
        Snapshot for the simulation.  Needed to estimate the outer orbital
        period.  
        IF a str, assumed to be a filename
        IF None, the file pointed to by param is used
    r_orb : float
        radius to calculate the outer orbital period at as a fraction of the
        radius of the farthest out particle.  Must be between 0 and 1
    n_orb : float
        number of outer orbital periods to run simulation for
    n_image : int or None
        Total number of frames to dump (ie, dDumpFrameStep)
        If None, defaults to n_snap
    n_snap : int
        Total number of simulation outputs
    n_check : int or None
        Total number of simulation checkpoints.  If None, defaults to n_snap
    """
    
    if (r_orb > 1) | (r_orb < 0):
        
        raise ValueError, 'r_orb must be between 0 and 1'
    
    if isinstance(snapshot, str):
        
        # A filename has been passed, not a tipsy snapshot
        snapshot = pynbody.load(snapshot)
        
    if isinstance(param, str):
        
        # A filename has been passed.  Load the dictionary
        param = configparser(param, 'param')
        
    else:
        
        # Copy so as to not overwrite the input dict
        param = copy.deepcopy(param)
        
    R_max = r_orb * snapshot.g['rxy'].max()
    M_star = snapshot.s['mass']
    
    # Read in .param stuff
    l_unit = '{} kpc'.format(param['dKpcUnit'])
    m_unit = '{} Msol'.format(SimArray(param['dMsolUnit'], 'Msol'))
    
    # Outer radius and star mass in simulation units
    r = float(R_max.in_units(l_unit))
    M = float(M_star.in_units(m_unit))
    
    # Calculate the number of time steps to use
    dt = param['dDelta']
    period = 2*np.pi*np.sqrt(r**3/M)
    N = int(np.round(n_orb * period/dt))
    param['nSteps'] = N
    
    # Calculate how often to output snapshots, frames, checkpoints
    if n_check is None:
        
        n_check = n_snap
        
    if n_image is None:
        
        n_image = n_snap
        
    param['dDumpFrameStep'] = int(N/n_image)
    param['iOutInterval'] = int(N/n_snap)
    param['iCheckInterval'] = int(N/n_check)
    
    return param
            
def make_submission_script(param_name, directory=None, nodes=1, walltime=12, changa='ChaNGa_uw_mpi', jobname='changasim', scriptname='subber.sh', backfill=True):
    """
    Creates a submission script for qsub.  This is highly platform dependent
    """
    
    # Set up simulation directory
    if directory is None:
        
        directory = os.getcwd()
        
    # Load param file        
    param = configparser(param_name, 'param')
    fprefix = param['achOutName']
    
    # Format walltime for qsub
    seconds = int(walltime*3600)
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    walltime_str = '{0:02d}:{1:02d}:{2:02d}'.format(h,m,s)
    
    # Format walltime for changa
    walltime_min = int(walltime*60)
    
    # Open submission script
    subber = open(scriptname,'w')
    
    # Write to submission script
    subber.write('#!/bin/bash\n\
#PBS -N {0}\n\
#PBS -j oe\n\
#PBS -m be\n\
#PBS -M ibackus@gmail.com\n\
#PBS -l nodes={1}:ppn=12,feature=12core\n\
#PBS -l walltime={2}\n\
#PBS -V\n'.format(jobname, nodes, walltime_str))
    
    if backfill:
        
        subber.write('#PBS -q bf\n')
        
    subber.write('module load gcc_4.4.7-ompi_1.6.5\n')
    subber.write('export MX_RCACHE=0\n')
    
    subber.write('workdir={0}\n'.format(directory))
    subber.write('cd $workdir\n')
    subber.write('changbin=$(which {0})\n'.format(changa))
    subber.write('if [ -e "lastcheckpoint" ]\n\
then\n\
    echo "lastcheckpoint exists -- restarting simulation..."\n\
    last=`cat lastcheckpoint`\n\
    mpirun --mca mtl mx --mca pml cm $changbin +restart {0}.chk$last +balancer MultistepLB_notopo -wall {2} $workdir/{1} >> $workdir/{0}.out 2>&1\n\
else\n\
    echo "lastcheckpoint doesnt exist -- starting new simulation..."\n\
    mpirun --mca mtl mx --mca pml cm $changbin -D 3 +consph +balancer MultistepLB_notopo -wall {2} $workdir/{1} >& $workdir/{0}.out\n\
fi\n\
'.format(fprefix, param_name, walltime_min))
    
    subber.close()
    
    # Make submission script executable
    os.system('chmod a+rwx {}'.format(scriptname))

def make_param(snapshot, filename=None):
    """
    Generates a default param dictionary.  Can be saved using isaac.configsave
    
    EXAMPLE
    
    snapshot = pynbody.load('snapshot.std')  # Load snapshot
    param_dict = isaac.make_param(snapshot)  # Make default param dict
    isaac.configsave(param_dict, 'snapshot.param', ftype='param') # Save
    
    Optionally, the user can set the snapshot filename manually
    """
    fname_def = os.path.join(self_dir, 'default.param')
    param = configparser(fname_def, ftype='param')
    
    if filename is not None:
        
        param['achInFile'] = filename
        param['achOutName'] = os.path.splitext(filename)[0]
        
    elif snapshot.filename != '<created>':
        
        param['achInFile'] = snapshot.filename
        param['achOutName'] = os.path.splitext(snapshot.filename)[0]
        
    # Set up the length units
    param['dKpcUnit'] = snapshot['pos'].units.ratio('kpc')
    # Set up the mass units
    param['dMsolUnit'] = snapshot['mass'].units.ratio('Msol')
    # Set the mean molecular mass
    param['dMeanMolWeight'] = snapshot.gas['mu'][0]
    
    return param
    
def make_director(sigma_min, sigma_max, r, resolution=1200, filename='snapshot'):
    """
    Makes a director dictionary for ChaNGa runs based on the min/max surface
    density, maximum image radius, and image resolution for a gaseous
    protoplanetary disk.  The created dictionary can be saved with
    isaac.configsave
    
    The method is to use an example director file (saved as default.director)
    which works for one simulation and scale the various parameters accordingly.
    default.director should have a commented line in it which reads:
        #sigma_max float
    where float is the maximum surface density of the simulation in simulation
    units.
    
    **ARGUMENTS**
    
    sigma_min : float
        The surface density that corresponds to 0 density on the image (ie the
        minimum threshold).  Required for setting the dynamic range
    sigma_max : float
        Maximum surface density in the simulation
    r : float
        Maximum radius to plot out to
    resolution : int or float
        Number of pixels in image.  The image is shape (resolution, resolution)
    filename : str
        prefix to use for saving the images.  Example: if filename='snapshot',
        then the outputs will be of form 'snapshot.000000000.ppm'
        
    **RETURNS**
    
    director : dict
        A .director dictionary.  Can be saved with isaac.configsave
    """
    # -----------------------------------------------------------
    # Parse defaults to get scale factor for c
    # -----------------------------------------------------------
    defaults = configparser(os.path.join(self_dir, 'default.director'))
    if '#sigma_max' not in defaults:
        
        raise KeyError,'Default .director file should have a line e.g. << #sigma_max 0.01 >>'
    
    sigma_max0 = defaults['#sigma_max']
    c0 = defaults['colgas'][3]
    n0 = defaults['size'][0]
    r0 = defaults['eye'][2]
    A = (c0 * float(n0)**2)/(sigma_max0 * r0**2)
    
    # -----------------------------------------------------------
    # Create new director dictionary
    # -----------------------------------------------------------
    director = copy.deepcopy(defaults)
    director.pop('#sigma_max', None)
    
    logscale_min = sigma_min/sigma_max
    
    if pynbody.units.has_units(logscale_min):
        
        logscale_min = float(logscale_min.in_units('1'))
        
    c = A * float(sigma_max * r**2 /float(resolution)**2)
    
    director['colgas'][3] = c
    director['size'] = [resolution, resolution]
    director['eye'][2] = r
    director['file'] = filename
    
    return director
    
def match_units(x, y):
    """
    Matches the units of x to y and returns x and y in the same units.
    
    IF x and y don't have units, they are unchanged
    
    IF one of x or y has units, the unitless quantity is returned as a 
    SimArray (see pynbody.array.SimArray) with the units of the other quantity.
    
    IF both have units, then an attempt is made to convert x into the units of
    y.  If this is not possible, an error is raised, for example if x is in
    units of 'au' and y is in units of 'Msol'
    
    x, y can be: scalar, array, SimArray, pynbody unit (eg pynbody.units.G),
        or a unit string (eg 'Msol a**-2')
    
    
    *** RETURNS ***
    
    x, y are returned as a tuple
    """
    # ----------------------------------------------
    # Check if either is a string
    # ----------------------------------------------
    if isinstance(x, str):
        
        x = SimArray(1.0, x)
        
    if isinstance(y,str):
        
        y = SimArray(1.0, y)
    
    # ----------------------------------------------
    # Check if one is a pynbody unit
    # ----------------------------------------------
    # If one is a named unit (eg pynbody.units.G), convert to SimArray
    if isinstance(x, pynbody.units.UnitBase):
        
        x = SimArray(1.0, x)
        
    if isinstance(y, pynbody.units.UnitBase):
        
        y = SimArray(1.0, y)
        
    # ----------------------------------------------
    # Check the units
    # ----------------------------------------------
    # If both have units, try to convert x to the units of y
    if (pynbody.units.has_units(x)) & (pynbody.units.has_units(y)):
        
        x_out = (x.in_units(y.units))
        y_out = y
    
    # If only x has units, make y a SimArray with the units of x
    elif (pynbody.units.has_units(x)):
        
        y_out = SimArray(y, x.units)
        x_out = x
        
    # If only y has units, make x a SimArray with the units of y
    elif (pynbody.units.has_units(y)):
        
        x_out = SimArray(x, y.units)
        y_out = y
    
    # Otherwise, neither has units
    else:
        
        x_out = x
        y_out = y
        
    # Try to copy so that changing x_out, y_out will not change x,y
    try: 
        
        x_out = x_out.copy()
        
    except AttributeError:
        
        pass
    
    try: 
        
        y_out = y_out.copy()
        
    except AttributeError:
        
        pass
    
    return x_out, y_out
    
def digitize_threshold(x, min_per_bin = 0, bins=10):
    
    """
    Digitizes x according to bins, similar to numpy.digitize, but requires
    that there are at least min_per_bin entries in each bin.  Bins that do not
    have enough entries are combined with adjacent bins until they meet the
    requirement.
    
    **ARGUMENTS**
    
    x : array_like
        Input array to be binned.  Must be 1-dimensional
    min_per_bin : int
        Minimum number of entries per bin.  Default = 0
    bins : int or sequence of scalars, optional
        [same as for np.histogram]
        If bins is an int, it defines the number of equal-width bins in the 
        given range (10, by default). If bins is a sequence, it defines the 
        bin edges, including the rightmost edge, allowing for non-uniform bin 
        widths.
        
    **RETURNS**
    
    A tuple containing:
    ind : array_like
        Indices of the bin each element of x falls into, such that:
        bin_edges[i] <= x[i] < bin_edges[i+1]
        (See np.digitize, this uses the same convention)
    bin_edges: array_like
        The edges of the bins
    """
    
    # Find number in each bin
    N, bin_edges = np.histogram(x, bins)
    
    if N.sum() < min_per_bin:
        
        raise RuntimeError,'Not enough particles within the bin range'
        
    n_bins = len(bin_edges) - 1
    
    # Find out which binedges to delete
    edge_mask = np.ones(len(bin_edges), dtype='bool')    
    
    for i in range(n_bins - 1):
        # Work forwards
        
        if N[i] < min_per_bin:
            
            # Set mask to not use the right bin edge
            edge_mask[i+1] = False
            # Combine the particles in current and next bin
            N[i] += N[i+1]
            N[i+1] = N[i]
            
    bin_mask = edge_mask[1:]
    N = N[bin_mask]
    bin_edges = bin_edges[edge_mask]
    edge_mask = np.ones(len(bin_edges), dtype='bool')
    n_bins = len(bin_edges) - 1
    
    for i in range(n_bins-1, 0, -1):
        # Work backwards
        
        if N[i] < min_per_bin:
            
            # Set mask to not use the left bin edge
            edge_mask[i] = False
            # Combine the particles in current and next bin
            N[i] += N[i-1]
            N[i-1] = N[i]
            
    bin_edges = bin_edges[edge_mask]
    ind = np.digitize(x, bin_edges)
    
    return ind, bin_edges

def binned_mean(x, y, bins=10, nbins=None, binedges = None, weights=None,\
weighted_bins=False, ret_bin_edges=False):
    """
    Bins y according to x and takes the average for each bin.  
    
    bins can either be an integer (the number of bins to use) or an array of
    binedges.  bins will be overridden by nbins or binedges

    Optionally (for compatibility reasons) if binedges is specified, the 
    x-bins are defined by binedges.  Otherwise the x-bins are determined by 
    nbins
    
    If weights = None, equal weights are assumed for the average, otherwise
    weights for each data point should be specified
    
    y_err (error in y) is calculated as the standard deviation in y for each
    bin, divided by sqrt(N), where N is the number of counts in each bin
    
    IF weighted_bins is True, the bin centers are calculated as a center of
    mass
    
    NaNs are ignored for the input.  Empty bins are returned with nans
    
    RETURNS a tuple of (bin_centers, y_mean, y_err) if ret_bin_edges=False
    else, Returns (bin_edges, y_mean, y_err)
    """
    if (isinstance(bins, int)) and (nbins is None):
        
        nbins = bins
        
    elif (hasattr(bins, '__iter__')) and (binedges is None):
        
        binedges = bins
        
    if binedges is not None:
        
        nbins = len(binedges) - 1
        
    else:
        
        binedges = np.linspace(x.min(), (1 + np.spacing(2))*x.max(), nbins + 1)
        
    if weights is None:
        
        weights = np.ones(x.shape)

    weights = strip_units(weights)
    
    # Pre-factor for weighted STD:
    A = 1/(1 - (weights**2).sum())
    
    
    # Initialize
    y_mean = np.zeros(nbins)
    y_std = np.zeros(nbins)
    # Find the index bins for each data point
    ind = np.digitize(x, binedges) - 1
    # Ignore nans
    nan_ind = np.isnan(y)
    N = np.histogram(x, binedges)[0]
    
    # Initialize bin_centers (try to retain units)
    bin_centers = 0.0*binedges[1:]
    
    for i in range(nbins):
        
        #Indices to use
        mask = (ind==i) & (~nan_ind)
        # Set up the weighting
        w = weights[mask].copy()
        w /= w.sum()
        A = 1/(1 - (w**2).sum())
        #y_mean[i] = np.nanmean(y[mask])
        y_mean[i] = (w * y[mask]).sum()
        var = A*(w*(y[mask] - y_mean[i])**2).sum()
        y_std[i] = np.sqrt(var)
        #y_std[i] = np.std(y[use_ind])
        
        if weighted_bins:
            # Center of mass of x positions
            bin_centers[i] = (w*x[mask]).sum()
        
    y_mean = match_units(y_mean, y)[0]
    y_err = y_std/np.sqrt(N)
    y_err = match_units(y_err, y)[0]

    y_mean[N==0] = np.nan
    y_err[N==0] = np.nan
    
    if not weighted_bins:
        
        bin_centers = (binedges[0:-1] + binedges[1:])/2.0
        binedges = match_units(binedges, x)[0]
        bin_centers = match_units(bin_centers, x)[0]
        
    else:
        
        bin_centers[N==0] = np.nan
    
    if ret_bin_edges:
        
        return binedges, y_mean, y_err
        
    else:
        
        return bin_centers, y_mean, y_err
    
def heatmap(x, y, z, bins=10, plot=True, output=False):
    """
    Creates a pcolor heatmap for z evaluated at (x,y).  z is binned and
    averaged according to x and y.  x, y, and z should be 1-D arrays with the
    same length.
    
    IF bins = N, a pcolor plot of shape (N,N) is returned
    IF bins = (M,N) [a tuple], a pcolor plot of shape (M,N) is returned
    
    IF plot = True (default) a plot is created.
    
    *** RETURNS ***
    IF output = False, nothing is returned (default)
    
    IF output = True:
    
    Returns x_mesh, y_mesh, z_binned
    
    x_mesh, y_mesh are the meshgrid x,y edges z is evaluted in.  z_binned is
    the average of z for each bin.
    """
    
    N, x_binedges, y_binedges = np.histogram2d(x, y, bins = bins)
    x_ind = np.digitize(x, x_binedges) - 1
    y_ind = np.digitize(y, y_binedges) - 1
    nx_bins = len(x_binedges) - 1
    ny_bins = len(y_binedges) - 1
    z_binned = np.zeros([nx_bins, ny_bins])
    
    for i in range(nx_bins):
        
        for j in range(ny_bins):
            
            z_binned[i,j] = z[(x_ind==i) & (y_ind==j)].mean()
            
    x_mesh, y_mesh = np.meshgrid(x_binedges, y_binedges, indexing = 'ij')
    
    if plot:
        
        cmap = copy.copy(matplotlib.cm.jet)
        cmap.set_bad('w',1.)
        masked_z = np.ma.array(z_binned, mask=np.isnan(z_binned))
        plt.pcolormesh(x_mesh, y_mesh, masked_z, cmap = cmap)
        plt.colorbar()

    if output:
        
        return x_mesh, y_mesh, z_binned     
    

def configparser(fname,ftype='auto'):
    """
     --------------------------------------------------
        parameters = configparser(fname,ftype='auto')
        
    Tries to parse ChaNGa configuration files
    
    ftype can be 'auto', 'param', or 'director'.  If auto, config parser will
    try to determine the filetype on its own.
    
    returns:
        dictionary 'parameters'.  The keys are the names of the parameters and 
        the values are the values defined in the file fname
     --------------------------------------------------
     """
    types = np.array(['param','director'])
    ftype = ftype.lower()
    param = {}
    if ftype == 'auto':
        # Try using extension to determine file type
        a = fname.split('.')
        ftype = a[-1].lower()
    if np.sum(types == ftype) == 0:
        # Could not find file type
        print ('Could not determine config filetype...exiting')
        return param
        # Try to determine filetype
    # --------------------------------------------------
    # Parse param file
    # --------------------------------------------------
    if ftype == 'param':
        farray = np.genfromtxt(fname,delimiter='=',dtype=None)
        for n in range(len(farray)):
            param[farray[n,0].strip()] = str2num(farray[n,1].strip())
    # --------------------------------------------------
    # Parse director file
    # --------------------------------------------------
    elif ftype == 'director':
        f = open(fname,'r')
        f.seek(0)
        dummy = 0
        for line in f:
            a = line.strip().split()
            if len(a) == 1:
                # we're dealing with a flag
                param[a[0]] = str2num(a[0])
            elif len(a) > 1:
                param[a[0]] = str2num(a[1:])
            else:
                # This is an empty line
                pass
        f.close()
    # --------------------------------------------------
    # Throw warning, return 'param' as empty
    # --------------------------------------------------
    else:
        warnings.warn('Still cannot determine filetype.')
    return param

def configsave(param,filename,ftype='auto'):
    """
     --------------------------------------------------    
    Saves parameters defined by param (see configparser) to filename.
    Possible ftypes are 'director' and 'param'.  If set to auto, configsave
    tries to guess file type from the extension.
     --------------------------------------------------
     """
    f = open(filename,'w')
    types = np.array(['param','director'])
    ftype = ftype.lower()
    if ftype == 'auto':
        # Try to figure out filetype
        a = filename.split('.')
        ftype = a[-1].lower()
    if ftype == 'param':
        pars = sorted(param.iteritems())
        for n in range(len(pars)):
            f.write('{0:25s}= {1}\n'.format(pars[n][0],pars[n][1]))
    elif ftype == 'director':
        values = param.values()
        keys = param.keys()
        for n in range(len(keys)):
            outstr = keys[n]
            if outstr == values[n]:
                # We just have a flag
                pass
            elif isinstance(values[n],(float,int,str)):
                outstr = outstr + ' {0}'.format(values[n])
            else:
                outstr = outstr + ' ' + ' '.join(map(str,values[n]))
            f.write('{0}\n'.format(outstr))
    else:
        #no file type
        warnings.warn('no such filetype {0}\nCould not save'.format(ftype))
    f.close()
    
def extrap1d(x,y):
    """
    Calculates a linear interpolation of x and y and does a linear
    extrapolation for points outside of x and y.
    Uses scipy.interpolate.interp1d
    """
    # Ignore nans
    ind = (~np.isnan(x)) & (~np.isnan(y))
    x = x[ind]
    y = y[ind]
    # calculate interpolation
    yspline = interp.interp1d(x,y,kind='linear')
    
    def fcn(x0):
        
        if hasattr(x0,'__iter__'):
            
            mask1 = x0 < x.min()
            mask2 = x0 > x.max()
            out = np.zeros(len(x0))
            out[mask1] = y[0] +  (x0[mask1] - x[0])*(y[1]-y[0])/(x[1]-x[0])
            out[mask2] = y[-1] + (x0[mask2] - x[-1])*(y[-1] - y[-2])/(x[-1] - x[-2])
            mask3 = (~mask1) & (~mask2)
            out[mask3] = yspline(x0[mask3])
            
        else:
            
            if x0 < x.min():
                
                out = y[0] +  (x0 - x[0])*(y[1]-y[0])/(x[1]-x[0])
                
            elif x0 > x.max():
                
                out = y[-1] + (x0 - x[-1])*(y[-1] - y[-2])/(x[-1] - x[-2])
                
            else:
                
                out = yspline(x0)
                
            # Don't return an array with one element
            out = float(out)
        
        return out        
    
    return fcn
            
def smoothstep(x,degree=5,rescale=False):
    """
    Calculates a smooth step function y(x) evaluated at the data points x.
    x should be a numpy array or float.  
    
    y(x) is a polynomial of order 'degree' (default is 5).  degree must be an
    odd number between 3 and 25 (inclusive).  The higher the order, the 
    sharper the step is.
    
    y(x) is defined by:
        y(0) = 0
        y(1) = 1
        The first (degree - 1)/2 derivatives are 0 at y = 0,1
        
    *** ARGUMENTS ***
    
    * x * Points at which to evaluate the smoothstep
    
    * degree * Degree of the smooth step.  Must be odd number between 3 and 25
        default = 5
        
    * rescale *  Rescale x to be between 0 and 1.  Default = False.  If True,
        x MUST be an array (greater than length 1)
    
    
    *** RETURNS ***
    
    """
    # -----------------------------------------------------------
    # Load up the hermite spline (polynomial) coefficients
    # -----------------------------------------------------------
    fname = os.path.join(self_dir,'hermite_spline_coeffs.dat')
    f =open(fname,'r')
    
    coeffs_list = []
    order_list = []
    
    for line in f:
        
        l = line.strip().split(',')
        order_list.append(int(l[0]))
        
        for n in range(len(l)):
            
            l[n] = float(l[n].strip())
            
        coeffs_list.append(np.array(l[1:],dtype='float'))
    
    order = np.array(order_list)
    coeffs = coeffs_list[(order==degree).argmax()]
    # -----------------------------------------------------------
    # Calculate the smooth step function y(x)
    # -----------------------------------------------------------
    n_coeffs = len(coeffs)
    
    if rescale:
        
        try:
            x = (x - x.min())/(x.max() - x.min())
        except:
            raise RuntimeError,'Could not rescale x.  Make sure x is an array'
    
    if isinstance(x, (int, long, float, complex)):
        
        # x is a number, handle accordingly
        y = 0.0
        
        if (x > 0) & (x < 1):
            # If 0<x<1, calculate the smooth step
            for n in range(n_coeffs):
                
                y += coeffs[n] * x**(degree - n)
                
        elif x <= 0:
            
            y = 0.0
            
        else:
            
            y = 1.0
        
    else:
        
        # Assume x is a numpy array
        y = np.zeros(x.shape)
        ind = (x > 0) & (x < 1)
        
        for n in range(n_coeffs):
            
            y[ind] += coeffs[n] * x[ind]**(degree-n)
            
        y[x >= 1] = 1
    
    return y
    

def str2num(string):
    """
     --------------------------------------------------
     Tries to see if 'string' is a number
     
     If 'string' is a string, returns:
       int(string) for integers
       float(string) for floats
       'string' otherwise
    
     If 'string' is a float or an integer, returns:
       string
    
     If none of the above, treats it like a list or tuple
     and returns for each entry of 'string' a float,int,
     or str as required.  Returns as a list
     --------------------------------------------------
     """
    if isinstance(string,int):
        output = string
    elif isinstance(string,float):
        output = string
    elif not isinstance(string,str):
        output = []
        for a in string:
            try:
                output.append(int(a))
            except:
                try:
                    output.append(float(a))
                except:
                    output.append(a)
        if len(output) == 1:
            output = output[0]
    else:
        output = string
        try:
            output = int(string)
        except:
            try:
                output = float(string)
            except:
                pass
    return output

def loadhalos(fname=''):
    """
     Load halo (.grp) file generated from fof
     Should be an ascii list of numbers, where the first row contains the
     total number of particles (gas+star+dark) and the remaining rows define
     which halo each particle belongs to
     """
    if fname == '':
        # Empty filename
        pass
    grp = np.loadtxt(fname,dtype=np.uint16)
    grp = grp[1:]   # (ignore the number of particles)
    
    return grp

def fof(fFilter,saveDir='',minMembers=8,linklen=0.01):
    """
     --------------------------------------------------
     A simple script that allows you to loop through calls to fof
     for many files in one directory
     --------------------------------------------------
     """
    flist = np.sort(glob.glob(fFilter))
    nfiles = len(flist)
    if (saveDir != '') and nfiles > 0:
        if ~os.path.isdir(saveDir):
            os.makedirs(saveDir)
    
    for n in range(nfiles):
        fname = flist[n]
        outname = os.path.join(saveDir,fname)
        os.system('totipnat < {0} | fof -g -m {1} -e {2} -o {3}'.format(fname,minMembers,linklen,outname))
        
def pbverbosity(cmd=None):
    """
    Changes and returns pynbody verbosity.  Works for different versions
    of pynbody.
    
    **ARGUMENTS**
    
    cmd
        -If None (default) current verbosity level is returned, nothing is done
        -If 'off', pynbody is silenced
        -If 'on', pynbody verbosity is set on
        -If something else, cmd is assumed to be a verbosity level
        
    **RETURNS**
    
    current_verbosity
        pynbody verbosity level before any changes were made
        
    **EXAMPLES**
    
    *Toggle pynbody verbosity*
    
        current_verbosity = pbverbosity('off')
        ...
        do stuff
        ...
        pbverbosity(current_verbosity)
    """
    
    # -----------------------------
    # Get current verbosity level
    # -----------------------------
    if hasattr(pb, 'logger'):
        # As of v0.30, pynbody uses python's logging to handle verbosity
        logger = True
        current_verbosity = pb.logger.getEffectiveLevel()
        pb.logger.setLevel(logging.ERROR)
        
    else:
        
        # For pynbody version < 0.3, verbosity is handled in the config
        logger = False
        current_verbosity = pb.config['verbose']
        
    # -----------------------------
    # Change verbosity
    # -----------------------------
    if cmd is None:
        # Don't change verbosity.  just return the current verbosity
        pass
        
    elif cmd == 'off':
        # Toggle verbosity off
        if logger:
            
            pb.logger.setLevel(logging.ERROR)
            
        else:
            
            pb.config['verbose'] = False
        
    elif cmd == 'on':
        # Toggle verbosity on
        if logger:
            
            pb.logger.setLevel(logging.DEBUG)
        
        else:
            
            pb.config['verbose'] = True
        
    else:
        # Set verbosity to the verbosity level specified by cmd
        if logger:
            
            pb.logger.setLevel(cmd)
            
        else:
            
            pb.config['verbose'] = cmd
        
    # Return the verbosity level before any changes were made
    return current_verbosity
