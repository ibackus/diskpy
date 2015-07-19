# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 15:42:56 2015

@author: ibackus
"""
import warnings
import pynbody as pb
SimArray = pb.array.SimArray
import numpy as np

from diskpy.utils import match_units
from diskpy.pdmath import binned_mean

def Q(snapshot, molecular_mass = 2.0, bins=100, use_velocity=False, \
use_omega=True):
    """Calculates the Toomre Q as a function of r, assuming radial temperature
    profile and kappa ~= omega
    
    Parameters
    ----------
    
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

    Returns
    -------
    
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
            p = pb.analysis.profile.Profile(snapshot, bins=r_edges)
            omega = p['omega']

        kappa_calc = omega

    else:

        if use_velocity:
            # Calculate directly from particle's velocities
            kappa_calc, dummy = kappa(snapshot, r_edges)

        else:
            # Estimate, from forces, using pynbody
            p = pb.analysis.profile.Profile(snapshot, bins=r_edges)
            kappa_calc = p['kappa']

    return (kappa_calc*c_s/(np.pi*G*sig)).in_units('1'), r_edges

def Qeff(snapshot, molecular_mass = 2.0, bins=100, use_velocity=False, \
use_omega=True, alpha=0.18, beta=2.2):
    """Estimates the effective Toomre Q as a function of r, defined as:
    
    .. math:: Q_{eff} = \\beta Q (h/R)^{\\alpha}
        
    See Q and h for the estimates of Q and h
    
    Parameters
    ----------
    
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
    alpha : float
        Powerlaw for height dependence
    beta : float
        Normalization such that disks fragment for Qeff = 1
    
    Returns
    -------
    
    Qeff : array
        Effective Toomre Q as a function of r
    r_edges : array
        Radial bin edges
    
    """
    Qcalc, r_edges = Q(snapshot, molecular_mass, bins, use_velocity, use_omega)
    dummy, h = height(snapshot, r_edges, center_on_star=False)
    r = (r_edges[1:] + r_edges[0:-1])/2.
    Qeff = beta * ((h/r).in_units('1'))**alpha

    return Qeff, r_edges
    
def kappa(f, bins=100):
    """Estimate the epicyclic frequency from velocity
    
    Parameters
    ----------

    f : TipsySnap
        `f` is a Simulation snapshot
    bins : int or array-like
        Either the number of bins to use or the bin edges

    Returns
    -------

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

            warnings.warn('Bins not uniformly spaced')

    r = f.g['rxy']
    v = f.g['vt']

    r_edges, v_mean, dummy = binned_mean(r, v, bins=bins, ret_bin_edges=True)
    dummy, rv_mean, dummy2 = binned_mean(r, r*v, bins=r_edges)
    r_cent = (r_edges[1:] + r_edges[0:-1])/2
    dr = r_edges[[1]] - r_edges[[0]]
    drv_dr = np.gradient(rv_mean, dr)

    kappa = np.sqrt(2*v_mean*drv_dr)/r_cent

    return kappa, r_edges
    
def height(snapshot, bins=100, center_on_star=True):
    """
    Calculates the characteristic height (h) of a flared disk as a function
    of cylindrical radius (r).

    Parameters
    ----------

    snapshot : TipsySnap
        Simulation snapshot for a flared disk
    bins : int or array_like
        Specifies the bins to use.  If int, specifies the number of bins.  If
        array_like, specifies the bin edges
    center_on_star : bool
        If true (DEFAULT), cylindrical r is calculated relative to the star

    Returns
    -------
    
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
    """Calculates surface density vs r (relative to the center of mass)
    
    Parameters
    ----------
    
    snapshot : tipsy snapshot
    bins : int, list, array...
        Either the number of bins to use or the binedges to use
    
    Returns
    -------
    
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