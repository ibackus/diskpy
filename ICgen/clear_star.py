# -*- coding: utf-8 -*-
"""
Deletes particles within a given radius of the star and recenters the snapshot
on the center of mass of the system.


Created on Tue Feb 25 10:27:39 2014

@author: ibackus
"""
import pynbody
SimArray = pynbody.array.SimArray
import numpy as np

def clstar(f,r_cut=None,r_max=None, center_on_star=True):
    """
    Deletes particles within a given radius of the star and recenters the snapshot
    on the center of mass of the system.
    
    *** ARGUMENTS ***

    * f *       Either a filename pointing to a snapshot to be loaded -OR- 
        a pynbody SimSnap    
    
    * r_cut *   Remove particles at r < r_cut from the star.  If None, defaults
        to 0.02 AU
    
    * r_max *   Remove particles at r > r_max.  If None, r_max is infinity
    
    * center_on_star *  True or False.  If True, the simulation is re-centered
        on the star.
    
    *** RETURNS ***
    
    Returns a re-centered SimSnap with particles removed
    """
    if r_cut is None:
        # Default radial cut-off for particles.  IE, delete all particles
        # at r < r_cut
        r_cut = SimArray(0.02,'au')
    if isinstance(f,basestring):
        # Assume f is a file name for a snapshot
        fName = f
        f = pynbody.load(fName)
    else:
        # Assume f is a pynbody SimSnap
        fName = f.filename
    # Center around star
    starPos = f.star['pos'][0].copy()
    for n in range(3):
        f['pos'][:,n] -= starPos[[n]]
    # Remove particles at r < r_cut (except the star)
    r = f.gas['rxy']
    use_ind = np.ones(len(f),dtype='bool')
    if r_max is None:
        use_ind[0:-1] = (r >= r_cut)
    else:
        use_ind[0:-1] = ((r >= r_cut) & (r < r_max))
    
    n_cut = (~use_ind).sum()
    
    if n_cut > 0:
        
        # Transfer momentum to the star
        m_g = f.g['mass'][[0]]
        m_g_tot = m_g * (~use_ind).sum()
        p_g = (m_g * f['vel'][~use_ind]).sum(0)
        p_s = f.s['mass']*f.s['vel']
        v_s = (p_s + p_g)/(f.s['mass'] + m_g_tot)
        f.s['vel'] = v_s.in_units(f.s['vel'].units)
        f.s['mass'] += m_g_tot
        
        # Delete particles
        f = f[use_ind]
        print 'Deleted {0} particles'.format(n_cut)

    # If we're not centering on the star, add it's initial position back in
    if not center_on_star:
        
        for n in range(3):
            f['pos'][:,n] += starPos[[n]]
            
    return f