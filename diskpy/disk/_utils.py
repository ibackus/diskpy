# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 15:43:08 2015

@author: ibackus
"""
import pynbody as pb
import numpy as np

def centerdisk(snapshot):
    """Centers a disk on the center of mass puts it into rest frame
    """
    cm = pb.analysis.halo.center_of_mass(snapshot)
    vel = pb.analysis.halo.center_of_mass_velocity(snapshot)

    snapshot['pos'] -= cm
    snapshot['vel'] -= vel

    return
    
def snapshot_defaults(snapshot):
    """Applies various defaults to tipsy snapshots of protoplanetary disk 
    simulations. These include:

        | -Sets nice units        
        | -Calculates particle smoothing lengths using mass and rho 
        |   (if available)        
        | -Centers on snapshot center-of-mass
        | -Puts in rest frame

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

    # Center the disk
    centerdisk(snapshot)

    return