# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 15:15:23 2015

@author: ibackus
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import copy

def gridplot(nrows, ncols=1):
    """
    Creates a grid of tightly-packed subplots and returns them as a numpy array,
    shape (nrows,ncols).  If nrows=ncols=1, a single subplot is made.
    Currently not fully implemented for a 2D array
    """
    # Just make a single subplot
    if (nrows == 1) & (ncols == 1):

        return plt.subplot(1,1,1)

    # Create a grid
    grid = mpl.gridspec.GridSpec(nrows,ncols)
    grid.update(wspace=0., hspace=0.)

    # Initialize subplots
    ax = np.zeros((nrows,ncols), dtype=object)
    counter = 0
    for i in range(nrows):
        for j in range(ncols):
            if i > 0:
                sharex = ax[0,j]
            else:
                sharex = None
            if j > 0:
                sharey = ax[i,0]
            else:
                sharey = None

            ax[i,j] = plt.subplot(grid[counter], sharex = sharex, sharey = sharey)
            counter += 1

    # Remove ticklabels inbetween plots
    for i in range(nrows-1):
        for j in range(ncols):
            plt.setp(ax[i,j].get_xticklabels(), visible=False)
    for i in range(nrows):
        for j in range(1,ncols):
            plt.setp(ax[i,j].get_yticklabels(), visible=False)

    # If this is a 1-D grid, flatten ax
    if (ncols == 1) or (nrows == 1):

        ax = ax.flatten()

    return ax
    
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

        cmap = copy.copy(mpl.cm.jet)
        cmap.set_bad('w',1.)
        masked_z = np.ma.array(z_binned, mask=np.isnan(z_binned))
        plt.pcolormesh(x_mesh, y_mesh, masked_z, cmap = cmap)
        plt.colorbar()

    if output:

        return x_mesh, y_mesh, z_binned