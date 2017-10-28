# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 15:15:23 2015

@author: ibackus
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import LogNorm
import numpy as np
import copy
import os
from pynbody.plot.sph import image

from diskpy.disk import powerspectrum_t
from diskpy.utils import get_units

def savefig(savename, exts=('pdf', 'png'), verbose=True, ext_kwargs=None, **kwargs):
    """
    A wrapper around matplotlib.pyplot.savefig to save a figure as multiple
    file types (.png, .pdf, etc...)
    
    Parameters
    ----------
    savename : str
        Filename WITHOUT extension
    exts : iterable (list, tuple, etc...)
        File extensions to save as.  Format is inferred from the extension.
        default is ('pdf', 'png')
    verbose : bool
        verbosity flag
    ext_kwargs : iterable, optional
        A list (or tuple) of dicts, each specifying keyword arguments
        to pass to matplotlib.pyplot.savefig for each format.  Should be
        length of exts.  Overrides **kwargs
    **kwargs : extra args
        Extra args are passed to matplotlib.pyplot.savefig
    """
    if not hasattr(exts, '__iter__') or len(exts) < 1:
        raise ValueError, 'Cannot iterate over provided extensions'
    if ext_kwargs is None:
        ext_kwargs = [{}] * len(exts)
    
    for ext, extra_args in zip(exts, ext_kwargs):
        fullname = savename + '.' + ext
        save_args = copy.deepcopy(kwargs)
        save_args.update(extra_args)
        plt.savefig(fullname, **save_args)
        if verbose:
            print "saved:", fullname

def colorCodeArray(x1, x2, cmap='Reds'):
    """
    Get a color-coded array from x1 (brightess) and x2 (color).  x1 and x2 
    are assumed to be normalized (0 to 1).  x1 defines the brightness, x2
    defines the color.  Returns an RGB array
    """
    # Get color image for the color-coded qty
    cm = plt.get_cmap(cmap)
    colors = cm(x2)
    # Drop alpha to retrieve RGB
    colors = colors[:, :, 0:-1]
    # get hsv (Hue-Saturation-Value)
    hsv = mpl.colors.rgb_to_hsv(colors)
    # Set value (brightness) to the first gray-scale image
    hsv[:, :, 2] = x1
    
    rgb = mpl.colors.hsv_to_rgb(hsv)
    return rgb

def normalizeArray(x, log=True, vmin=None, vmax=None):
    """
    Normalize an array.  Values below/above the thresholds of vmin/vmax are
    set to 0/1
    """
    if log:
        Norm = mpl.colors.LogNorm
    else:
        Norm = mpl.colors.Normalize
    norm = Norm(vmin=vmin, vmax=vmax)
    x = norm(x)
    x[x < 0] = 0
    x[x > 1] = 1
    
    return x

def colorCodedImage(sim, qtys, width = 60, resolution = 500, z_camera = None,
                    av_zs = (False, False), cmap = 'Reds', vmins = (None, None), 
                    vmaxs = (None, None), logs = (True, True)):
    """
    Generates an image as with pynbody.plot.sph.image with brightness set by
    one quantity and color set by a second quantity.
    
    Parameters
    ----------
    sim : SimSnap
        Snapshot or sub-snap to render
    qtys : list or tuple
        A tuple/list of 2 quantities (strings) to render.  The quantity are
        SimSnap keys, e.g. 'rho', 'temp', 'mass', etc.  The first quantity
        sets the brightness of the image, the second sets the color
    width : number
        Width of the plot to render
    resolution : int
        Number of pixels to render along each axis
    z_camera : number
        z-camera location (see pynbody.plot.sph.image)
    av_zs : tuple
        Average along line-of-sight (see pynbody.plot.sph.image).  If True,
        volume averating is used.  If a string averaging is over that qty.
        E.g., mass averaging can be acheived by using av_zs = ('rho', 'rho')
    cmap : str or cmap
        colormap to use for color-coding
    vmins, vmaxs : list/tuples
        vmin, vmax for the quantities (one for each quantity)
    logs : list/tuple
        List or tuple of booleans.  Defines whether color should use a log
        scale.
    
    Returns
    -------
    im : array
        The rendered RGB image used
    """
    cmaps = ('gray', cmap)
    ims = []
    # Generate the images
    for i in range(2):
        im = image(sim, qtys[i], width, resolution, av_z=av_zs[i], 
                           z_camera=z_camera, log=logs[i], vmin=vmins[i],
                           vmax=vmaxs[i], cmap=cmaps[i], ret_im=True)
        cb = plt.gcf().axes[-1]
        ims.append(im)
    
    # Grayscale image arrays
    grayims = [im.get_array() for im in ims]
    for i in range(2):
        grayims[i] = normalizeArray(grayims[i], logs[i], vmins[i], vmaxs[i])
    # Convert to rgb
    rgb = colorCodeArray(grayims[0], grayims[1], cmap=cmap)
    # Plot
    #rgb = mpl.colors.hsv_to_rgb(hsv)
    top = 0.5 * width
    plt.gca().imshow(rgb, extent = (-top, top, -top, top))
    # Add colorbar
    #ax = im.get_axes()
    cb = plt.colorbar(ims[1], orientation='vertical', label=qtys[1])
    cb = plt.colorbar(ims[0], orientation='horizontal', label=qtys[0])
    plt.ylabel('y ($' + sim['y'].units.latex() + '$)')
    plt.tight_layout()
    return rgb



def gridplot(nrows, ncols=1, square=False):
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
            
            if square:
                
                ax[i,j].set(adjustable='box-forced', aspect='equal')
                
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

def lineanimate(x, y, *args, **kwargs):
    """
    Makes an animation of a line vs time and saves as .mp4 file.  Uses 
    matplotlib.animation.
    
    Parameters
    ----------
    
    x : 1D array like
        x points (same at every time step)
    y : 2D array like
        y points (change at each time step).  y[i] is plotted at frame i
    *args
        Optional additional arguments passed to matplotlib.plot
    **kwargs
        Optional keyword args.  The Some are used in movieplot, all others
        are passed to matplotlib
        
    Keyword-arguments
    -----------------
    
    fname : str
        Filename to save movie to
    xlim : list, arraylike
        x limits of the plot
    ylim : list, arraylike
        y limits of the plot
    fps : int
        Frames per second
    ax : matplotlib axes instance
        axes instace to plot to.  Can be used to set up the plot beforehand 
        (log scale, plot labels, title, etc)
    verbose : bool
        verbosity (true or false)
    
    Returns
    -------
    
    None
    """
    
    # Parse kwargs and set defaults
    keys = ['fname', 'xlim', 'ylim', 'fps', 'ax','verbose']
    filename = kwargs.get('fname', 'movie.mp4')
    fps = kwargs.get('fps', 25)
    verbose = kwargs.get('verbose', False)
    
    xlim = kwargs.get('xlim')
    ylim = kwargs.get('ylim')
    ax = kwargs.get('ax')
    
    if ax is None:
        
        ax = plt.gca()
        
    if xlim is None:
        
        xlim = [x.min(), x.max()]
        
    if ylim is None:
        
        ylim = [y.min(), y.max()]
        
    # Ignore keys not meant for plot calls
    for key in keys:
        if key in kwargs:
            del kwargs[key]
            
    # Set up axes    
    fig = ax.figure
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    # Draw an empty line to update
    line, = ax.plot([], [], *args, **kwargs)
    
    # Initialization function for matplotlib animate.  Makes line empty
    def init():
        line.set_data([],[])
        return line,
    
    # Frame drawing function
    def animate(i):
        
        line.set_data(x, y[i])
        plt.title(i)
        
        if verbose:
            
            print i
            
        return line, 
    
    iMax = len(y)
    # Make animation
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=iMax,\
    interval=10, blit=True)
    # Render frames and save
    if os.path.exists(filename):
        
        os.remove(filename)
        
    anim.save(filename, fps=fps)
    print 'Movie saved to ' + filename
    
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
        
def waterfall(power, t=None, m=None, log=True, normalize=True, 
              cmap='cubehelix', vmin=None, vmax=None, colorbar=True):
    """
    Generates a waterfall plot for power spectrum vs time for a simulation
    
    Parameters
    ----------
    
    power : array, SimArray
        2D shape (nt, nm) array of power spectrum
    t : array, SimArray
        (optional) time points.  1D array
    m : array, SimArray
        (optional) Fourier mode numbers.  1D array
    log: bool
        logscale the colors
    normalize : bool
        Normalize by the DC component at each time step
    cmap : str or colormap
        colormap to use
    vmin, vmax : float
        limits
    colorbar : bool
        Display colorbar
    
    Examples
    --------
    
    >>> flist = diskpy.pychanga.get_fnames('snapshot')
    >>> m, power = diskpy.disk.powerspectrum_t(flist, spacing='log')
    >>> waterfall(power)
    
    """
    # Initialize m and t
    nt, nm = power.shape
    
    if m is None:
        
        m = np.arange(nm)
        
    if t is None:
        
        t = np.arange(nt)
        
    mMesh, tMesh = np.meshgrid(m, t)
    
    # Normalize
    if normalize:
        
        # Power in DC component
        p0 = power[:, 0, None]
        # Make a 2D array (same shape as power)
        p0 = np.dot(p0, np.ones([1, nm]))
        # normalize
        power = power/p0
    
    # Set up log-scale for plot
    if log:
        
        norm = LogNorm(vmin, vmax)
        
    else:
        
        norm = None
        
    plt.pcolormesh(tMesh[:, 1:], mMesh[:, 1:], power[:, 1:], norm=norm, \
    cmap = cmap)
    
    plt.ylim(m[1], m[-1])
    plt.ylabel('m')
    
    tUnits = get_units(t)
    if tUnits == 1:
        
        plt.xlabel('time step')
        
    else:
        
        plt.xlabel(' time $(' + tUnits.latex() + ')$')
    
    if colorbar:
        
        cbar = plt.colorbar()
        if normalize:
            
            text = r'$A_m/A_0$'
            
        else:
            
            units = get_units(power)
            text = r'$A_m (' + units.latex() + r') $'
            
        cbar.set_label(text)

