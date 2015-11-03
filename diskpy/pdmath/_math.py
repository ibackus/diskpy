# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 14:26:07 2015

@author: ibackus
"""

import os
import numpy as np
import scipy.interpolate as interp
interp1d = interp.interp1d
import pynbody as pb
SimArray = pb.array.SimArray

from diskpy.utils import strip_units, match_units

def _loadcoeffs(fname):
    """
    Loads hermite polynomial coefficients stored in fname and returns them
    as a dictionary, where the keys are the degree of the polynomial and the
    values are the coefficients
    """
    
    # Load up the hermite spline (polynomial) coefficients
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
    coeffs = {}
    for i in range(len(order)):
        
        coeffs[order[i]] = coeffs_list[i]
        
    return coeffs

# Load hermite coefficients
_dir = os.path.dirname(os.path.realpath(__file__))
_coeffsfile = os.path.join(_dir, 'hermite_spline_coeffs.dat')
hermite_coeffs = _loadcoeffs(_coeffsfile)

def dA(redges, thetaedges):
    """
    Calculates the area of bins in cylindrical coordinates as
    
    .. math:: dA = r(\\Delta r) (\\Delta \\theta)
    
    on a grid of r, theta values
    
    Parameters
    ----------
    
    redges, thetaedges : array like
        1D arrays of the binedges in r, theta.
    
    Returns
    -------
    
    dA : array
        2D array of dA values over r, theta.  dA[i,j] corresponds to r[i]
        theta[j]
    """
    
    dr = redges[1:] - redges[0:-1]
    dtheta = thetaedges[1:] - thetaedges[0:-1]
    r = (redges[1:] + redges[0:-1])/2.
    rdr = r * dr
    return np.dot(rdr[:, None], dtheta[None, :])
    
def setupbins(x, bins=10):
    """
    Sets up bins for data x (similar to numpy.histogram).  If bins is an 
    integer, the min/max are set to include all data
    
    Parameters
    ----------
    
    x : arraylike
        Data to be binned
    bins : int or array-like
        Number of bins or binedges to use
    
    Returns
    -------
    
    binedges : array
        Bin edges
    """
    
    # If bins is not iterable, it is the number of bins
    if not hasattr(bins, '__iter__'):
        
        dtype = x.dtype
        if isinstance(dtype, int):
            
            xmin = x.min() - 0.5
            xmax = x.max() + 0.5
            
        else:
            
            eps = np.finfo(x.dtype).eps
            xmin = x.min() * (1 - 2*eps)
            xmax = x.max() * (1 + 2*eps)
            
        binedges = np.linspace(xmin, xmax, bins + 1)
        
    else:
        # Bins is already binedges.  do nothing
        binedges = bins
    
    return binedges

def resolvedbins(x, y, minbins=200, ftol=None):
    """
    Generates a sub-view of x that allows y to be resolved up to ftol.  Since
    f will only be sampled at x, x should be very (overly) high resolution.
    
    Parameters
    ----------
    
    x : array or SimArray
        initial bin positions to try.  should be many more than required
    y : array or SimArray
        y(x), same shape as x.  the function values to resolve
    minbins : int
        Minimum number of bins to return.  IF the ftol is met but there are 
        not enough bins, ftol is decreased
    ftol : SimArray, float
        Function tolerance.  Maximum difference between the interpolated, 
        subsampled y(x) and the full-ly resolved version
    
    Returns
    -------
    
    x2 : array or SimArray
        A subview of x that provides the required resolution
    
    """
    
    # Initialize
    if ftol is None:
        ftol = y.max() * 1e-5
    nx = len(x)
    # Begin with bin edges at both ends
    binind = [0, nx-1]
    
    maxiter = int(1e5)
    
    # Main loop
    for i in range(maxiter):
        
        # Linear spline interpolation
        yspl = interp1d(x[binind], y[binind], kind='linear')
        # Difference between interpolated and actual values (on the grid)
        diff = abs(y - yspl(x))
        
        # Find if any new bins need to be made
        nbins = len(binind) - 1    
        newbins = []            
        for i in range(nbins):
            
            # Get the indices of the current bin's edges
            iLo = binind[i]
            iHi = binind[i+1]
            
            # If they are separated by more than one, we may need another bin
            if (iHi-iLo) > 1:
                
                # Find maximum error in the bin
                err = diff[iLo+1:iHi].max()
                
                # If the error  exceeds the tolerance, add a new bin
                if err > ftol:
                    
                    newbin = int((iLo + iHi)/2)
                    newbins.append(newbin)
        
        # Check if there are new bins
        if len(newbins) == 0:
            
            # No new bins - check if we have enough
            if nbins < minbins:
                
                # Decrease ftol if we don't have enough bins
                ftol /= 1.5
                
            else:
                
                # we have enough bins, break
                break
        
        else:
            # There are new bins, keep going
            binind.extend(newbins)
            binind.sort()
    
    binind = np.array(binind)
    return x[binind]
    
def bin2dsum(x, y, z, xbins=10, ybins=10):
    """
    Bins x,y using bin2d and sums z in those bins
    
    Parameters
    ----------
    
    x, y, z: array-like
        x, y, and z values.  Bins are in x-y, z values are summed
    xbins, ybins: int or array-like
        (see bin2d) Number of bins or bin edges
    
    Returns
    -------
    
    zbinned : array-like
        2D array of z-values, summed in the bins.  zbinned[i,j] gives the value
        of z summed in xbin[i], ybin[j]
    xedges, yedges : array
        1D arrays of binedges in x,y
    """
    
    ind, xedges, yedges = bin2d(x, y, xbins, ybins)
    xind = ind[0]
    yind = ind[1]
    nx = len(xedges) - 1
    ny = len(yedges) - 1
    
    zbinned = np.zeros([nx, ny])
    
    if pb.units.has_units(z):
        
        zbinned = SimArray(zbinned, z.units)
    
    for i in range(nx):
        
        xmask = (xind == i)
        for j in range(ny):
            
            mask = xmask & (yind == j)
            zbinned[i, j] = z[mask].sum()
            
    return zbinned, xedges, yedges
            

def bin2d(x, y, xbins=10, ybins=10):
    """
    2-dimensional binning of x, y
    Works as a 2-D extension of numpy.digitize but also automatically sets-up
    binedges
    
    Parameters
    ----------
    
    x, y : array-like
        x, y values to bin according to
    xbins, ybins : int OR list/array like
        Either the number of bins or the binedges to use
    
    Returns
    -------
    
    ind : list of arrays
        The x,y bin indices each entry belongs to.  ind[0][i] gives the x-bin
        of the ith entry.  ind[1][i] gives the y-bin of the ith entry
    xedges, yedges: arrays
        Bin edges used
    """
    
    xedges = setupbins(x, xbins)
    yedges = setupbins(y, ybins)
    
    xind = np.digitize(x, xedges) - 1
    yind = np.digitize(y, yedges) - 1
    ind =[xind, yind]
    return ind, xedges, yedges

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
    
def meshinterp(xedges, y, z, kind='linear', bounds_error=False, fill_value=0,
               assume_sorted=True):
    """
    Generates a 2D interpolating function for z defined on a non-uniform mesh
    Handles units
    
    Parameters
    ----------
    
    xedges : array
        1D array defining the x bin edges, monotonically increasing
    y : array
        2D array defining y values.  shape (nx, ny), where nx is the number
        of xedges and ny is the number of y-points at each x-bin
        So, y[i, :] are the monotonically increasing y values at xedges[i]
    z : array
        2D array of z(x,y).  shape (nx, ny) = y.shape
    kind : str
        (optional) Sets the kind of interpolation to perform
        [see scipy.interpolate.interp1d]
    bounds_error : bool
        (optional) Flag to raise error if values outside of y are called
        [see scipy.interpolate.interp1d]
    fill_value : float
        (optional) Sets the value to fill with if bounds_error = True
        [see scipy.interpolate.interp1d]
    assume_sorted : bool
        [see scipy.interpolate.interp1d]
    
    Returns
    -------
    
    meshspline(x, y): callable interpolation function
        Function which can be called on x, y pairs to give the interpolated
        value of z.  Values outside of the range of y are set to fill_value.
        x values outside the range of xedges are set to the boundary of xedges
    """
    # Check shapes
    
    if z.shape != y.shape:
        
        raise ValueError, 'y and z must have same shape'
        
    if len(xedges) != len(y):
        
        raise ValueError, 'x and y must have same len'
        
    # Handle units
    pos = [xedges, y, z]
    units = []
    for a in pos:
        
        if pb.units.has_units(a):
            
            units.append(a.units)
            
        else:
            
            units.append(None)
    
    xedges, y, z = strip_units(pos)
    # Setup bin information
    binsize = xedges[1:] - xedges[0:-1]
    xmin = xedges[0]
    xmax = xedges[-1]
    nbins = len(xedges) - 1
    
    # set up spliness
    splines = []
    
    for i in range(nbins+1):
        
        # perform interpolation to make spline
        splines.append(interp1d(y[i], z[i], kind=kind, \
        bounds_error=bounds_error, fill_value=fill_value, \
        assume_sorted=assume_sorted))
        
    # Define the callable interplation function to return
    def meshspline(x1, y1):
        """
        Callable interpolation function, interoplates the value of z at
        points (x1, y1)
        
        Parameters
        ----------
        
        x1, y1 : array
            x and y points to evaluate z at.  Must be the same shape.  ie,
            x1[i], y1[i] define a point (x, y).
            If @x1 or @y1 have no units, they are assumed to have the units of
            the nodes used to make the interpolator.  Otherwise they are
            converted to the proper units
            
        Returns
        -------
        
        z(x1, y1) : array
            z evaluated at @x1, @y1
        """
        # Handle units
        x1 = strip_units(match_units(x1, units[0])[0])
        y1 = strip_units(match_units(y1, units[1])[0])
        
        # Setup x and y points to estimate z at
        x1 = np.asarray(x1).copy()
        y1 = np.asarray(y1)
        
        if len(x1.shape) < 1:
            
            x1 = x1[None]
            
        if len(y1.shape) < 1:
            
            y1 = y1[None]
            
        # Flatten arrays
        shape = x1.shape
        nElements = np.prod(shape)
        x1 = np.reshape(x1, [nElements])
        y1 = np.reshape(y1, [nElements])
            
        # Deal with xs outside of boundaries
        x1[x1 < xmin] = xmin
        x1[x1 > xmax] = xmax
        # Find bin indices
        ind = np.digitize(x1, xedges) - 1
        ind[ind < 0] = 0
        ind[ind > nbins - 1] = nbins - 1
        
        # Get bin info for every point
        xlo = xedges[ind]
        xhi = xedges[ind + 1]
        dx = binsize[ind]
        
        # Get weights for bins (distance from bin edges)
        wlo = (xhi - x1)/dx
        whi = (x1 - xlo)/dx
        
        # Get function values at left and right xedges
        flo = np.zeros(x1.shape)
        fhi = np.zeros(x1.shape)
        
        for i in range(nbins):
            
            # Select everything in bin i
            mask = (ind == i)
            
            if np.any(mask):
                
                # Retrieve function values
                flo[mask] = splines[i](y1[mask])
                fhi[mask] = splines[i+1](y1[mask])
        
        # Take a weighted average of the function values at left and right
        # bin edges
        fout = wlo*flo + whi*fhi
        
        # Unflatten fout:
        fout = np.reshape(fout, shape)
        
        return SimArray(fout, units[2])
    
    return meshspline
    
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
    coeffs = hermite_coeffs[degree]
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