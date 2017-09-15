#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This implements a simple routine for writing/reading tables of arrays/simarrays
in a human readable ascii form in a way that preserves dtypes and units.

See write_table() and read_table()


Created on Thu Sep 14 11:45:35 2017

@author: ibackus
"""
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray

_COLPAD = 2

def get_units(x):
    """
    Gets the units of x.  Returns None if x has no units
    """
    try:
        return x.units
    except AttributeError:
        return None
    
def _line_slicer(colwidth):
    """
    Makes a function to slice a line string at the column points defined
    by colwidth
    
    colwidth is list-like, defining width of columns.  The standard column
    padding is assumed
    """
    i0 = 0
    iStart = []
    iEnd = []
    for width in colwidth:
        i1 = i0 + width
        iStart.append(i0)
        iEnd.append(i1)
        i0 = i1
    slice_ind = zip(iStart, iEnd)
    def slice_line(line, strip=False):
        """
        split line into different columns
        
        if strip == True, all the columns will have .strip() applied
        """
        line = line.strip('\n')
        out = [line[i0:i1] for i0, i1 in slice_ind]
        if strip:
            out = [x.strip() for x in out]
        return out
    return slice_line

def _line_length(colwidth):
    
    length = np.sum(colwidth)
    return length
    
def _read_header(f):
    """
    Reads a sim array table header from an open file buffer.
    NOTE: this reads starting at the current file pointer location.
    """
    import json
    header = [f.readline() for _ in range(6)]
    # First character in the header should be a comment, replace it with a 
    # space
    header = [' ' + line[1:] for line in header]
    # Parse header
    colwidth, names, dtypes, units, ind, hline = header
    colwidth = json.loads(colwidth)['colwidth']
    slice_line = _line_slicer(colwidth)
    names = slice_line(names, strip=True)
    # Names are formated as: '   "array name"   '
    names = [name.strip()[1:-1] for name in names]
    dtypes = slice_line(dtypes, strip=True)
    units = slice_line(units, strip=True)
    for i, unit in enumerate(units):
        if unit in ('None', 'none', 'NoUnit()', 'NoUnit'):
            units[i] = None
    ind = slice_line(ind, strip=True)
    ind = np.array([json.loads(i) for i in ind], dtype=int)
    # Read the last line, should just be a bunch of --- lines
    hline = hline.strip()
    if hline != '-'*_line_length(colwidth):
        raise RuntimeError, 'header format not what is expected'
    return colwidth, names, dtypes, units, ind    


def _read_file(f):
    """
    Loads the data from a table file
    """
    close_file = False
    if isinstance(f, str):
        close_file = True
        # Try loading the file
        f = open(f, 'r')
    # Read the data
    try:
        colwidth, names, dtypes, units, ind = _read_header(f)
        data_str = f.readlines()
        nrow = len(data_str)
        n_input_col = len(colwidth)
        data_array = np.zeros([nrow, n_input_col], dtype='|S{0}'.format(max(colwidth)))    
        slice_line = _line_slicer(colwidth)
        for i, line in enumerate(data_str):
            data_array[i] = slice_line(line, strip=True)
    finally:
        if close_file:
            f.close()
            
    return colwidth, names, dtypes, units, ind, data_array
    
def read_table(f):
    """
    Read a SimArray table stored by write_table()
    
    f should be a filename or a buffer
    
    returns a dict
    """
    colwidth, names, dtypes, units, ind, data_array = _read_file(f)
    nrow, n_input_col = data_array.shape
    
    # Infer array shapes
    num_cols = np.bincount(ind[:,0])
    num_arrays = len(num_cols)
    # Initialize arrays
    array_dtype = num_arrays * [None]
    array_units = num_arrays * [None]
    array_names = num_arrays * [None]
    for i, dtype, unit, name in zip(ind[:,0], dtypes, units, names):
        array_dtype[i] = dtype
        array_units[i] = unit
        array_names[i] = name
    arrays = []
    for ncol, dtype, unit in zip(num_cols, array_dtype, array_units):
        array = SimArray(np.zeros([nrow, ncol], dtype=dtype), unit)
        arrays.append(array)
    # Map data to arrays
    for input_col in range(n_input_col):
        dtype = dtypes[input_col]
        iArray, iCol = ind[input_col]
        column = data_array[:, input_col]
        if dtype.startswith('|S'):
            # handle strings
            column = [val.strip()[1:-1] for val in column]
        arrays[iArray][:, iCol] = column
        
    return dict(zip(array_names, arrays))

def write_table(f, x, *args, **kwargs):
    """
    write_table(f, x, (x1, x2, ...), **kwargs) saves a human-readable ASCII table
    from arrays/SimArrays.
    
    Parameters
    ----------
    f : str or buffer
        Filename or open buffer to save to
    x, (x1, x2, ...) : array-like or SimArray or dict
        Data to save.  Must be 1-D or 2-D.  The length of all arrays must
        be the same.  Must be shape (nrow,) or (nrow, ncol).
        x can optionally be a dict of arrays.  The optional args will be
        ignored.
        
    kwargs 
    ------
    f : str or filepointer
        File or filename to save the table to.
    col_labels : list-like, optional
        Names of the columns.  If x is a dict, then x.keys() is used.  
        Defaults to 0, 1, 2...
    
    Notes
    -----
    This is NOT optimized for large arrays
    """
    import json
    
    # Parse kwargs
    col_labels = kwargs.get('col_labels', None)
    
    # Initial data formatting
    if isinstance(x, dict):
        col_labels = x.keys()
        arrays = x.values()
    else:
        arrays = [x] + list(args)
        if col_labels is None:
            col_labels = np.arange(len(arrays))
    
    # Deal with multiple-column data
    for i, array in enumerate(arrays):
        array = np.asanyarray(array)
        if np.ndim(array) == 1:
            arrays[i] = array[:,None]
        elif np.ndim(array) != 2:
            raise ValueError, 'Can only write tables for 1-D or 2-D arrays'
    save_arrays = []
    save_labels = []
    array_ind = []
    for i, array in enumerate(arrays):
        # Concatenate along columns
        col_label = col_labels[i]
        for j in range(array.shape[1]):
            save_labels.append('"{0}"'.format(col_label))
            save_arrays.append(array[:, j])
            array_ind.append(json.dumps([i, j]))
            
    dtypes = [array.dtype for array in save_arrays]
    dtypes = [dtype.str for dtype in dtypes]
    # Get units
    units = [get_units(array) for array in save_arrays]
    units = [str(unit) for unit in units]
    # Convert arrays to strings
    header = zip(save_labels, dtypes, units, array_ind)
    header = [[str(a) for a in header_i] for header_i in header]
    array_str = []
    for array, dtype in zip(save_arrays, dtypes):
        if '|S' in dtype:
            # we have a string
            string = ['"{0}"'.format(a) for a in array]
        else:
            string = [str(a) for a in array]
        array_str.append(string)
    
    colwidth = [max([len(a) for a in array_str_i]) for array_str_i in array_str]
    headcolwidth = [max([len(a) for a in array_str_i]) for array_str_i in header]
    colwidth = [max(c, h) for c, h in zip(colwidth, headcolwidth)]
    colwidth = [width + _COLPAD for width in colwidth]
    colwidth[0] += 1
    
    print_str = ''
    ncols = len(save_arrays)
    # Format header
    for i in xrange(len(header[0])):
        line = ''
        for j in xrange(ncols):
            width = colwidth[j]
            line += header[j][i].center(width)
        line = '#' + line[1:]
        print_str += line + '\n'
    hline = '#' + '-' * np.sum(colwidth) + '\n'
    print_str += hline
    print_str = '# ' + json.dumps({'colwidth': colwidth}) + '\n' + print_str
    # Format arrays
    for i in xrange(len(array_str[0])):
        for j in xrange(ncols):
            width = colwidth[j]
            print_str += array_str[j][i].center(width)
        print_str += '\n'
    
    if isinstance(f, str):
        with open(f, 'w') as fp:
            fp.write(print_str)
    else:
        fp.write(print_str)