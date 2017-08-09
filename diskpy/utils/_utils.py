# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:19:03 2015

@author: ibackus
"""
import warnings
import logging
import os
import glob
import fnmatch
import numpy as np
import pynbody as pb
SimArray = pb.array.SimArray
import logging

def snap_param(snapshot):
    """
    Converts the param dict in a pynbody SimSnap to a more useable format 
    (same as configparser)
    
    Parameters
    ----------
    snapshot : SimSnap
        SimSnap, should be loaded from disk, tipsy format assumed
    
    Returns
    -------
    param_dict : dict
        Same format as output of configparser
    """
    if not hasattr(snapshot, '_paramfile'):
        raise RuntimeError, "Snapshot does not have attribute _paramfile"
    
    pnew = {}
    for k, v in snapshot._paramfile.iteritems():
        pnew[k] = str2num(v)
    
    if 'filename' in pnew:
        pnew.pop('filename', None)
    
    return pnew

def as_simarray(x, assumed_units='1'):
    """
    Return x as a SimArray.  x can be a number, a numpy array, a SimArray, 
    an interpretable string (e.g. '1 mm' or '2e3 au').
    
    Parameters
    ----------
    x : object
        object to make a simarray from.  If already a simarray, x is returned
    assumed_units : unit-like
        The units used if x has no units.  Can be a string that pynbody units 
        takes (e.g. 'au' or 'Msol yr**-1') or a pynbody unit.  If None, no 
        units will be assumed.
        
    Returns
    -------
    xsimarray : SimArray
        x as a SimArray
    """
    
    if isinstance(x, (str, pb.units.UnitBase)):
        
        x1 = SimArray(1.0, x)
        
    elif isinstance(x, SimArray):
        
        x1 = x
        
    else:
        
        x1 = SimArray(x)
        
    if not pb.units.has_units(x) and assumed_units is not None:
        
        x1.units = pb.units.Unit(assumed_units)
        
    return x1

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
        farray = np.genfromtxt(fname,delimiter='=',dtype='|S256')
        if farray.size==2:
            # Handle the case that there is only 1 line the .param file
            farray = farray.reshape([1, 2])
        for n in range(len(farray)):
            param[farray[n,0].strip()] = str2num(farray[n,1].strip())
    # --------------------------------------------------
    # Parse director file
    # --------------------------------------------------
    elif ftype == 'director':
        f = open(fname,'r')
        f.seek(0)
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
    
def logparser(fname, verbose=False):
    """
    Parses a ChaNGa log file to find run-time parameters.  Also returns the 
    header under the key 'header'
    
    Parameters
    ----------
    
    fname : str
        Filename of the log file to open
    verbose : bool
        (optional) If True, prints all the parameters
    
    Returns
    -------
    
    param : dict
        Dictionary of the parameters
        
    See Also
    --------
    
    configparser
    """
    header = []
    with open(fname,'r') as f:
        
        # Parse until finding parameters
        found = False
        while not found:
            
            line = f.readline().strip()
            if line[0] != '#':
                
                raise RuntimeError, 'Could not find parameters'
                
            line = line.strip('#').strip()
            
            if line == 'Parameters:':
                
                found = True
                
            else:
                
                header.append(line)
        
        # Now read in parameters
        done = False
        param = {'header': header}
        while not done:
            
            line = f.readline().strip()
            if line[0] != '#':
                
                raise RuntimeError, 'Expected # at beginning of line: ' + line
                
            if ':' not in line:
                
                done = True
                
            else:
                
                line = line.strip('#').strip()
                k, v = line.split(': ')
                v = str2num(v)
                param[k] = v
                
                if verbose:
                    
                    print k, v
                    
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
    
def units_from_param(param):
    raise RuntimeError, "diskpy.utils.units_from_param is depcreated.  "\
    "use diskpy.pychanga.units_from_param"
        
def get_units(x):
    """
    Retrieves the units of x if:
        * x has a unit
        * x is a unit
        * x is a unit string (see pynbody)
    
    Else, returns None
    """
    
    if pb.units.is_unit(x):
        
        unit = x
        
    elif pb.units.has_unit(x):
        
        unit = x.units
        
    elif isinstance(x, str):
        
        unit = pb.units.Unit(x)
        
    else:
        
        unit = None
        
    return unit
    
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
            
            x_i = np.asarray(x_i)

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
        
        x = np.asarray(x)
        
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

            if pb.units.has_units(x_i):

                x_i_array = strip_units(x_i)
                x_out.append(SimArray(x_i_array, units[i]))

            else:

                x_out.append(SimArray(x_i, units[i]))

    else:

        if pb.units.has_units(x):

            x_array = strip_units(x)
            x_out = SimArray(x_array, units)

        else:

            x_out = SimArray(x, units)

    return x_out
    
def match_units(x, y):
    """
    Matches the units of x to y and returns x and y in the same units.

    IF x and y don't have units, they are unchanged

    IF one of x or y has units, the unitless quantity is returned as a
    SimArray (see pb.array.SimArray) with the units of the other quantity.

    IF both have units, then an attempt is made to convert x into the units of
    y.  If this is not possible, an error is raised, for example if x is in
    units of 'au' and y is in units of 'Msol'

    x, y can be: scalar, array, SimArray, pynbody unit (eg pb.units.G),
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
    # If one is a named unit (eg pb.units.G), convert to SimArray
    if isinstance(x, pb.units.UnitBase):

        x = SimArray(1.0, x)

    if isinstance(y, pb.units.UnitBase):

        y = SimArray(1.0, y)

    # ----------------------------------------------
    # Check the units
    # ----------------------------------------------
    # If both have units, try to convert x to the units of y
    if (pb.units.has_units(x)) & (pb.units.has_units(y)):

        x_out = (x.in_units(y.units))
        y_out = y

    # If only x has units, make y a SimArray with the units of x
    elif (pb.units.has_units(x)):

        y_out = SimArray(y, x.units)
        x_out = x

    # If only y has units, make x a SimArray with the units of y
    elif (pb.units.has_units(y)):

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
        useLogger = True
        logger = logging.getLogger('pynbody')
        current_verbosity = logger.getEffectiveLevel()

    else:

        # For pynbody version < 0.3, verbosity is handled in the config
        useLogger = False
        current_verbosity = pb.config['verbose']

    # -----------------------------
    # Change verbosity
    # -----------------------------
    if cmd is None:
        # Don't change verbosity.  just return the current verbosity
        pass

    elif cmd == 'off':
        # Toggle verbosity off
        if useLogger:

            logger.setLevel(logging.ERROR)

        else:

            pb.config['verbose'] = False

    elif cmd == 'on':
        # Toggle verbosity on
        if useLogger:

            logger.setLevel(logging.DEBUG)

        else:

            pb.config['verbose'] = True

    else:
        # Set verbosity to the verbosity level specified by cmd
        if useLogger:

            logger.setLevel(cmd)

        else:

            pb.config['verbose'] = cmd

    # Return the verbosity level before any changes were made
    return current_verbosity
    
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
    
def get_module_names(fname):
    """
    An import utility that returns the module names in the directory of file.  
    Ignores filenames beginning with an underscore.
    
    Parameters
    ----------
    
    fname : str
        Filename
    
    Returns
    -------
    
    modulenames : list
        A list of the modules
    """
    directory = os.path.dirname(os.path.realpath(fname))
    searchstr = os.path.join(directory, '*.py')
    fullpaths = glob.glob(searchstr)
    
    fnames = []
    
    for fullpath in fullpaths:
        
        f = os.path.split(fullpath)[-1]
        
        if f[0] is not '_':
            
            fnames.append(fullpath)
    
    modulenames = []
    
    for fname in fnames:
        
        modulename = os.path.split(os.path.splitext(fname)[0])[-1]
        modulenames.append(modulename)
        
    return modulenames
    
def which(cmd):
    """
    Pythonic equivalent of UNIX which
    
    Parameters
    ----------
    
    cmd : str
        Command to find
    
    Returns
    -------
    
    path : str -or- None
        Full path to command which would be executed in current shell by cmd
        If the command cannot be found, None is returned
    """
    
    cmd_loc = None
    for path in os.getenv('PATH').split(os.path.pathsep):
        
        fullpath = os.path.join(path, cmd)
        if os.path.exists(fullpath):
            
            cmd_loc = fullpath
            break
        
    return cmd_loc

def deepreload(module):
    """
    Convenience package for reloading diskpy modules when working.
    
    Parameters
    ----------
    
    module : str
        Module to reload
    
    Returns
    -------
    
    command : str
        Command to exec which will reload the module and its parents
        
    Examples
    --------
    
    >>> command = deepreload('diskpy.disk._spirals')
    >>> exec(command)
    
    Or, equivalently:
    
    >>> exec(deepreload('diskpy.disk._spirals'))
    
    This is equivalent to:
    
    >>> reload(diskpy.disk._spirals)
    >>> reload(diskpy.disk)
    >>> reload(diskpy)
    """
    
    if not isinstance(module, str):
        
        raise ValueError('module must be a str')
        
    modules = module.split('.')
    nLevels = len(modules)
    
    command = ''
    for i in reversed(range(nLevels)):
        
        x = '.'.join(modules[0:i+1])
        command += 'reload(' + x + '); '
    
    return command
    
class logPrinter():
    """
    A simple class to print both to screen and to a log file
    
    logPrinter(verbose=True, filename=None, overwrite=True)
    
    Parameters
    ----------
    verbose : bool
        print to screen
    filename : str
        File to save to (if None, don't save to file)
    overwrite : bool
        Overwrite file if it exists
        
    printing can be used via the printer method or by calling the logprinter
    """    
    def __init__(self, verbose=True, filename=None, overwrite=True):
        
        if filename is not None:
            
            self.savetofile = True
            
            if overwrite:
                self.logfile = open(filename, 'w')
            else:
                self.logfile = open(filename, 'a')
        else:
            
            self.savetofile = False
                
        if verbose:
            if self.savetofile:
                self.printer = self._print
            else:
                self.printer = self._print_screen
        elif self.savetofile:
            self.printer = self._print_log
        else:
            self.printer = self._pass
                

    def __call__(self, string):
        
        return self.printer(string)
        
    def _pass(self, string):
        
        pass
    
    def _print_log(self, string):
        
        if not string.endswith('\n'):
            string += '\n'
        self.logfile.write(string)
        
    def _print_screen(self, string):
        
        print string
        
    def _print(self, string):
        
        self._print_screen(string)
        self._print_log(string)
                
    def close(self):
        """
        close the log file
        """
        if self.savetofile:
            self.logfile.close()
        
    def __del__(self):
        
        self.close()
    