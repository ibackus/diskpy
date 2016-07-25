# -*- coding: utf-8 -*-
"""
Defines the CustomList class, version of the python list which only accesses
data through __getitem__ and __setitem__, thus making it easy to subclass 
and control data types.

Created on Tue Jul 12 22:10:59 2016

@author: ibackus
"""
class CustomList():
    """
    customlist = CustomList(iterable)
    
    A modified implementation of a python list where all data access (i.e. 
    retreiving and setting items) takes place via _setOneItem and __getOneItem.
    This strategy makes this Class easy to subclass and control the behavior
    of setting/getting items.  To return a list of the same SubClass type,
    make sure to override the __getitem__ method (see __getitem__ in the code
    for CustomList as an example)
    
    All normal public list methods are implemented.  Currently the operands
    + and * are not implemented
    
    Printing (e.g. __repr__ and __str__) methods are implemented.  A custom
    class name can be printing by setting customlist._className = "class name"
    
    Examples
    --------
    Create a sub-class to apply some formatting before and after 
    setting/retrieving data:
    
    .. code-block:: python

        Class SubCustomList(CustomList):
            def __init__(self, iterable=[]):
                CustomList.__init__(self, iterable)
                self._className = "SubCustomList"
            def _setOneItem(ind, value):
                value = str(value).lower()
                CustomList._setOneItem(ind, value)
            def _getOneItem(ind):
                value = CustomList._getOneItem(ind)
                value = "{}: {}".format(ind, value)
                return value
    """
    def __init__(self, iterable=[]):
        
        # used for printing the CustomList
        self._className = "CustomList"
        # used for storing the data
        self._list = []
        self.extend(iterate(iterable))
            
    def __contains__(self, value):
        
        contains = False
        
        for item in self:
            
            if item == value:
                
                contains = True
                break
            
        return contains
        
    def __delitem__(self, ind):

        self._list.__delattr__(ind)
        
    def __getitem__(self, ind):
        
        if isinstance(ind, slice):
            
            indices = ind.indices(len(self))
            return CustomList([self._getOneItem(i) for i in range(*indices)])
            
        else:
            
            return self._getOneItem(ind)
        
    def __setitem__(self, ind, value):
        
        if isinstance(ind, slice):
            
            indices = ind.indices(len(self))
            
            for i, v in zip(range(*indices), value):
                
                self._setOneItem(i, v)
            
        else:
            
            self._setOneItem(ind, value)
        
    def __len__(self):
        
        return len(self._list)
        
    def __repr__(self):
        
        return reprIterable(self, self._className)
        
    def __str__(self):
        
        return strIterable(self, self._className)
        
    def _getOneItem(self, ind):
        """
        Handles all getting of data from the custom list.  Overwrite this to 
        change how a subclass handles getting data.
        """
        return self._list[ind]
        
    def _setOneItem(self, ind, value):
        """
        Handles all setting of data in the custom list.  Overwrite this to 
        change how a subclass handles setting data.
        """
        self._list[ind] = value
        
    def _addEmpty(self):
        
        self._list.append(None)
        
    def append(self, item):
        """
        Append item to end
        """
        self._addEmpty()
        self.__setitem__(-1, item)
        
    def count(self, value):
        """
        Return integer number of occurences of value
        """
        n = 0
        
        for item in self:
            
            n += (item == value)
            
        return n
        
    def extend(self, iterable):
        """
        Extend list by appending elements from the iterable
        """
        for item in iterable:
            
            self.append(item)
            
    def index(self, value):
        """
        Return first index of value.  Raises ValueError if value is not present
        """
        for i, item in enumerate(self):
            
            if value == item:
                
                return i
        
        raise ValueError, "Not in list"
        
    def insert(self, index, item):
        """
        Inserts item at index
        """
        if index < 0:
            
            index += len(self)
            
        if index < 0:
            
            index = 0
        
        elif index > len(self):
            
            index = len(self)
            
        self._addEmpty()
        self[index+1:] = self[index:-1]
        self[index] = item
            
    def pop(self, index):
        """
        Remove and return item located at index.  Raises IndexError if list is
        empty or index is out of range
        """
        value = self.__getitem__(index)
        self.__delitem__(index)
        return value
    
    def remove(value):
        """
        Remove first occurence of value.  Raises ValueError if the value is
        not present
        """
        index = self.index(value)
        dummy = self.pop(index)
        return
        
    def reverse(self):
        """
        Reverse items IN PLACE
        """
        self[:] = self[::-1]
        return
        
    def sort(self):
        """
        Sort items IN PLACE
        """
        sortedList = sorted(self)
        
        for index, value in enumerate(sortedList):
            
            self[index] = value

def reprIterable(iterable, classname="iterable"):
    """
    A function which returns a list-like string representation of an interable,
    with an class name
    """
    # Set up class name at beginning of string
    reprstr = ""
    if classname is not None:
        
        reprstr += "<" + classname + ": "
    
    # print items
    itemlist = list(iterable)
    reprstr += repr(itemlist) + ">"
    
    return reprstr
    
def strIterable(iterable, classname="iterable"):
    """
    A function which returns a list-like string representation of an interable,
    with an class name
    """
    # Set up class name at beginning of string
    string = ""
    if classname is not None:
        
        string += "<" + classname + ": "
    
    # print items
    itemlist = list(iterable)
    string += str(itemlist) + ">"
    
    return string
    
def iterate(x):
    """
    If x is not iterable, returns list(x).  Else, returns iter(x)
    """
    try:
            
        iterator = iter(x)
        
    except TypeError:
        
        iterator = [x]
        
    return iterator