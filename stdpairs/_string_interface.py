
#### Variable Transformations
"""
    string
        _np2d_to_string(numpy 2D array):
            return hashstring of the numpy 2D array.

    numpy 2D vector
        _string_to_np2d(string):
            Given string formed by _np2d_to_string(array) method,
            return corresponding numpy 2D array.
        _json_dump_overlap_classes(cover Ov)
            Change a cover of overlapped classes Ov as its string form, and
            return the json file containing a dictionary of strings.
        _json_dump_cover(cover C):
            Change a cover of standard pairs C as its string form, and
            return the json file containing a dictionary of strings.
        _unique_np_arrays
            Given a list of numpy 2d array with the same shape, return the unique ones

            
"""

import numpy as np
import json

from . import properpair

def _json_dump_overlap_classes(ov_cover):
    """Return string converted cover."""
    if not isinstance(ov_cover, dict):
        raise ValueError("[Error]: Given instance is not a dictionary")
    str_cover = {}
    for key,value in ov_cover.items():
        str_cover[str(key)] = []
        for ov_class in value:
            str_cover[str(key)].append([pair._hashstr() for pair in ov_class])
    return json.dumps(str_cover)

def _json_dump_cover(cover):
    """Return string converted cover."""
    if not isinstance(cover, dict):
        raise ValueError("[Error]: Given instance is not a dictionary")
    str_cover = {}
    for key,value in cover.items():
        str_cover[str(key)]=[pair._hashstr() for pair in value]
    return json.dumps(str_cover)

def _np2d_to_string(twod_array):
    if not isinstance(twod_array,np.ndarray):
        raise ValueError("[Error]: This instance is not numpy.ndarray type")
    if twod_array.shape == (0,):
        return ""
    if len(twod_array.shape) !=2:
        raise ValueError("[Error]: This array is not 2-Dimensoinal")
    return_str = ""
    for row in twod_array:
        for val in row:
            return_str = return_str + str(val)+","
        return_str = return_str[:-1] + "|"
    return return_str[:-1]

def _string_to_np2d(twod_string):
    if not isinstance(twod_string,str):
        raise ValueError("[Error]: This does not give string object")
    result= twod_string.split('|')
    new_result = [ list(map(int, item.split(','))) for item in result]
    return np.array(new_result).astype('int64')

def _unique_np_arrays(list_of_np_array):
    """Given a list of numpy 2d array with the same shape, return the unique ones"""
    if not isinstance(list_of_np_array, list):
        raise ValueError("[Error]: Given argument is not a list")
    for ind in list_of_np_array:
        if not isinstance(ind,np.ndarray):
            raise ValueError("[Error]: List contains non-np.array")
    if list_of_np_array ==[]:
        return []
    sh = list_of_np_array[0].shape
    for ind in list_of_np_array:
        if ind.shape != sh:
            raise ValueError("[Error]: Given numpy array has distinct shapes")
    list_of_strings = list(set([_np2d_to_string(item) for item in list_of_np_array]))
    return [_string_to_np2d(item) for item in list_of_strings]
