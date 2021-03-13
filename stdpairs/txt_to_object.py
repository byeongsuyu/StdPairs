r"""
txt_to_object

A module saving work variables in ``stdpairs`` as a text to archive it. 
``txt_to_affinemonoid()`` and ``txt_to_monomialideal()`` methods appear for the end users.

AUTHORS:

- Byeongsu Yu (2021-02-25): initial version.

"""
#*****************************************************************************
#       Copyright (C) 2021 Byeongsu Yu <byeongsu.yu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; version 3 of
#  the License. 
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import numpy as np
import json
import ast
from sage.all import matrix
from sage.all import ZZ

from . import affinemonoid
from . import monomialideal
from . import properpair
from . import _string_interface

def txt_to_affinemonoid(text_info):
    r"""
    Given a string ``text_info`` generated by a method  ``save_txt()`` of an ``AffineMonoid`` object,
    this method returns an ``AffineMonoid`` object which is the same as the original ``AffineMonoid`` object.

    INPUT:

    - ``text_info`` -- A ``string`` objct generated by the method ``AffineMonoid.save_txt()``.

    OUTPUT:

    - An ``AffineMonoid`` object.

    EXAMPLE::

        sage: from stdpairs import AffineMonoid, txt_to_affinemonoid
        sage: R = AffineMonoid(matrix(ZZ, [[2,0,1],[0,1,1]])) 
        sage: R.save_txt()                                                                                                                              
        'Q\n2,0,1|0,1,1\n'
        sage: L = txt_to_affinemonoid('Q\n2,0,1|0,1,1\n')
        sage: R == L                                                                                                                                    
        True

    """
    stripped_info = text_info.split("\n")
    if stripped_info[0] != "Q":
        raise ValueError("[Error]: Given file is not a valid save file")
    if stripped_info[1] == '':
        return affinemonoid.AffineMonoid(matrix(ZZ,0))
    else:
        return affinemonoid.AffineMonoid(_string_interface._string_to_np2d(stripped_info[1]))

def txt_to_monomialideal(text_info):
    r"""
    Given a string ``text_info`` generated by ``MonomialIdeal.save_txt()`` method,
    returns an ``MonomialIdeal`` object which is the same as the original ``MonomialIdeal`` object.

    INPUT:

    - ``text_info`` -- A ``string`` objct generated by the method ``MonomialIdeal.save_txt()``.

    OUTPUT:

    - An ``MonomialIdeal`` object.

    EXAMPLE::

        sage: from stdpairs import AffineMonoid, MonomialIdeal, txt_to_monomialideal
        sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
        sage: Q = AffineMonoid(A)
        sage: I = MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q)
        sage: I.standard_cover()
        Calculate the standard cover of an ideal
        It takes a few minutes, depending on the system.
        Cover for 1  generator was calculated.  2  generators remaining. 
        Cover for 2  generators was calculated.  1  generators remaining. 
        Cover for 3  generators was calculated.  0  generators remaining. 
        {(0, 3): [[([[0], [0], [0]]^T,[[0, 0], [0, 1], [1, 1]])],
        [([[1], [1], [1]]^T,[[0, 0], [0, 1], [1, 1]]),
         ([[1], [0], [1]]^T,[[0, 0], [0, 1], [1, 1]])]]}
        sage: str_I = I.save_txt()                                                                                                                              
        sage: J = txt_to_monomialideal(str_I)
        sage: J == I                                                                                                                                    
        True
        sage: J.standard_cover()
        # This does not require re-calculation.                                                                                                                        
        {(0, 3): [([[1], [1], [1]]^T,[[0, 0], [0, 1], [1, 1]]),
          ([[0], [0], [0]]^T,[[0, 0], [0, 1], [1, 1]]),
          ([[1], [0], [1]]^T,[[0, 0], [0, 1], [1, 1]])]}
    """
    stripped_info = text_info.split("\n")
    if stripped_info[0] != "I":
        raise ValueError("[Error]: Given string is not a valid text storing MonomialIdeal object")
    if stripped_info[1] == "Q":
        ambient_monoid = affinemonoid.AffineMonoid(_string_interface._string_to_np2d(stripped_info[2]))
    else:
        raise ValueError("[Error]: Given file is not a valid save file")
    if stripped_info[4] == "1":
        return monomialideal.MonomialIdeal(np.array([]).astype('int64'), ambient_monoid)
    an_ideal = monomialideal.MonomialIdeal(_string_interface._string_to_np2d(stripped_info[6]), ambient_monoid) 
    remaining = stripped_info[7:]
    # Now deal with information part.
    while (len(remaining) >0):
        info_text = remaining[0]
        del remaining[0]
        if info_text == "__is_std_cover_calculated":
            if remaining[0] == "1":
                del remaining[0]
                str_cover = json.loads(remaining[0])
                del remaining[0]
                cover = {}
                for key, value in str_cover.items():
                    cover[ast.literal_eval(key)] = [properpair.ProperPair(_string_interface._string_to_np2d(pair.split("&")[0]),ast.literal_eval(pair.split("&")[1]),an_ideal) for pair in value]
                an_ideal._manually_set_attributes('__dict_standard_pairs', cover)
            else:
                del remaining[0]
        if (info_text == "__is_overlap_calculated") or (info_text == "__is_max_overlap_calculated"):
            if remaining[0] == "1":
                del remaining[0]
                str_cover = json.loads(remaining[0])
                del remaining[0]
                ovClasses = {}
                for key, value in str_cover.items():
                    ovClasses[ast.literal_eval(key)] =[]
                    for ov_class in value:
                        ovClasses[ast.literal_eval(key)].append([properpair.ProperPair(_string_interface._string_to_np2d(pair.split("&")[0]),ast.literal_eval(pair.split("&")[1]),an_ideal) for pair in ov_class])
                if (info_text == "__is_overlap_calculated"):
                    an_ideal._manually_set_attributes('__dict_overlap_classes', ovClasses)
                else:
                    an_ideal._manually_set_attributes('__dict_max_overlap_classes',ovClasses)
        if (info_text == "__is_ass_prime_calculated"):
            if remaining[0] == "1":
                del remaining[0]
                str_ass_primes =json.loads(remaining[0])
                del remaining[0]
                ass_primes={}
                for key, value in str_ass_primes.items():
                    ass_primes[ast.literal_eval(key)] = txt_to_monomialideal(value)
                an_ideal._manually_set_attributes('__is_ass_prime_calculated', ass_primes)
                
        if info_text == "__is_irr_decom_prime_calculated":
            if remaining[0] == "1":
                del remaining[0]
                str_irrDecoms = json.loads(remaining[0])
                del remaining[0]
                irrDecoms = [txt_to_monomialideal(txt_ideal) for txt_ideal in str_irrDecoms]
                an_ideal._manually_set_attributes('__list_irreducible_decom_ideals', irrDecoms)
                
    return an_ideal

def _cover_of_stdpairs_to_txt(cover, ambient_monoid):
    r"""
    Given a ``dictionary`` object ``cover``, whose keys are faces of ``ambient_monoid`` and whose values are ``ProperPair`` objects,
    return a string which stores the information of ``cover``. This can be recovered by ``txt_to_cover_of_stdpairs()`` method.
    Notes that an ambient ideal of the pair in ``cover`` was lost. In that case, use ``save_txt()`` in ``MonomialIdeal`` object instead.

    INPUT:
    
    - ``cover`` -- A ``dictionary`` object whose keys are tuples representing face of ``ambient_monoid`` and whose values are ``list`` of ``ProperPair`` objects.
    - ``ambient_monoid`` -- An ``AffineMonoid`` object.

    OUTPUT:

    - A ``string`` object.

    EXAMPLE::

        sage: from stdpairs import AffineMonoid, MonomialIdeal, cover_of_stdpairs_to_txt, 
        ....: txt_to_cover_of_stdpairs
        sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
        sage: Q = AffineMonoid(A)
        sage: I = MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q)
        sage: C=I.standard_cover()
        Calculate the standard cover of an ideal
        It takes a few minutes, depending on the system.
        Cover for 1  generator was calculated.  2  generators remaining. 
        Cover for 2  generators was calculated.  1  generators remaining. 
        Cover for 3  generators was calculated.  0  generators remaining.
        sage: C 
        {(0, 3): [[([[0], [0], [0]]^T,[[0, 0], [0, 1], [1, 1]])],
        [([[1], [1], [1]]^T,[[0, 0], [0, 1], [1, 1]]),
         ([[1], [0], [1]]^T,[[0, 0], [0, 1], [1, 1]])]]}
        sage: C[(0,3)][0].ambient_ideal()                                                 
        An ideal whose generating set is 
        [[2 2 2]
         [0 1 2]
         [2 2 2]]
        sage: str_C = cover_of_stdpairs_to_txt(C, I.ambient_monoid())   
        sage: str_C                                                                       
        'C\nQ\n0,1,1,0|0,0,1,1|1,1,1,1\n{"(0, 3)": ["1|1|1&(0, 3)", "0|0|0&(0, 3)", "1|0|1&(0, 3)"]}'
        sage: D = txt_to_cover_of_stdpairs(str_C)
        sage: D
        [An affine semigroup whose generating set is 
         [[0 1 1 0]
          [0 0 1 1]
          [1 1 1 1]],
         {(0, 3): [([[1], [1], [1]]^T,[[0, 0], [0, 1], [1, 1]]),
           ([[0], [0], [0]]^T,[[0, 0], [0, 1], [1, 1]]),
           ([[1], [0], [1]]^T,[[0, 0], [0, 1], [1, 1]])]}]
        sage: D[1][(0,3)][0].ambient_ideal() 
        An ideal whose generating set is 
        []
        # An information about ambient ideal was deleted;
        # 

    """
    """Save the cover; all pairs will be saved with proper pair of zero ideal.
    If you know that covers are standard, then use the save ideal instead of it."""
    zero_ideal = monomialideal.MonomialIdeal(matrix(ZZ,0),ambient_monoid)
    if not isinstance(cover, dict):
        raise ValueError("[Error]: Given 1st instance is not a dictionary")
    if not isinstance(ambient_monoid, affinemonoid.AffineMonoid):
        raise ValueError("[Error]: Given 2nd instance is not an AffineMonoid object")
    return "C\nQ\n"+_string_interface._np2d_to_string( ambient_monoid.gens())+"\n"+_string_interface._json_dump_cover(cover)

def _txt_to_cover_of_stdpairs(text_info):
    r"""
    Given a string generated by ``_cover_of_stdpairs_to_txt()`` method,
    return a ``dictionary`` object ``cover``, whose keys are faces of an affine monoid and whose values are ``ProperPair`` objects.
    
    INPUT:

    - A ``string`` object generated by ``_cover_of_stdpairs_to_txt()`` method.

    OUTPUT:
    
    - ``cover`` -- A ``dictionary`` object whose keys are tuples representing face of ``ambient_monoid`` and whose values are ``list`` of ``ProperPair`` objects.


    EXAMPLE::

        sage: from stdpairs import AffineMonoid, MonomialIdeal, cover_of_stdpairs_to_txt, 
        ....: txt_to_cover_of_stdpairs
        sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
        sage: Q = AffineMonoid(A)
        sage: I = MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q)
        sage: C=I.standard_cover()
        Calculate the standard cover of an ideal
        It takes a few minutes, depending on the system.
        Cover for 1  generator was calculated.  2  generators remaining. 
        Cover for 2  generators was calculated.  1  generators remaining. 
        Cover for 3  generators was calculated.  0  generators remaining.
        sage: C 
        {(0, 3): [[([[0], [0], [0]]^T,[[0, 0], [0, 1], [1, 1]])],
        [([[1], [1], [1]]^T,[[0, 0], [0, 1], [1, 1]]),
         ([[1], [0], [1]]^T,[[0, 0], [0, 1], [1, 1]])]]}
        sage: C[(0,3)][0].ambient_ideal()                                                 
        An ideal whose generating set is 
        [[2 2 2]
         [0 1 2]
         [2 2 2]]
        sage: str_C = _cover_of_stdpairs_to_txt(C, I.ambient_monoid())   
        sage: str_C                                                                       
        'C\nQ\n0,1,1,0|0,0,1,1|1,1,1,1\n{"(0, 3)": ["1|1|1&(0, 3)", "0|0|0&(0, 3)", "1|0|1&(0, 3)"]}'
        sage: D = _txt_to_cover_of_stdpairs(str_C)
        sage: D
        [An affine semigroup whose generating set is 
         [[0 1 1 0]
          [0 0 1 1]
          [1 1 1 1]],
         {(0, 3): [([[1], [1], [1]]^T,[[0, 0], [0, 1], [1, 1]]),
           ([[0], [0], [0]]^T,[[0, 0], [0, 1], [1, 1]]),
           ([[1], [0], [1]]^T,[[0, 0], [0, 1], [1, 1]])]}]
        sage: D[1][(0,3)][0].ambient_ideal() 
        An ideal whose generating set is 
        []
        # An information about ambient ideal was deleted;
        # 

    """
    stripped_info = text_info.split("\n")
    ambient_monoid = affinemonoid.AffineMonoid(_string_interface._string_to_np2d(stripped_info[2]))
    str_cover = json.loads(stripped_info[3])
    zero_ideal = monomialideal.MonomialIdeal(np.array([]).astype('int64'), ambient_monoid)
    cover = {}
    for key, value in str_cover.items():
        cover[ast.literal_eval(key)] = [properpair.ProperPair(_string_interface._string_to_np2d(pair.split("&")[0]),ast.literal_eval(pair.split("&")[1]),zero_ideal) for pair in value]
    return [ambient_monoid, cover]
