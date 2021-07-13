r"""
txt_to_object

A module provides methods to communicate with Macaulay2. 

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
from sage.all import macaulay2
from sage.all import matrix
from sage.all import ZZ

from . import affinemonoid
from . import monomialideal

def from_macaulay2(var_name):
    r"""
    Given a ``macaulay2`` type object ``mac2_mon_subalg`` which saves the `MonomialSubalgebra <http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.12/share/doc/Macaulay2/Normaliz/html/___Monomial__Subalgebra.html>`__ object,
    this method returns ``AffineMonoid`` object corresponding to the `MonomialSubalgebra <http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.12/share/doc/Macaulay2/Normaliz/html/___Monomial__Subalgebra.html>`__. 

    INPUT:

    - ``var_name`` -- A ``string`` object which equal to the variable name of a `MonomialSubalgebra` class variable created by ``Normaliz`` package in `Macaulay2 <http://www.math.uiuc.edu/Macaulay2/>`__.

    OUTPUT:

    - An ``AffineMonoid`` object.

    EXAMPLE::

        sage: from stdpairs import from_macaulay2
        sage: R = macaulay2.eval('ZZ[x,y,z]')
        sage: temp=macaulay2.needsPackage('"Normaliz"')
        sage: temp=macaulay2.eval('S=createMonomialSubalgebra {x^5*y, y*z^2, z^3}')
        sage: Q = from_macaulay2('S')
        sage: Q
        An affine semigroup whose generating set is
        [[5 0 0]
         [1 1 0]
         [0 2 3]]

    """
    test = macaulay2.eval("L = gens "+var_name)
    numgens = eval(macaulay2.eval("#L"))
    gens=[]
    for idx in range(numgens):
        gens.append(macaulay2.eval("exponents L#"+str(idx)).split('\n')[0].replace("{","[").replace("}","]"))
    np_gens = np.concatenate([np.array(ast.literal_eval(item), dtype='int64').T for item in gens ],axis=1)
    return affinemonoid.AffineMonoid(np_gens)

def to_macaulay2(monomial_ideal, ring_name="R", ideal_name="I", std_cover_name = "SC"):
    r"""
    Given a ``MonomialIdeal`` object ``monomial_ideal``, 
    this method returns a ``dictionary`` object containing the ambient monoid, the given monomial ideal, and its standard cover as ``string`` objects.
    Users can access those ``Macaulay2`` variables using ``ring_name``, ``ideal_name``, and ``std_cover_name``.
    Notes that the ambient monoid and the monomial ideal will be delivered as
    ``macaulay2`` type object, especially with `MonomialSubalgebra <http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.12/share/doc/Macaulay2/Normaliz/html/___Monomial__Subalgebra.html>`__  type.
    See `MonomialSubalgebra <http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.12/share/doc/Macaulay2/Normaliz/html/___Monomial__Subalgebra.html>`__ for detail.
    The standard cover will be delivered as a ``list`` object in Macaulay 2. Also, this method will executes ``monomial_ideal.standard_cover()`` if it was not executed previously.

    INPUT:

    - ``monomial_ideal`` -- A ``MonomialIdeal`` type object.

    OUTPUT:

    - A ``dictionary`` object as below.

    ::

        { 'AffineSemigroupRing':mac2_mon_subalg, 'MonomialIdeal':mac2_mon_sub_alg_ideal, 'StandardCover':mac2_std_cover}

    where ``mac2_mon_subalg``, ``mac2_mon_sub_alg_ideal``, and ``mac2_std_cover`` are ``string`` object.

    EXAMPLE::

        sage: from stdpairs import AffineMonoid,MonomialIdeal, to_macaulay2
        sage: Q=AffineMonoid(matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])) 
        sage: I=MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q) 
        sage: S=to_macaulay2(I)
        sage: S                                                                           
        {'AffineSemigroupRing': ZZ[c, a*c, a*b*c, b*c]
        <BLANKLINE>
        monomial subalgebra of PolyRing,
        'MonomialIdeal':   2 2   2   2   2 2 2
        {a c , a b*c , a b c }
        <BLANKLINE>
        List,
        'StandardCover': {{1, {c, b*c}}, {a*c, {c, b*c}}, {a*b*c, {c, b*c}}}
        <BLANKLINE>
        List}
        sage: # To access values in dictionary via Macaulay2,
        sage: mac_ring=macaulay2.eval("R")                                                       
        sage: mac_ideal=macaulay2.eval("I")                                                       
        sage: mac_sc=macaulay2.eval("SC")                                                      

    """
    
    cover = monomial_ideal.standard_cover()
    if not isinstance(monomial_ideal, monomialideal.MonomialIdeal):
        raise ValueError("[Error]: 1st argument should be type of MonomialIdeal.")
    affine_monoid = monomial_ideal.ambient_monoid()
    # Number of attributes we need.
    numvar = (affine_monoid.gens()).shape[0]
    #Generate polynomial ring
    gens_poly_ring='PolyRing=ZZ[vars (0..'+str(numvar-1)+')]'
    temp = macaulay2.eval(gens_poly_ring)
    #Load normaliz in Macaulay2
    temp = macaulay2.needsPackage('"Normaliz"')
    #Generate monomial subalgebra in Macaulay2.
    mon_sub_alg = ring_name+'=createMonomialSubalgebra {'
    list_of_gens =[]
    for column in (affine_monoid.gens()).T:
        eq =''
        for jdx in range(numvar):
            eq = eq + 'PolyRing_'+str(jdx)+'^'+str(column[jdx])+'*'
        eq = eq[:-1]+','
        mon_sub_alg = mon_sub_alg + eq
        list_of_gens.append(eq[:-1])
    if mon_sub_alg != '{':
        mon_sub_alg = mon_sub_alg[:-1]+'}'
    else:
        mon_sub_alg = mon_sub_alg+'}'
    mac2_mon_subalg=macaulay2.eval(mon_sub_alg)

    # Generate monomial ideal of monomial subalgebra in Macaulay2 as a list
    mon_sub_alg_ideal =ideal_name+'={'
    for column in (monomial_ideal.gens()).T:
        eq =''
        for jdx in range(numvar):
            eq = eq + 'PolyRing_'+str(jdx)+'^'+str(column[jdx])+'*'
        eq = eq[:-1]+','
        mon_sub_alg_ideal = mon_sub_alg_ideal + eq
    if mon_sub_alg_ideal != '{':
        mon_sub_alg_ideal = mon_sub_alg_ideal[:-1]+'}'
    else:
        mon_sub_alg_ideal = mon_sub_alg_ideal+'}'
    mac2_mon_sub_alg_ideal = macaulay2.eval(mon_sub_alg_ideal)   
    
    
    total_list =std_cover_name+'={'
    for key, val in sorted(cover.items(), key=lambda e: e[0][0]):
        faces = [list_of_gens[idx] for idx in key]
        string_faces = '{'+ ','.join(faces) +'}'
        for pair in val:
            gens = ''
            for kdx in range(numvar):
                gens = gens+'PolyRing_'+str(kdx)+'^'+str(pair.monomial()[kdx,0])+'*'
            gens = gens[:-1]
            std_pair = '{'+gens+','+string_faces+'}'
            total_list = total_list + std_pair +',' 
    total_list = total_list[:-1]+'}'
    mac2_std_cover =macaulay2.eval(total_list)
    return { 'AffineSemigroupRing':mac2_mon_subalg, 'MonomialIdeal':mac2_mon_sub_alg_ideal, 'StandardCover':mac2_std_cover}

