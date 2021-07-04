r"""
stdpairs

Hidden module for main computation. Only ``pair_difference()`` method appears for the end users, which has its own documentation below.

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
import warnings
import itertools
import functools
from sage.all import ZZ
from sage.all import matrix
from sage.all import macaulay2
from sage.all import DiGraph
from sage.combinat.posets.posets import FinitePoset

from . import affinemonoid 
from . import monomialideal
from . import properpair
from . import _zsolve_interface
from . import _string_interface


#### Methods on covers of pairs.
"""
    Here, cover means a dictionary whose keys are faces of an affine semigroup,
    and whose values are pairs lie in the affine semigroup.

    cover
        _cover_difference(cover A, cover B):
            return difference A-B. In this case, pairs are the same if they have the same monomials and faces; in other words, we don't care their ambient ideal.
            # Moreover, pairs in cover B was preserved; not in A.
        _cover_intersection(cover A, cover B):
            return intersection of covers A and B.
"""

def _cover_difference(cover,new_cover):
    difference_cover = {}
    empty_key=[]
    for key in cover:
        if key in new_cover:
            # Change hash string not including ideal information
            hash_cover =[pair._hashstr() for pair in cover[key]]
            hash_new_cover=[pair._hashstr() for pair in new_cover[key]]
            if set(hash_cover) != set(hash_new_cover):
                diff_list = list(set(hash_cover).difference(set(hash_new_cover)))
                ind_list = [hash_cover.index(hashstr) for hashstr in diff_list]
                ind_list.sort()
                difference_cover[key] = [cover[key][ind] for ind in ind_list]
                if difference_cover[key] == []:
                    empty_key.append(key)
        else:
            difference_cover[key] = cover[key]
    # Check and delete empty.
    for key in empty_key:
        del difference_cover[key]
    return difference_cover
def _cover_intersection(dictone,dicttwo):
    """Given two dictionary with standard pairs, return their intersection."""
    if not (isinstance(dictone,dict) and isinstance(dicttwo,dict)):
        raise ValueError("[Error]:One of the argument is not dictionary.")
    inter = list(dictone.keys() & dicttwo.keys())
    returning_dict = {}
    for key in inter:
        hashs_of_dictone = set([pair._hashstr() for pair in dictone[key]])
        hashs_of_dicttwo = set([pair._hashstr() for pair in dicttwo[key]])
        intersection_hashs = list(hashs_of_dictone & hashs_of_dicttwo)
        list_of_pairs = _unique_pairs(dictone[key]+dicttwo[key])
        returning_pairs =[]
        for pair in list_of_pairs:
            if (pair.hashstring in intersection_hashs):
                returning_pairs.append(pair)
        if returning_pairs:
            returning_dict[key]=returning_pairs
    return returning_dict

#### Methods for finding overlap classes from a cover of pairs.
"""
    Here, ovCover means a dictionary whose keys are faces of an affine semigroup,
    and whose values are overlap class (i.e. list of ProperPairs overlapping each other)
     lie in the affine semigroup.

    list (of pairs)
        _overlap_pairs(ProperPair (a,F), list L(of ProperPairs)):
            return pairs in L overlapping with (a,F).
    ovCover
        _find_overlap_classes(cover A):
            return a ovCover based on the given cover A.

        _find_maximal_overlap_classes_in_cover(ovCover A)
            return new ovCover, containing the maximal overlap classes
            w.r.t. divisibility.
"""
def _overlap_pairs(a_pair,list_of_pairs):
    """Given a pair P and list of pairs L, return a list of pairs in L which overlap P."""
    if not isinstance(a_pair, properpair.ProperPair):
        raise ValueError('[Error]: Given the first instance is not a proper pair')
    if not isinstance(list_of_pairs,list):
        raise ValueError('[Error]: Given the second instance is not a list')
    if not list_of_pairs:
        #list of pairs is empty
        return([a_pair])
    for pair in list_of_pairs:
        if (not isinstance(pair, properpair.ProperPair)) or ((pair.ambient_ideal()).ambient_monoid() != (a_pair.ambient_ideal()).ambient_monoid()):
            raise ValueError('[Error]: Given list of pairs contain an element which is either not proper pair or a proper pair with distinct ambient monoid.')
    # If a_pair is pair with zero face, then return it; it cannot be overlap with other standard pairs.
    if (a_pair.face() == ()):
        return [a_pair]
    overlap_classes=[a_pair] #Initialize the returnign list.
    for pair in list_of_pairs:
        # Compare two pairs.
        # Here we care about whether two pairs has overlapped part or not.
        result = ((a_pair.ambient_ideal()).ambient_monoid())._intersection_of_pairs(a_pair.monomial(),a_pair.face(), pair.monomial(), pair.face())
        if result:
            # Result is not empty;
            overlap_classes.append(pair)
    return overlap_classes 

def _find_overlap_classes(cover):
    """ This is the implementation of overlap classes; used only 
    for internally in monomialideal.MonomialIdeal class object or standard cover object
    Do not use otherwise, since it lacks error checking code."""
    new_cover ={}
    for key, value in cover.items():
        if len(value) <= 1:
            new_cover[key]=[value]
        else:
            list_of_overlaps =[]
            # Copy attributes to delete it freely.
            new_list = [ind for ind in value]
            while len(new_list)>0:
                fixed_pair = new_list[0]
                del new_list[0]
                overlap_class =_unique_pairs(_overlap_pairs(fixed_pair,new_list))
                for pair in overlap_class:
                    if pair in new_list:
                        new_list.remove(pair)
                list_of_overlaps.append(overlap_class)
            new_cover[key]=list_of_overlaps
    return new_cover

def _find_maximal_overlap_classes_in_cover(classes):
    """ This is an implementation of maximalOverlap Classes over a given cover; 
    used only for _stdcover_to_ideal."""
    dict_classes = {}
    for key,value in classes.items():
        max_list = [];
        # It may be poset; that's why we need O(n^2) algorithm
        for pair_a in value:
            parity = True
            for pair_b in value:
                #Check whether pair_b is divisible by pair_a
                #In case they are not the same.
                if pair_a[0] != pair_b[0]:
                    if properpair.div_pairs(pair_a[0],pair_b[0]):
                        parity = False
                        break;
            if parity == True:
                max_list.append(pair_a)
        dict_classes[key]=max_list
    return dict_classes
   

#### Theorem 4.1
"""
    list (of ProperPair objects)
        pair_difference(pair (a,F), pair (b,G)):
            1) face F should be face of G.
            2) return list of proper pairs (c,F') such that
                the union of corresponding submonoids c+NF' is equal to
                (a+NF)-(b+NG)
    list (of lists of integers)
        _std_pairs_poly_ring(int n, matrix gens):
            1) n = gens.shape[1], gens is a matrix representing generator of ideal in the polynomial ring.
            2) return its standard pairs as integer notation.
            This method is used for internal communication between Macaulay2
            and other classes.

"""
def pair_difference(pair_a, pair_b):
    r"""
    Given two pairs ``pair_a`` and ``pair_b`` representing :math:`(a,F)` and :math:`(b,G)` such that
    :math:`F \subseteq G`, return a collection :math:`C` such that :math:`(a+\mathbb{N}F) \setminus (b+\mathbb{N}G) = \bigcup_{(c,H) \in C}(c+\mathbb{N}H)`.
    This is an implementation of Theorem 4.1 in https://arxiv.org/abs/2005.10968
    
    INPUT:

    - ``pair_a``, ``pair_b`` -- ``ProperPair`` objects. ``pair_a.face()`` must contain ``pair_b.face()`` as a subtuple.

    OUTPUT:

    - A ``dictionary`` object whose keys are face, and whose values are ``ProperPair`` objects having the same face with its key.

    EXAMPLE::

        sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair, pair_difference
        sage: Q = AffineMonoid(matrix(ZZ, [[2,0,1],[0,1,1]])) 
        sage: I = MonomialIdeal(matrix(ZZ,0),Q)
        sage: P= ProperPair(matrix(ZZ,[[0],[0]]), (0,1,2), I ) 
        sage: T= ProperPair(matrix(ZZ,[[0],[2]]), (0,1,2), I ) 
        sage: diff_PT = pair_difference(P,T)
        sage: diff_PT.keys()
        dict_keys([(0,)])
        sage: sorted(diff_PT[(0,)],key=str)
        [([[0], [0]]^T,[[2], [0]]),
         ([[0], [1]]^T,[[2], [0]]),
         ([[1], [1]]^T,[[2], [0]]),
         ([[1], [2]]^T,[[2], [0]])]
    """
    
    if (not isinstance(pair_a,properpair.ProperPair)) or (not isinstance(pair_b,properpair.ProperPair)):
        raise ValueError("[Error]: At least one of the argument is not a proper pair.")
    if (pair_a.ambient_ideal()).ambient_monoid() != (pair_b.ambient_ideal()).ambient_monoid():
        raise ValueError("[Error]: Two pairs have different ambient monoid")
        # Notes that ambient ideal may be distinct.
    if  (not (((pair_a.ambient_ideal()).ambient_monoid()).face_lattice()).is_lequal(pair_a.face(),pair_b.face())):
        raise ValueError("[Error]: The face of 1st argument is not less than or equal to that of 2nd argument.")
    if (((pair_a.ambient_ideal()).ambient_monoid()).is_pointed() == False):
        raise ValueError("[Error]: An affine monoid is not pointed, hence it does not work.")
    ### Before going on, if faces of both pairs are zero,
    if (pair_a.face() == ()):
        #Check whether pair_a.monomial is the member of pair_b or not
        if np.array(pair_b.is_element(pair_a.monomial())).astype('int64').size >0 :
            return {}
        else:
            # return pair_a
            return {pair_a.face():[pair_a]}

    # If pair_a = (a,F), pair_b=(b,G), with nonzero F,G, then solve
    # a+u*F = b+v*G for u and v
    # as a form
    # [F;-G][u;v] = b-a
    minsol = ((pair_a.ambient_ideal()).ambient_monoid())._intersection_of_pairs(pair_a.monomial(), pair_a.face(), pair_b.monomial(), pair_b.face())
    #Determine whether transpose it or not.
    shape_of_face_a=(((pair_a.ambient_ideal()).ambient_monoid()).face(pair_a.face())).shape
    #Reduce minsol = [F;-G] to [F] if F is nonzero
    if shape_of_face_a[1] == 0:
        # In case of pair_a.face() == (),
        if minsol:
            #minsol exists => pair_a.monomial is contained in pair_b
            return {}
        else:
            return {pair_a.face():[pair_a]}
    else:
        minsol=minsol[range(minsol.nrows()),range(shape_of_face_a[1])]
    if minsol.nrows()>0: # If minimal solution exists
        # Get a standard pair of J in k[N^pair_a.face]
        stdpair_j=_std_pairs_poly_ring(minsol.ncols(), minsol)
        # Result, a cover in paper, is a dictionary whose key is a face of ambient monoid
        # and whose value is set of pairs sharing faces.
        result={face:[] for face in (((pair_a.ambient_ideal()).ambient_monoid()).face_lattice()).list()}
        for idx in stdpair_j:
            # generate b+ Gu in the paper
            u= idx[0][np.newaxis].T
            g=((pair_a.ambient_ideal()).ambient_monoid()).face(pair_a.face())
            newmonomial=pair_a.monomial()+np.matmul(g,u)
            # generate N{a_i: i in sigma} in the paper
            newface=tuple(sorted(tuple([(pair_a.face())[index] for index in idx[1]])))
            result[newface].append(properpair.ProperPair(newmonomial, newface, pair_a.ambient_ideal(), True))
        for face in (((pair_a.ambient_ideal()).ambient_monoid()).face_lattice()).list():
            if not result[face]:
                del result[face]
        # Change lists in result on set
        return {key:_unique_pairs(val) for key,val in result.items()}
    else:
        return({pair_a.face():[pair_a]})

def _std_pairs_poly_ring(numvar, genset):
    """ Given a generating set of a polynomial ring, return standard pairs using Macaulay2
    This is a subroutine of Theorem 4.1 in https://arxiv.org/abs/2005.10968
        Arguments:
            int numvar : Number of attributes of the polynomial ring
            numpy (2d)Array or sage matrix 
                genset : A matrix whose row represent an exponent notation of the given polynomial ring.
        Return
            list(list) int_result: List of standard pairs.
    """
    if not isinstance(numvar,int):
        raise ValueError("1st argument is not an integer.")
    #Generate Poly ring with numvar attributes
    gens_poly_ring='R=ZZ[vars (0..'+str(numvar-1)+')]'
    macaulay2(gens_poly_ring)
    # Generate ideal whose generator corresponds to a row of generating set.
    gens_ideal='I=monomialIdeal(';
    for row in genset:
        eq=''
        for idx in range(numvar):
            eq=eq+'R_'+str(idx)+'^'+str(row[idx])+'*'
        # Delete last +
        eq = eq[:-1]
        gens_ideal=gens_ideal + eq+','
    gens_ideal=gens_ideal[:-1] + ')'
    macaulay2(gens_ideal)
    # Find standard pair of ideal.
    macaulay2('L=standardPairs I')
    # Translate lists in Macaulay2 for python.
    result=list(macaulay2('T=apply(L, i -> {(exponents i_0)_0, apply(i_1, j -> index j)})'))
    result = [list(i) for i in result]
    result = [ [list(i[0]), list(i[1])] for i in result]
    int_result=[]
    for i in result:
        temp=[];
        for j in i[0]:
            temp.append(int(j))
        temp = np.array(temp).astype('int64')
        tempt = [];
        for j in i[1]:
            tempt.append(int(j))
        int_result.append([temp,tempt])
    return int_result

#### Lemma 4.2
"""
    sage.matrix
        _minimal_holes(vector a, tuple face, affinemonoid.AffineMonoid Q):
            return the minimal solutoin of the intersection between (a+RF) and Q.

"""
def _minimal_holes(col_vector, face, affine_monoid):
    """ Given a vector and face, find the minimal element of intersection of (vector+ RR* face) and affine_monoid 
    This is an implementation of Lemma 4.2 in https://arxiv.org/abs/2005.10968
        Arguments:
            numpy (2D array) vector: a column integer vector
            affinemonoid.AffineMonoid affine_monoid: an affine monoid having the vector as an element.
            tuple face : a face of the affine_monoid
        
        Return:
            sage.matrix
                result[0]: Minimal solution of the intersection between (vector+ RR* face) and affine_monoid.
            
    
    """
    #Genset may numpy 2d array or sage 2d matrix with integer.
    #Preprocess
    if len(col_vector.shape) !=2:
        raise ValueError("1st argument must be a column vector of 2D array.")
    if col_vector.shape[1] !=1:
        raise ValueError("1st argument must be a column vector of 2D array.")
    if not isinstance(affine_monoid,affinemonoid.AffineMonoid):
        raise ValueError("3rd argument is not a class of AffineMonoid.''")
    if not( face in (affine_monoid.face_lattice()).list()):
        raise ValueError("2nd argument is not a face of 3rd argument.")
    
    hvectors = (affine_monoid.integral_support_vectors())[face]
    hvector_times_gens = np.matmul(hvectors, affine_monoid.gens())
    # Make a matrix B=[A;A]
    #RHS : [hvectors times given vector;hvectors times given col_vector]
    return _zsolve_interface._zsolve_exact_solver(hvector_times_gens, np.matmul(hvectors,col_vector))

#### Proposition 4.4
"""
    cover
        _cover_to_std_pairs(cover A, MonomialIdeal I, loop_max=100):
            Given a cover A of the ideal I, generate new cover consisting of 
            standard pairs of I.
            (loop_max is needed for preventing the infinite loop.)
        _czero_to_cone(cover A, MonomialIdeal I):
        _cone_to_ctwo(cover A, MonomialIdeal I):
            Each methods do refinement C_0->C_1 and C_1 -> C_2 
            described in Proposition 4.4.
        _max_standard_pairs(cover A, MonomialIdeal I):
            Given a cover A of the ideal I, return cover consisting of 
            standard pairs of I from A. (Notes that it do not refine the pairs itself.)
"""
def _cover_to_std_pairs(cover, ambient_ideal):
    """This is an implementation of Proposition 4.4 of https://arxiv.org/abs/2005.10968"""
    new_cover={}
    loop_num=0
    parity=0
    while (parity == 0):
        new_cover=_czero_to_cone(cover)
        new_cover=_cone_to_ctwo(new_cover,ambient_ideal)
        # Now compare two covers.
        loop_num = loop_num+1
        if len(_cover_difference(cover,new_cover))==0:
            parity=1
        else:
            cover = new_cover
        if (loop_num % 1000) == 0:
            warnings.warn("[Warning]: # of loops exceeds "+str(loop_num), RuntimeWarning)
    #If stabilizes, then take the maximal elements.
    max_cover = _max_standard_pairs(new_cover,ambient_ideal)
    return max_cover

#C_0 to C_1
def _czero_to_cone(cover):
    """This is the 1st subroutine of Proposition 4.4 of https://arxiv.org/abs/2005.10968"""
    new_cover={}
    for face, pair_set in cover.items():
        pair_list = list(pair_set)
        new_list=[]
        for pair in pair_list:
            minsol=_minimal_holes(pair.monomial(), pair.face(), ((pair.ambient_ideal()).ambient_monoid()))
            result=np.matmul(((pair.ambient_ideal()).ambient_monoid()).gens(), minsol.transpose().numpy()).T                
            # Use monomials in minsol which divides the original monomial
            new_result=[ col_vector[np.newaxis].T for col_vector in result if ((pair.ambient_ideal()).ambient_monoid()).is_element(pair.monomial()-col_vector[np.newaxis].T).nrows() >0 ]
            result=new_result
            for column in result:
                if not np.array_equal( column, pair.monomial()):
                    new_list.append(properpair.ProperPair(column, face, pair.ambient_ideal(),True))
                else:
                    new_list.append(pair)
        new_cover.update({face:_unique_pairs(new_list)})
    return new_cover
#C_1 to C_2
def _cone_to_ctwo(new_cover, ambient_ideal):
    """This is the 2nd subroutine of Proposition 4.4 of https://arxiv.org/abs/2005.10968"""
    renew_cover={}
    # This is new pair list, which we will store for return
    new_pair_list=[]
    used_faces = set([])
    # To avoid repeated calculation, just sort all items in the cover by monomials and its related face.
    monomial_face_dict={}
    for face, pair_list in new_cover.items():
        for pair in pair_list:
            if _string_interface._np2d_to_string(pair.monomial()) not in monomial_face_dict.keys():
                monomial_face_dict[_string_interface._np2d_to_string(pair.monomial())]=[face]
            else:
                monomial_face_dict[_string_interface._np2d_to_string(pair.monomial())].append(face)
    # Now use that monomial_face_dict to make it much more conveniently.
    #print("C1->C2, given new_cover:",new_cover)
    #print("C1->C2, monomial_face_dict:",monomial_face_dict)

    for monomial, list_of_faces in monomial_face_dict.items():
        #Find complement of downsets.
        complement = _complement_of_down_set(list_of_faces, ambient_ideal.ambient_monoid())
        # Sort it by length of tuples
        complement.sort(key=len, reverse = True)
        # Store possible faces for avoiding duplicates
        possible_faces=[]
        # Now try to generate proper pair of the ambient_ideal for each faces in the complement.
        for face in complement:
            # Check that whether the given face is contained in a face which are possible to make properpair
            is_face_contained_in_possible_face=False
            for item in possible_faces:
                if set(face).issubset(set(item)):
                    is_face_contained_in_possible_face=True
            # If the given face is not contained in one of possible faces,
            if (not is_face_contained_in_possible_face):
                try:
                    temp_pair=properpair.ProperPair(_string_interface._string_to_np2d(monomial),face,ambient_ideal)
                    possible_faces.append(face)
                    del temp_pair
                except Warning:
                    pass
        # Add original faces in the possible_faces
        possible_faces = possible_faces+list_of_faces
        # Find maximal elements of possible faces as a subposet.
        maximal_faces=(((ambient_ideal.ambient_monoid()).face_lattice()).subposet(possible_faces)).maximal_elements()
        used_faces = used_faces.union(set(maximal_faces))
        # Now save it in new_pair_list
        for face in maximal_faces:
            new_pair_list.append(properpair.ProperPair(_string_interface._string_to_np2d(monomial), face, ambient_ideal))
            
    # So now we generate all new_pair_lists. 
    # Sort it as a form of cover dictionary.
    for item in new_pair_list:
        if item.face() not in renew_cover.keys():
            renew_cover.update({item.face(): [item]})
        else:
            renew_cover[item.face()].append(item)
    #Lastly delete repeated ones.
    returning_cover = {key:_unique_pairs(val) for key, val in renew_cover.items()}
    #print("C1->C2, returning_cover:",returning_cover)
    return returning_cover

def _max_standard_pairs(cover, ambient_ideal):
    """ This is 3rd subroutine of Proposition 4.4 of https://arxiv.org/abs/2005.10968
        Argument:
            dict cover {face: set(ProperPair...)}
                A cover of ProperPairs.
            MonomialIdeal ambient_ideal
                An ambient ideal which has all proper pairs in the cover as proper pair.
        
        Return
            dict new_dict {face:list(ProperPair...)}
                A dictionary of all standard pairs in the given cover.
        
    """
    if not isinstance(cover, dict):
        raise ValueError("[Error]: The instance is not a type of dictionary")
    new_dict={}
    for face, pair_set in cover.items():
        pair_list = list(pair_set)
        if not (face in ((ambient_ideal.ambient_monoid()).face_lattice()).list()):
            raise ValueError("[Error]: Given face in a dictionary is not a face of ambient monoid")
        for i in pair_list:
            if i.ambient_ideal() != ambient_ideal:
                raise ValueError("[Error]: \n",i ,"\n is not a pair of the given ambient ideal.")
        hvectors=(ambient_ideal.ambient_monoid()).integral_support_vectors()[face]
        h_coordinates=[np.matmul(hvectors,i.monomial()).flatten() for i in pair_list]
        # Find unique H coordinates.
        unique_h_coord=[ np.array(i).astype('int64')  for i in set(map(tuple, h_coordinates))]
        # List for non-standard pairs
        updated_list=[];
        for arr in unique_h_coord:
            #Find a duplicaed list
            tf_list = []
            # Now given truth check of H-coordinates
            for i in (h_coordinates==arr):
                #If all H-coordinates do not agree with given coordinate,
                if i.all()==False:
                    tf_list.append(0)
                else:
                    tf_list.append(1)
            if sum(tf_list)>1:
                #If there is at least two elements whose H-coordinates are the same,
                #Generate indices of those elements,
                idx=[i for i, x in enumerate(tf_list) if x == 1]
                # Generate new list of these pairs with the same H-coordinates.
                specific_pairs = [pair_list[i] for i in idx]
                #Here we generate a digraph whose vertices are index of specific_pairs, and
                #whose edges are order relation.
                list_of_vertices = [i for i in range(len(specific_pairs))]
                list_of_edges = []
                #Storing attribute for pairs within the same overlap classes.
                overlap_classes=[] # We order it by jdx.
                #Generate combination of elements.
                for jdx in range(len(specific_pairs)-1):
                    for kdx in range(len(specific_pairs)-jdx-1):
                        #Compare two sets.
                        # Notes that if a pair_a contains pair_b, then _intersection_of_pairs(pair_a,pair_b)
                        # returns [nonzero half part, zero half part]. (Likewise, if pair_b contains pair_b, then
                        # it returns [zero half part, nonzero half part].)
                        n_kdx = kdx+jdx+1
                        #Compare two pairs
                        first = specific_pairs[jdx].is_element(specific_pairs[n_kdx].monomial())
                        second = specific_pairs[n_kdx].is_element(specific_pairs[jdx].monomial())

                        if (first.numpy().size > 0) and (second.numpy().size >0): 
                            check_whether_stored=False
                            for ldx in range(len(overlap_classes)):
                                # If there is already a list of overlapped class which contains one of the pair, append another.
                                if (jdx in overlap_classes[ldx]):
                                    overlap_classes[ldx].append(n_kdx)
                                    check_whether_stored=True
                                if (n_kdx in overlap_classes[ldx]):
                                    overlap_classes[ldx].append(jdx)
                                    check_whether_stored=True
                            if check_whether_stored == False:
                                # If the above process do not store these overlapped pairs, make new overlapped class.
                                overlap_classes.append([jdx, n_kdx]) #Here order matters, for deleting cycles.
                                #Comment: if we have j_0 <j_1< cdots < j_k pairs in the same overlap classes,
                                # then loop will be given as (j_0, j_1), cdots (j_0, j_1), (j_1,j_2), cdots
                                # Thus, j_0 should be in the first element of a sublist of overlap_classes.
                                # This fact will be used to contracting digraph.
                                list_of_edges.append((n_kdx,jdx))
                                list_of_edges.append((jdx,n_kdx))
                        elif (first.numpy().size > 0) and (second.numpy().size == 0):
                            #If first matrix is nonempty, i.e., pair[jdx] contains pair[n_kdx] i.e., pair[n_kdx] divides pair[jdx]
                            list_of_edges.append((n_kdx,jdx))
                        elif (first.numpy().size == 0) and (second.numpy().size > 0):
                            #If second matrix is nonempty, i.e., pair[n_kdx] contains pair[jdx] i.e., pair[jdx] divides pair[n_kjdx]
                            list_of_edges.append((jdx,n_kdx))
                        else:
                            # If two pairs are the in the same H-coordinates but not overlapped,
                            # Then we do not need to make an edge between them
                            pass
                #Using vertices and edges, generate sage digraph
                g = DiGraph([list_of_vertices,list_of_edges])
                # Contract digraph by overlap_classes
                for sublist in overlap_classes:
                    for ldx in range(len(sublist)-1):
                        g.contract_edge((sublist[0],sublist[ldx+1]))
                # From above, we just delete overlap classes, so we do not have error on temp_poset.
                temp_poset = FinitePoset(g)
                max_elements = temp_poset.maximal_elements()
                # Now we add deleted overlap pairs if it is a maximal element.
                overlapped_max_elements=[]
                ldx = 0
                while ldx < len(overlap_classes):
                    udx =0
                    while udx <len(max_elements):
                        if max_elements[udx] in overlap_classes[ldx]:
                            overlapped_max_elements = overlapped_max_elements +overlap_classes[ldx]
                            ldx=ldx+1
                            udx=0
                            
                        udx=udx+1
                    ldx = ldx+1
                    udx=0
                updated_list = updated_list + [specific_pairs[i] for i in overlapped_max_elements] + [specific_pairs[i] for i in max_elements]
            else:
                updated_list.append( pair_list[tf_list.index(1)] )
        new_dict.update({face:updated_list})
    # Last, delete repeated dictionary
    renew_dict = { key:_unique_pairs(val) for key, val in new_dict.items()}
    return renew_dict

#### Theorem 4.5
"""
    cover
        _standard_pairs(MonomialIdeal I, verbose = True):
            return standard pairs of I as a cover. If verbose != True, then show its progress on how much calculation was completed.
"""
def _standard_pairs(an_ideal,verbose=True):
    """Return the standard cover of the given ideal.
    This is an implementation of Theorem 4.5 in https://arxiv.org/abs/2005.10968
        Argument:
            MonomialIdeal an_ideal
                An ideal of an affine monoid.
        
        Return
            dict new_dict {face:list(ProperPair...)}
                A dictionary of standard cover.
        
    """
    if an_ideal._is_std_cover_calculated():
        return an_ideal.standard_cover()
    if not isinstance(an_ideal, monomialideal.MonomialIdeal):
        raise ValueError("[Error]: The parameter is not a type of MonomialIdeal")
    if ((an_ideal.ambient_monoid()).is_pointed() == False):
        raise ValueError("[Error]: An ambient affine monoid is not pointed.")
    #if _is_permuation_matrix(an_ideal.ambient_monoid.mingens):
        #This is the case when an affine semigroup ring is a polynomial ring.
    #    transp_generator = matrix(an_ideal.gens).transpose()
    #    stdpair_j=_std_pairs_poly_ring(transp_generator.ncols(),transp_generator)
    #    # All generator matrix of a polynomial ring should contain a standard basis.
        # Figure out which number of columns are standard pairs.
    #    attributes_list = [];
    #    temp_col = matrix(an_ideal.ambient_monoid.gens).columns()
    #    for col in temp_col:
    #        length =0
    #        for col_idx in col:
    #            length = length +col_idx*col_idx
    #        if length == 1:
    #            attributes_list.append(list(col).index(1))
    #        else:
    #            attributes_list.append(-1)
        
    #    # Result, a cover in paper, is a dictionary whose key is a face of ambient monoid
    #    # and whose value is set of pairs sharing faces.
    #    result={face:[] for face in an_ideal.ambient_monoid.face_lattice.list()}
    #    for idx in stdpair_j:
    #        # generate b+ Gu in the paper
    #        u= idx[0][np.newaxis].T
    #        face = tuple(sorted([attributes_list.index(poly_var) for poly_var in idx[1]]))
    #        result[face].append(ProperPair(u,face,an_ideal, True))
    #    for face in an_ideal.ambient_monoid.face_lattice.list():
    #        if not result[face]:
    #            del result[face]
    #    # Change lists in result on set
    #    return {key:_unique_pairs(val) for key,val in result.items()}
    
    zero_ideal = monomialideal.MonomialIdeal( np.array([], dtype='int64'), an_ideal.ambient_monoid())
    # For the first element, find initial standard cover for principal ideal.
    pair_c = properpair.ProperPair((an_ideal.gens())[:,[0]] * 0, (an_ideal.ambient_monoid()).index_of_face((an_ideal.ambient_monoid()).gens()), zero_ideal)
    pair_d = properpair.ProperPair((an_ideal.gens())[:,[0]], (an_ideal.ambient_monoid()).index_of_face((an_ideal.ambient_monoid()).gens()), zero_ideal)
    cover = pair_difference(pair_c,pair_d)
    if an_ideal.is_principal():
        return _cover_to_std_pairs(cover, an_ideal)
    else:
        if (verbose != True):
            print("Calculate the standard cover of an ideal")
            print("It takes a few minutes, depending on the system.")
        count_column=1
        number_of_columns = (an_ideal.gens()).shape[1]
        if (verbose != True):
            print("Cover for", str(count_column) ," generator was calculated. ", str(number_of_columns-1)," generators remaining. ")
        for idx in range((an_ideal.gens()).shape[1]-1): # For loop by # of columns
            # Since we already use 0th column, idx will be used from 1 to (n-1), i.e., idx+1
            new_pair = properpair.ProperPair((an_ideal.gens())[:,[idx+1]], (an_ideal.ambient_monoid()).index_of_face((an_ideal.ambient_monoid()).gens()), zero_ideal)
            temp_ideal = monomialideal.MonomialIdeal(an_ideal.gens()[:,0:(idx+2)], an_ideal.ambient_monoid())
            new_cover={}
            for list_pair in cover.values():
                for pair in list_pair:
                    pair_diff = pair_difference(pair,new_pair)
                    for key, val in pair_diff.items():
                        if key in new_cover.keys():
                            new_cover[key] = _unique_pairs(new_cover[key] + val)
                        else:
                            new_cover[key] = _unique_pairs(val)
            # Now update our cover by prop 4.4
            
            cover = _cover_to_std_pairs(new_cover, temp_ideal)
            count_column = count_column+1
            number_of_columns = number_of_columns-1
            if (verbose != True):
                print("Cover for", str(count_column) ," generators was calculated. ", str(number_of_columns-1)," generators remaining. ")
        return cover

#### Theorem 4.8
"""
    MonomialIdeal
        _stdcover_to_ideal(cover A, affinemonoid.AffineMonoid Q):
            given a standard cover A whose pairs are over Q,
            return the ideal whose standard cover is A.
    [bool, MonomialIdeal]
        _initiate_stdcover_to_ideal(cover A, affinemonoid.AffineMonoid Q)
            1st routine of the _stdcover_to_ideal(...), which approximate
            the generators first time.
            return False when it is successful, True when the solution is zero ideal.
        _coordination_cover_to_ideal(cover A, affinemonoid.AffineMonoid Q, MonomialIdeal candidate)
            2nd routine of the _stdcover_to_ideal(...). Based on the approximation
            in the previous method, calculate the difference of covers in the
            candidate, and add all monomials which are not in the given cover
            to return the better estimation of the ideal.
            Boolean value indicates whether the returned ideal is exactly what we found
            or not.
    bool
        _is_cover_valid(dict A):
            check whether dictionary A is a valid cover or not.


"""
def _stdcover_to_ideal(cover,ambient_monoid):
    """ Given a standard cover of an ideal, return the ideal as an 'MonomialIdeal' object
    This is an implementation of Theorem 4.8 in https://arxiv.org/abs/2005.10968
        Arguments:
            dictionary 
                cover: a dictionary which is a standard cover of some unknown ideal.

                affinemonoid.AffineMonoid ambient_monoid: an affine monoid where the cover lives.
        
        Return:
            MonomialIdeal
                an_ideal: an ideal having cover as its standard cover.
           """
    result = _initiate_stdcover_to_ideal(cover,ambient_monoid)
    while (result[0] == False):
        result = _coordination_cover_to_ideal(cover, ambient_monoid,result[1])
    return result[1]

def _is_cover_valid(cover, ambient_monoid):
    """Check the validity of a cover of proper pairs in some ambient monoid
    This is a subroutine of Theorem 4.8 for avoiing error."""
    if not isinstance(cover,dict):
        raise ValueError("[Error]: the 1st argument is not a dictionary.")
    if not isinstance(ambient_monoid,affinemonoid.AffineMonoid):
        raise ValueError("[Error]: the 2nd argument is not an AffineMonoid object.")
    for key, value in cover.items():
        if key not in (ambient_monoid.face_lattice()).list():
            raise ValueError("[Error]: One of the key is not in the face lattice")
        for pair in value:
            if not isinstance(pair,properpair.ProperPair):
                raise ValueError("[Error]: One of the pair is not the ProperPair object")
    return True

def _tracing_of_overlap_class(ov_class,face, level, unusable_rays):
    r"""
    Given an overlap class and level as intege, find all elements in overlap classes of pairs which can be obtained by
    adding ``level`` many consecutive elements in the given face of the overlap class to its generators.
    Notes that len(unusable_rays) should be the same as len(ov_class).
    """
    if ov_class == []:
        return []

    for pair in ov_class:
        if not isinstance(pair,properpair.ProperPair):
            raise ValueError("[Error]: One of the pair is not the ProperPair object")
    if not isinstance(face, tuple):
        raise ValueError("[Error]: One of the pair is not the ProperPair object")
    # Change monomial as 1D array object.
    gens_of_pairs = [pair.monomial().sum(axis=1) for pair in ov_class]
    base_monoid = (ov_class[0].ambient_ideal()).ambient_monoid()
    # Choose mingens index of faces
    result = []
    for kdx,item in enumerate(gens_of_pairs):
        #Now we need to construct it using unusable_rays
        mingens = [tuple(item) for item in base_monoid.mingens().T]
        min_face = [mingens.index(tuple(base_monoid.gens()[:,jtem])) for jtem in face if (tuple(base_monoid.gens()[:,jtem]) in mingens) and (jtem not in unusable_rays[kdx])]
        min_face = list(set(min_face))
        min_face.sort()
        min_face = tuple(min_face)
        mingens = np.concatenate([np.array(item)[np.newaxis].T for item in mingens],axis=1)

        # Find all possible combinations of rays in faces with repetition.
        possible_ways = list(itertools.product(min_face,repeat=level))
        # Change those combinations into sum of columns represented by combinations.
        vectors = [mingens[:, item].sum(axis=1) for item in possible_ways]
        del possible_ways
        for jtem in vectors:
            result.append( tuple(item+jtem))
        del vectors
    result = list(set(result))
    result = [np.array(item)[np.newaxis].T for item in result]
    return result

def _initiate_stdcover_to_ideal(cover,ambient_monoid):
    r"""
    This is 1st subroutine of Theorem 4.8 in https://arxiv.org/abs/2005.10968
    """
    if not _is_cover_valid(cover, ambient_monoid):
        raise ValueError("[Error]: Given cover is not a valid cover of the given monoid.")    
    if ambient_monoid.face_lattice().top() in cover:
        return [True, monomialideal.MonomialIdeal(matrix(ZZ,0),ambient_monoid)]
    max_overlap_classes =_find_maximal_overlap_classes_in_cover(_find_overlap_classes(cover))
    # All generators will be saved in final_gens
    num_cols_of_monoid_mingens = ambient_monoid.mingens().shape[1]
    final_gens=[]
    if () in cover:
        candid_of_gens_from_zero_face_pairs=[]
        for zero_pair in cover[()]:
            for ind in range(num_cols_of_monoid_mingens):
                candid_of_gens_from_zero_face_pairs.append(zero_pair.monomial()+ambient_monoid.mingens()[:,(ind,)])
        # Now make all items in candidates unique
        candid_of_gens_from_zero_face_pairs = list(set([tuple(item[:,0]) for item in candid_of_gens_from_zero_face_pairs]))
        candid_of_gens_from_zero_face_pairs = [np.array(item,dtype='int64')[np.newaxis].T for item in candid_of_gens_from_zero_face_pairs]
        # Then check whether it is in other standard covers.
        for monomial in candid_of_gens_from_zero_face_pairs:
            is_monomial_in_cover = False
            for face,list_of_pairs in cover:
                for pair in list_of_pairs:
                    if np.array(pair.is_element(monomial)).size>0:
                        is_monomial_in_cover = True
                        break
                if is_monomial_in_cover == True:
                    break
            if is_monomial_in_cover == False:
                final_gens.append(monomial)
    single_list_of_pairs = functools.reduce(lambda x,y: x+y, list(cover.values()))
    single_list_of_pairs_in_max = functools.reduce(lambda x,y: x+y,(functools.reduce(lambda x,y: x+y, list(max_overlap_classes.values()))))
    # First of all, find combination
    pairs_of_the_list = list(itertools.combinations(range(len(single_list_of_pairs)), r=2))
    # Save all combinations such that representing pairs have nonempty intersection.

    tuples_of_intersection =[item for item in pairs_of_the_list if np.array(ambient_monoid._intersection_of_pairs(single_list_of_pairs[item[0]].monomial(),single_list_of_pairs[item[0]].face(), single_list_of_pairs[item[1]].monomial(),single_list_of_pairs[item[1]].face())).size >0 ]
    tuples_of_intersection.sort()
    
    for pair in single_list_of_pairs_in_max:
        count = single_list_of_pairs.index(pair)
        # Skip the cases when face == () since it was done by above code.
        if pair.face() == ():
            continue
        # Also, choose pairs which already have nonempty intersections.
        pairs_for_comparing = [single_list_of_pairs[  list(set(item).difference(set([count])))[0] ] for item in tuples_of_intersection if count in item]
        tracing_monomials_in_pair = [pair.monomial()]
        candid_monomials=[]
        
        while len(tracing_monomials_in_pair)>0:
            # This start_leveling controls next tracing.
            start_leveling=1
            result_boolean=[]
            #unusable_rays keep which rays in pair.face() cannot be usable to tracing, 
            #since it overlap with some other_pair whos face also share the ray.
            unusable_rays=[]
            for candid in tracing_monomials_in_pair:
                parity_for_candid = False
                temp_rays=[]
                for other_pair in pairs_for_comparing:
                    if np.array(other_pair.is_element(candid)).size>0:
                        parity_for_candid=True
                        #Instead of break, check whether the other pair's face overlap with pair's face.
                        temp_rays=temp_rays+[jdx for jdx in other_pair.face()]
                temp_rays = list(set(list(temp_rays)))
                temp_rays.sort()
                unusable_rays.append(tuple(temp_rays))
                result_boolean.append(parity_for_candid)
            # Delete monomials which is in other standard pairs.
            result_index = [jdx for jdx,temp in enumerate(result_boolean) if temp == True]
            dead_monomials_as_pair =[]
            dead_monomials_unusable_rays=[]
            for jdx in sorted(result_index, reverse=True):
                dead_monomials_unusable_rays.append(unusable_rays[jdx])
                dead_monomials_as_pair.append(properpair.ProperPair(tracing_monomials_in_pair[jdx],pair.face(), pair.ambient_ideal()))
                del unusable_rays[jdx]
                del tracing_monomials_in_pair[jdx]
            # Recall that if tracing monomial is not in any of std pair,
            #then unusable_rays for that monomial is empty.
            # Now generate
            for item in tracing_monomials_in_pair:
                candid_monomials.append(item)
            tracing_monomials_in_pair = _tracing_of_overlap_class(dead_monomials_as_pair, pair.face(), start_leveling, dead_monomials_unusable_rays)
            while (len(candid_monomials)==0 and len(tracing_monomials_in_pair)==0):
                # This loop is needed since candid_monomials cannot be zero by construction of std pairs.
                start_leveling = start_leveling+1
                # To get monomials, delete the unusable_ray condition.
                relieved_unusable_rays = [() for item in dead_monomials_unusable_rays]
                tracing_monomials_in_pair = _tracing_of_overlap_class(dead_monomials_as_pair, pair.face(), start_leveling, relieved_unusable_rays)
            del dead_monomials_as_pair
            del result_boolean
            del result_index
        # Leave unique monomials
        candid_monomials = _string_interface._unique_np_arrays(candid_monomials)
        # To check redundant one, generate mingens
        mingens = [tuple(item) for item in ambient_monoid.mingens().T]
        monomials_in_face=[item for item in mingens if item in [tuple(item) for item in ambient_monoid.gens()[:,pair.face()].T ]]
        # Generate possible monomials in ideal.
        for candid in candid_monomials:
            for jdx in range(num_cols_of_monoid_mingens):
                if tuple(ambient_monoid.mingens()[:,jdx]) not in monomials_in_face:
                    final_gens.append(candid + ambient_monoid.mingens()[:,[jdx]])

        # To prevent memory leak,
        del pairs_for_comparing
        del candid_monomials
    # Again, leave unique monomials
    final_gens = _string_interface._unique_np_arrays(final_gens)
    if final_gens:
        candidate_ideal = monomialideal.MonomialIdeal(np.concatenate(final_gens,axis=1),ambient_monoid)
        return [False, candidate_ideal]
    else:
        #Return zero ideal
        return [True, monomialideal.MonomialIdeal(matrix(ZZ,0),ambient_monoid)]

def _coordination_cover_to_ideal(cover,ambient_monoid,candidate_ideal):
    """Coordination part; do not use it otherwise.
    This is 2nd subroutine of Theorem 4.8 in https://arxiv.org/abs/2005.10968"""
    num_cols_of_monoid_mingens = ambient_monoid.mingens().shape[1]
    # Generate pairs in new cover,which is not in original one.
    pairs_only_in_new_cover = _cover_difference(candidate_ideal.standard_cover(),cover)
    # If there is no difference, return it as true!
    if len(pairs_only_in_new_cover) == 0:
        return [True, candidate_ideal]
    # Otherwise, let's add monomials which are from pairs only in new cover.
    
    #Again,
    single_list_of_new_pairs = functools.reduce(lambda x,y: x+y, list(pairs_only_in_new_cover.values()))
    single_list_of_pairs = functools.reduce(lambda x,y: x+y, list(cover.values()))
    # Save all combinations such that representing pairs have nonempty intersection.
    pairs_of_pairs_of_above_two_vars = list(itertools.product(range(len(single_list_of_new_pairs)),range(len(single_list_of_pairs))))
    # Just find possible intersections.
    tuples_of_intersection =[item for item in pairs_of_pairs_of_above_two_vars if np.array(ambient_monoid._intersection_of_pairs(single_list_of_new_pairs[item[0]].monomial(),single_list_of_new_pairs[item[0]].face(), single_list_of_pairs[item[1]].monomial(),single_list_of_pairs[item[1]].face())).size >0 ]
    tuples_of_intersection.sort()
    #All generators will be saved in final_gens
    final_gens=[]
    for count, pair in enumerate(single_list_of_new_pairs):
        pairs_for_comparing = [single_list_of_pairs[item[1]] for item in tuples_of_intersection if count == item[0]]
        if (pair.face() == ()):
            if(len(pairs_for_comparing)==0): #If it is not an element of any original stdpairs
                final_gens.append(pair.monomial())
                continue
        tracing_monomials_in_pair = [pair.monomial()]
        candid_monomials=[]
        while len(tracing_monomials_in_pair)>0:
            # This start_leveling controls next tracing.
            start_leveling=1
            result_boolean=[]
            #unusable_rays keep which rays in pair.face() cannot be usable to tracing, 
            #since it overlap with some other_pair whos face also share the ray.
            unusable_rays=[]
            for candid in tracing_monomials_in_pair:
                parity_for_candid = False
                temp_rays=[]
                for other_pair in pairs_for_comparing:
                    if np.array(other_pair.is_element(candid)).size>0:
                        parity_for_candid=True
                        #Instead of break, check whether the other pair's face overlap with pair's face.
                        temp_rays=temp_rays+[jdx for jdx in other_pair.face()]
                temp_rays = list(set(list(temp_rays)))
                temp_rays.sort()
                unusable_rays.append(tuple(temp_rays))
                result_boolean.append(parity_for_candid)
            # Delete monomials which is in other standard pairs.
            result_index = [jdx for jdx,temp in enumerate(result_boolean) if temp == True]
            dead_monomials_as_pair =[]
            dead_monomials_unusable_rays=[]
            for jdx in sorted(result_index, reverse=True):
                dead_monomials_unusable_rays.append(unusable_rays[jdx])
                dead_monomials_as_pair.append(properpair.ProperPair(tracing_monomials_in_pair[jdx],pair.face(), pair.ambient_ideal()))
                del unusable_rays[jdx]
                del tracing_monomials_in_pair[jdx]
            # Recall that if tracing monomial is not in any of std pair,
            #then unusable_rays for that monomial is empty.
            # Now generate
            for item in tracing_monomials_in_pair:
                candid_monomials.append(item)
            del tracing_monomials_in_pair
            tracing_monomials_in_pair = _tracing_of_overlap_class(dead_monomials_as_pair, pair.face(), 1, dead_monomials_unusable_rays)
            while (len(candid_monomials)==0 and len(tracing_monomials_in_pair)==0):
                # This loop is needed since candid_monomials cannot be zero by construction of std pairs.
                start_leveling = start_leveling+1
                # To get monomials, delete the unusable_ray condition.
                relieved_unusable_rays = [() for item in dead_monomials_unusable_rays]
                tracing_monomials_in_pair = _tracing_of_overlap_class(dead_monomials_as_pair, pair.face(), start_leveling, relieved_unusable_rays)
            del dead_monomials_as_pair
            del dead_monomials_unusable_rays
            del result_boolean
            del result_index
        # Leave unique monomials
        candid_monomials = _string_interface._unique_np_arrays(candid_monomials)
        # Generate possible monomials in ideal.
        # Notes that now we don't need to add generators to candidate,
        # since candidates are already outside of the given cover.
        final_gens = final_gens+candid_monomials
        # To prevent memory leak,
        del pairs_for_comparing
        
    #Again, if final_gens are not empty,
    if final_gens:
        newly_gens = np.concatenate(final_gens,axis=1)
        candidate_ideal = monomialideal.MonomialIdeal(np.concatenate( [candidate_ideal.gens(), newly_gens], axis=1 ),ambient_monoid)
        return [False, candidate_ideal]
    else:
        #Return zero ideal
        warnings.warn("[Error]: Currently, we don't know why but two ideals not agree even if there is no added one. ")
        print("First of all, here is the original cover,")
        print(cover)
        print("And here is the standard cover by calculation")
        print(candidate_ideal.standard_cover())
        print("whose generators are", candidate_ideal)
        return [True, candidate_ideal]

def _complement_of_down_set(list_of_given_face, ambient_monoid):
    # Find list of faces
    face_lists = (ambient_monoid.face_lattice()).list()
    # Deal with trivial cases
    #print("list_of_given_face", list_of_given_face)
    if list_of_given_face == []:
        return face_lists
    face_lists.remove((-1,))    
    if list_of_given_face == [(-1,)]:
        return face_lists
    if (-1,) in list_of_given_face:
        list_of_given_face.remove((-1,))
        if () not in list_of_given_face:
            face_lists.remove(())
            
    new_lists =[]
    for item in face_lists:
        parity=False
        for given_face in list_of_given_face:
            if (set(given_face) >= set(item)): # If given face contains item as a subface,
                parity = True
        if parity == False:
            new_lists.append(item)
    return new_lists

#### Uniqueness and comparison of lists or covers.
"""
    list
        _unique_pairs(list of ProperPair objects):
            return list of (mathematically) unique ProperPair objects.
            Two pairs are the same if their monomial, their ambient ideal, and their faces are the same.

"""
def _unique_pairs(list_of_pairs):
    """ Added to recover a proper pair with original generators."""
    if not isinstance(list_of_pairs, list):
        raise ValueError('[Error]: Given instance is not a list.')
    for item in list_of_pairs:
        if not isinstance(item,properpair.ProperPair):
            raise ValueError('[Error]: The given list contains an element which is not a proper pair.')
    hash_list = [item._hashstr() for item in list_of_pairs]
    unique_hash_list = list(set(hash_list)) # Find unique hash lists.
    index_pairs = [hash_list.index(hashstring) for hashstring in unique_hash_list]
    return [list_of_pairs[ind] for ind in index_pairs]