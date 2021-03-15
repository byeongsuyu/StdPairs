r"""
ProperPair

This module provides the base class :class:`ProperPair` and a method :method:`div_pair()` in StdPairs.

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
from sage.all import ZZ
from sage.all import matrix

from . import monomialideal
from . import _zsolve_interface
from . import _string_interface

class ProperPair:
    r"""
    A class representing a proper pair with respect to given monomial ideal :math:`I`, which is a pair :math:`(a,F)` where :math:`a` is a monomial
    in the ambient affine monoid, and :math:`F` denotes a face such that :math:`a+\mathbb{N}F \cap I = \emptyset`.
    Or, abusively, this class also generated without checking the condition :math:`a+\mathbb{N}F \cap I = \emptyset` using ``properness`` argument.

    INPUT:

    - ``monomial`` -- A ``NumPy.ndarray`` object with 2-dimensional shape and integer elements or ``sage.matrix.matrix_integer_dense`` type variable. In any cases, this argument should have only one column.
    - ``face`` -- A ``tuple`` object representing a face in ``ambient_ideal.ambient_monoid()``. 
    - ``ambient_ideal`` -- A ``MonomialIdeal`` object.
    - ``properness=False`` -- A ``bool`` object. If this argument is ``True``, then this class does not give an error even if :math:`a+\mathbb{N}F \cap I` is nonempty.

    OUTPUT:

    - ``ProperPair`` object :math:`(a,F)` whose :math:`a` is ``monomial`` argument and :math:`F` is ``face`` argument.

    EXAMPLE:

    ::

        sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
        sage: import numpy
        sage: A = matrix(ZZ,[[1,2],[0,2]])
        sage: Q = AffineMonoid(A)
        sage: M = matrix(ZZ,[[4,6],[4,6]])
        sage: I = MonomialIdeal(M,Q)
        sage: #Using ``NumPy``
        sage: # Below ``P`` is a proper pair.
        sage: P = ProperPair(numpy.array([[2],[2]]), (0,), I)
        sage: # But below ``Not_P`` raise error that this is not a proper pair.
        sage: # Not_P =ProperPair(matrix(ZZ,[[2],[2]]), (1,), I)
        sage: # However, using ``properness`` as ``True`` the same construction does not give you error;
        sage: # Basically the module does not check whether the given pair is proper or not.
        sage: Not_P =ProperPair(matrix(ZZ,[[2],[2]]), (1,), I,True) 
        
    ::

        sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
        sage: A = matrix(ZZ,[[1,2],[0,2]])
        sage: Q = AffineMonoid(A)
        sage: M = matrix(ZZ,[[4,6],[4,6]])
        sage: I = MonomialIdeal(M,Q)
        sage: #Using ``NumPy``
        sage: # Below ``P`` is a proper pair.
        sage: P = ProperPair(matrix(ZZ,[[2],[2]]), (0,), I)

    TESTS::
            
        sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
        sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
        sage: Q = AffineMonoid(A)
        sage: I = MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q)
        sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,), I)
        sage: if (P.is_maximal() != False):
        ....:     raise SyntaxError("Method is_maximal() is problematic; please report it to the developer.")
        sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (3,), I)
        sage: if (P.is_maximal() != False):
        ....:     raise SyntaxError("Method is_maximal() is problematic; please report it to the developer.")
        sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,3), I)
        sage: if (P.is_maximal() != True):
        ....:     raise SyntaxError("Method is_maximal() is problematic; please report it to the developer.")

    """
    def __init__(self, monomial, face, ambient_ideal, properness=False):
        if not isinstance(ambient_ideal,monomialideal.MonomialIdeal):
            raise ValueError("[Error]: The 3rd parameter is not a type of 'MonomialIdeal'")
        type_mat = type(matrix(ZZ,0))
        if isinstance(monomial, type_mat):
            if monomial == matrix(ZZ,0):
                monomial = 0*(ambient_ideal.gens())[:,0][np.newaxis].T
            elif monomial.nrows() != (ambient_ideal.ambient_monoid()).gens().shape[0]:
                raise ValueError("[Error]: The monomial is not in the same dimension with the generators of the ambient monoid.")
            else:
                monomial = np.array(monomial).astype('int64')
        del type_mat
        if not properness:
            if (ambient_ideal.ambient_monoid()).is_element(monomial).nrows()==0:
                raise ValueError("[Error]: The monomial is not in the ambient monoid.")
            if not ambient_ideal.is_standard_monomial(monomial):
                raise ValueError("[Error]: The monomial is in the given monomial ideal.")
        if not(face in ((ambient_ideal.ambient_monoid()).face_lattice()).list()):
            raise ValueError("[Error]: Given tuple is not in the faces of the ambient monoid.")
        # Check whether ideal is empty, or we ignore properness.
        if ((ambient_ideal.is_empty()) or (properness)):
            self.__monomial = monomial
            self.__face = face
            self.__ambient_ideal = ambient_ideal
            self.__checking_well_defined = False
        else:
            self.__checking_well_defined = True
            #Test that whether it is proper pair or not.
            for column in (ambient_ideal.gens()).T:
                # Check that monomial + face *x = column + Monoid element
                if (ambient_ideal.ambient_monoid())._intersection_of_pairs(monomial, face, column[np.newaxis].T, ((ambient_ideal.ambient_monoid()).face_lattice()).top()).nrows():
                    # The reason for setting Warning is that it will be used as a status below.
                    raise Warning("[Error]: This is not a proper pair of the given ideal.")
            self.__monomial = monomial
            self.__face = face
            self.__ambient_ideal = ambient_ideal
        #Now set hashstring
        self.__hashstring = _string_interface._np2d_to_string(self.__monomial)+"&"+ str(self.__face)

        #Check maximality over the division
        self.__bool_is_maximal_with_divisibility_calculated= False;
        self.__bool_is_maximal_with_divisibility= False;
    def monomial(self):
        r"""
        returns the monomial of ``self``.
    
        OUTPUT:

        - ``monomial`` -- A ``NumPy.ndarray`` object with 2-dimensional shape, integer elements, and only one column.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
            sage: Q=AffineMonoid(matrix(ZZ,[[1,2],[0,2]]))
            sage: I=MonomialIdeal(matrix(ZZ,[[4,6],[4,6]]),Q)
            sage: P=ProperPair(matrix(ZZ,[[2],[2]]),(0,),I)
            sage: P.monomial()
            array([[2],
                   [2]])

        """
        return self.__monomial
    def face(self):
        r"""
        returns the face of ``self`` as a ``tuple`` object.
    
        OUTPUT:

        - ``gens`` -- A ``tuple`` object.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
            sage: Q=AffineMonoid(matrix(ZZ,[[1,2],[0,2]]))
            sage: I=MonomialIdeal(matrix(ZZ,[[4,6],[4,6]]),Q)
            sage: P=ProperPair(matrix(ZZ,[[2],[2]]),(0,),I)
            sage: P.face()
            (0,)
            sage: Q.face((0,))
            array([[1],
                   [0]])
        """
        return self.__face
    def ambient_ideal(self):
        r"""
        returns an ideal which decides whether ``self`` is proper (i.e., a translation of submonoid does not intersect with the ideal) or not. This was given when ``self`` was created.
    
        OUTPUT:

        - ``ambient_ideal`` -- A ``MonomialIdeal`` object.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
            sage: Q=AffineMonoid(matrix(ZZ,[[1,2],[0,2]]))
            sage: I=MonomialIdeal(matrix(ZZ,[[4,6],[4,6]]),Q)
            sage: P=ProperPair(matrix(ZZ,[[2],[2]]),(0,),I)
            sage: P.ambient_ideal()
            An ideal whose generating set is 
            [[4]
             [4]]
        """
        return self.__ambient_ideal
    def is_maximal(self):
        r"""
        returns ``True`` if and only if a pair is maximal with respect to divisibility. A proper pair :math:`(a,F)` is maximal with respect to divisibility
        if for any proper pair :math:`(b,G)` there is no such element :math:`c` in the ambient affine monoid such that :math:`a+c+\mathbb{N}F \subseteq b+\mathbb{N}G`.
        Notes that this function is not properly working if you generated the given pair with properness=`False` argument.
        
        OUTPUT:

        - A ``bool`` object.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
            sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
            sage: Q = AffineMonoid(A)
            sage: I = MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q)
            sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,), I)
            sage: P.is_maximal()
            False
            sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,3), I)
            sage: P.is_maximal()
            True

        """
        if self.__bool_is_maximal_with_divisibility_calculated:
            return self.__bool_is_maximal_with_divisibility
        if self.__checking_well_defined == False:
            warnings.warn("[Warning]: Since the pair itself may not proper, it is impossible to check maximality with respect to divisibility.")
        else:
            top_face = (((self.__ambient_ideal.ambient_monoid()).face_lattice()).maximal_elements())[0]
            parity = 0; # Check whether we have proper pair divisible by given std pair.
            for ind in top_face:
                if ind not in self.__face:
                    # Get that column in gens.
                    gen = ((self.__ambient_ideal.ambient_monoid()).gens())[:,ind][np.newaxis].T
                    try:
                        ProperPair(gen+self.__monomial, self.__face, self.__ambient_ideal)
                        self.__bool_is_maximal_with_divisibility = False
                        self.__bool_is_maximal_with_divisibility_calculated = True
                        return self.__bool_is_maximal_with_divisibility
                    except:
                        pass
            self.__bool_is_maximal_with_divisibility = True
            self.__bool_is_maximal_with_divisibility_calculated = True
            return self.__bool_is_maximal_with_divisibility
    def is_element(self, a_monomial):
        r"""
        Given ``a_monomial`` monomial, say :math:`b`, in the ambient monoid, find a matrix whose row :math:`u` is a solution of
        :math:`F*u =b-a` where :math:`(a,F)` is the give proper pair ``self``. In other words, if ``a_monomial`` in :math:`a+\mathbb{N}F`, return a matrix
        whose row :math:`u` satisfies :math:`a+F*u=b`. Otherwise, return an empty matrix.

        INPUT:

        - ``a_monomial`` -- A ``NumPy.ndarray`` object with 2-dimensional shape, integer elements, and only one vector, or ``sage.matrix.matrix_integer_dense`` object with one column.

        OUTPUT:

        - A ``sage.matrix.matrix_integer_dense`` object.

        EXAMPLE::
            
            sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
            sage: Q=AffineMonoid(matrix(ZZ,[[1,2],[0,2]]))
            sage: I=MonomialIdeal(matrix(ZZ,[[3,7],[2,0]]),Q)
            sage: P= ProperPair(matrix(ZZ,[[2],[2]]),(1,),I)
            sage: P.is_element(matrix(ZZ,[[8],[8]]))                                                                                                        
            [3]
            sage: # This means that [[8],[8]]- [[2],[2]] = 3*[[2],[2]]
            sage: P.is_element(matrix(ZZ,[[7],[7]]))                                                                                                        
            []
            sage: # This does not have solution since [[7],[7]] cannot be obtained by integer.
            
        """
        type_mat = type(matrix(ZZ,0))
        if isinstance(a_monomial, type_mat):
            a_monomial = np.array(a_monomial).astype('int64')
        del type_mat
        if not isinstance(a_monomial, np.ndarray):
            raise ValueError("[Error]: The given instance is not a numpy or sage.matrix object")
        if len(a_monomial.shape)!=2:
            raise ValueError("[Error]:The given instance is not a numpy 2D array")
        else:
            if a_monomial.shape[1]!=1:
                raise ValueError("[Error]:The given instance is not a column vector")
        # If face is zero, then just check whether it is the same as monomial or not.
        if self.__face == ():
            if np.any(self.__monomial-a_monomial) == False:
                # Then they are the same
                return matrix(ZZ,[0])
            else:
                return matrix(ZZ,0)
        # Check it for the submonoid gen by face.
        return _zsolve_interface._zsolve_exact_solver((((self.ambient_ideal()).ambient_monoid()).gens())[:,self.__face], a_monomial-self.__monomial)
    def _hashstr(self):
        r"""This is used for string saving."""
        return self.__hashstring
    def __str__(self):
        return "("+str(self.__monomial.tolist())+"^T,"+str((((self.ambient_ideal()).ambient_monoid()).face(self.__face)).tolist())+")"
    def __repr__(self):
        return str(self)
    def __eq__(self, other):
        r""" 
        returns ``True`` iff two ``ProperPair`` objects are equal, i.e., having the same generating monomial and the same face.

        
        EXAMPLE::

            sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair
            sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
            sage: Q = AffineMonoid(A)
            sage: I = MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q)
            sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,), I)            
            sage: T = ProperPair(matrix(ZZ,[[1], [1], [1]]), (3,), I)
            sage: P == T
            False
        """
        if not isinstance(other, type(self)): return NotImplemented
        if self.ambient_ideal() == other.ambient_ideal():
            if (self.__face == other.face()) and _zsolve_interface._check_swap_column_equivalence(self.__monomial,other.monomial()):
                return True
            else:
                return False
        else:
            warnings.warn("[Warning]: their ambient ideals are not the same.")
            return False
def div_pairs(pair_a, pair_b):
    r""" 
        Given two pairs ``pair_a`` and ``pair_b`` representing :math:`(a,F)` and :math:`(b,G)` respectively,
        this method returns a set of monomials :math:`c` such that :math:`a+c+\mathbb{N}F \subseteq b+\mathbb{N}G` as a matrix.
        
        INPUT:

        - ``pair_a``,``pair_b`` -- ``ProperPair`` objects.

        OUTPUT:

        - A ``sage.matrix.matrix_integer_dense`` object whose columns are minimal solutions of :math:`a+c+\mathbb{N}F \subseteq b+\mathbb{N}G`.
        
        EXAMPLE::

            sage: from stdpairs import AffineMonoid, MonomialIdeal, ProperPair, div_pairs
            sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
            sage: Q = AffineMonoid(A)
            sage: I = MonomialIdeal(matrix(ZZ,[[2,2,2],[0,1,2],[2,2,2]]),Q)
            sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,), I)            
            sage: T = ProperPair(matrix(ZZ,[[1], [1], [1]]), (3,), I)
            sage: div_pairs(P,T)
            []
            sage: TT = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,3), I)
            sage: div_pairs(P,TT)
            [0]
            [0]
            [0]
    
    """

    """If self=(a,F) divides other=(b,G), then return c where a+c+NF subseteq b+NG.

        If two pairs are not comparable, gives a warning message and return an empty matix.
           
        Warning: In case of two faces are the same, then we determine it by coordinatewise
        lexicographical order. However, this need proof.
    """
    if not isinstance(pair_a,ProperPair):
        raise ValueError("[Error]: The 1st parameter is not a type of 'ProperPair'")
    if not isinstance(pair_b,ProperPair):
        raise ValueError("[Error]: The 2nd parameter is not a type of 'ProperPair'")
    if pair_a.ambient_ideal() != pair_b.ambient_ideal():
        raise ValueError("[Error]: Two has different ambient ideals")
    tfone=(((pair_a.ambient_ideal()).ambient_monoid()).face_lattice()).is_lequal(pair_a.face(), pair_b.face())
    if (not tfone):
        #Not divisible
        return matrix(ZZ,0)
    else:
        #Get an index of the whole monoid as a face
        temp_face= ((pair_a.ambient_ideal()).ambient_monoid()).index_of_face(((pair_a.ambient_ideal()).ambient_monoid()).gens())
        numgen = len(temp_face)
        #Generate a pair (self.monomial, ambient_monoid, zero ideal)
        return matrix(ZZ,((pair_a.ambient_ideal()).ambient_monoid().gens()))*(((pair_a.ambient_ideal()).ambient_monoid())._intersection_of_pairs(pair_a.monomial(),temp_face, pair_b.monomial(),pair_b.face())[:,0:numgen]).T


