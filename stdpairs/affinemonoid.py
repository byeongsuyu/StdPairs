r"""
AffineMonoid

This module provides the base class :class:`AffineMonoid` in StdPairs.

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
from pathlib import Path
from sage.rings.integer_ring import Z as ZZ
from sage.matrix.constructor import Matrix as matrix
from sage.misc.persist import save as sobj_save
import sage.geometry.polyhedron.constructor as const

from . import _zsolve_interface
from . import _string_interface

class AffineMonoid():
    r"""
    A class representing an *affine semigroup*. An *affine semigroup* is a semigroup 
    generated by a finite subset of a free abelian group :math:`\mathbb{Z}^{n}`
    for some integer :math:`n`. The finite subset used for construction is represented as
    a matrix over :math:`\mathbb{Z}` and called the *generating matrix*.

    INPUT:

    - ``genset`` -- A ``NumPy.ndarray`` object with 2-dimensional shape and integer elements or a ``sage.matrix.matrix_integer_dense`` object.
    - ``is_normaliz`` -- A ``bool`` object. Default is ``False``. If it is ``True``, generate corresponding polyhedron using ``PyNormaliz`` package.

    OUTPUT:

    - An ``AffineMonoid`` object whose generating set is ``genset``.

    EXAMPLE:

    ::

        sage: #Using ``NumPy``
        sage: from stdpairs import AffineMonoid
        sage: import numpy
        sage: A = numpy.array([[1,2],[0,2]])
        sage: A
        array([[1, 2],
               [0, 2]])
        sage: Q = AffineMonoid(A)
        sage: Q
        An affine semigroup whose generating set is 
        [[1 2]
         [0 2]]

    ::

        sage: #Using ``sage.matrix.matrix_integer_dense``
        sage: from stdpairs import AffineMonoid
        sage: A = matrix(ZZ,[[1,2],[0,2]])
        sage: Q = AffineMonoid(A)
        sage: Q
        An affine semigroup whose generating set is 
        [[1 2]
         [0 2]]

    ::

        sage: #Turn on Normaliz and execute its method.
        sage: # Be careful that Normaliz does not distinguish 
        sage: # non-normal semigroup from its saturated one.
        sage: from stdpairs import AffineMonoid
        sage: A = matrix(ZZ,[[1,2],[0,2]])
        sage: Q = AffineMonoid(A, True)
        sage: Q
        An affine semigroup whose generating set is 
        [[1 2]
         [0 2]]
        sage: Q.poly().hilbert_series([1,0])
        1/(t^2 - 2*t + 1)
        sage: A = matrix(ZZ,[[1,1],[0,1]])                                              
        sage: Q = AffineMonoid(A, True)                                                 
        sage: Q.poly().hilbert_series([1,0])                                            
        1/(t^2 - 2*t + 1)

    TESTS::
            
        sage: from stdpairs import AffineMonoid
        sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
        sage: Q = AffineMonoid(A)
    """

    def __init__(self, genset, is_normaliz=False):
        """ Initialize ``self`` using the generating set ``genset``."""
        # Argument type handling.
        if not isinstance(is_normaliz, bool):
            raise ValueError("[Error] 2nd argument must be a boolean value for determining whether use PyNormliz or not.")
        type_mat = type(matrix(ZZ,0))
        if isinstance(genset, type_mat):
            genset = np.array(genset).astype('int64')
        del type_mat
        if not isinstance(genset,np.ndarray):
            raise ValueError("[Error]: 1st argument is not numpy.ndarray type")
        if len(genset.shape)!=2:
            raise ValueError("[Error] Not a valid affine semigroup ring; parameter should be 2D array.")
        if (genset.dtype != 'int32') and (genset.dtype != 'int64'):
            raise ValueError('[Error] genset should have dtype int32 or int64.')
        
        #Initialize the generators of affine monoid
        #Notes that we sort genset by lexicographical order of columns; this is required for equality.
        self.__gens=genset
            
        #Initialize corresponding real polyhedron.
        if ~(self.__gens.any(axis=0).any()): # When the matrix is just zero matrix,
            self.__is_zeromonoid_bool =True
            #Create a normal integral polyhedron with vertex 0.
            if is_normaliz:
                self.__poly = const.Polyhedron(vertices=list(map(tuple, self.__gens)),backend='normaliz', base_ring = ZZ)
            else:
                self.__poly = const.Polyhedron(vertices=list(map(tuple, self.__gens)), base_ring = ZZ)
        else:
            self.__is_zeromonoid_bool =False
            #Creat a normal integral polyhedron whose rays are columns of self.__gens
            if is_normaliz:
                self.__poly=const.Polyhedron(rays=(self.__gens.transpose()).tolist(),backend='normaliz', base_ring = ZZ)
            else:
                self.__poly=const.Polyhedron(rays=(self.__gens.transpose()).tolist(),base_ring = ZZ)
      
        # Determine whether affine monoid is pointed or not.
        if self.__poly.faces(0):
            self.__is_pointed_bool=True
        else:
            self.__is_pointed_bool=False
        #Initialize the dictionary of integral support vectors
        self.__integral_support_vectors={}

        # Generate Face Dictionary - Start
        # 1) Get a face lattice of real polyhedron.
        # 2) Set the top and bottom (especially 0-dim and -1-dim faces) as a tuple.
        #    -1 dim face will be replaced with (-1,)
        #     0 dim face will be replaced with () (empty tuple)
        #     top dim face will be replaced with self.__genset
        lattice_dict={}
        lattice=(self.__poly.face_lattice()).list()
        # For -1-dimensional face, set index (-1,)
        if self.__poly.faces(-1):
            lattice_dict.update({self.__poly.faces(-1)[0]:(-1,)})
            del lattice[0]
        # For any other faces,

        for face in lattice:
            # For given face, find an integral support vectors of face
            hvectors= np.array([ np.array(eq.A()) for eq in face.ambient_Hrepresentation()]).astype('int64')
            if hvectors.size == 0:
                hvectors == np.array([]).astype('int64')
            # Now find indices of face as a submatrix of self.__gens
            if hvectors.any():
                ind=tuple(np.where((np.matmul(hvectors, self.__gens)).sum(axis=0) == 0)[0])
            else:
                ind=tuple([i for i in range(self.__gens.shape[1])])
            #Sort indices
            ind = tuple(sorted(ind))
            #Update dictionaries
            lattice_dict.update({face:ind})
            self.__integral_support_vectors.update({ind:hvectors})
        # Relabel facelattice for an affine monoid.
        self.__face_lattice=(self.__poly.face_lattice()).relabel(lattice_dict)
        # Save index_to_face as {Ind:face} form.
        self.__index_to_face={v: k for k, v in lattice_dict.items()}
        del lattice_dict
        del lattice
        
        # Calculate minimal generators
        if self.__is_zeromonoid_bool:
            self.__mingens = self.__gens
        else:
            if not self.__is_pointed_bool:
                self.__mingens = _zsolve_interface._mingens_for_nonpointed(self.__gens)
            else:
                self.__mingens = _zsolve_interface._mingens(self.__gens, self.__gens)
        # This attribute is for hashstring. This is used for saving the generator.
        self.__hashstring = _string_interface._np2d_to_string(self.__mingens)
    def index_to_face(self):
        r"""
        returns a ``dictionary`` whose keys are ``tuple`` denoting column indices of the generating set and whose values are corresponding face of real polyhedron as ``sage.geometry.polyhedron.face.PolyhedronFace`` object.
    
        OUTPUT:

        - ``index_to_face`` -- ``dictionary`` object.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]]) 
            sage: Q = AffineMonoid(A)
            sage: Q.index_to_face()
            {(): A 0-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex,
             (-1,): A -1-dimensional face of a Polyhedron in ZZ^2,
             (0,): A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 1 ray,
             (0,
              1,
              2): A 2-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 rays,
             (1,): A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 1 ray}

        """
        return self.__index_to_face
    def integral_support_vectors(self):
        r"""
        returns a ``dictionary`` whose keys are ``tuple`` denoting face as column indices of the generating set, and whose values are an integer inner normal vector of the face. See [MY2020]_ for details about integral support function.
    
        OUTPUT:

        - ``integral support vectors`` -- A ``dictionary`` object.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]]) 
            sage: Q = AffineMonoid(A)
            sage: Q.integral_support_vectors()
            {(): array([[ 1, -1],
                    [ 0,  1]]),
             (0,): array([[0, 1]]),
             (0, 1, 2): array([], dtype=int64),
             (1,): array([[ 1, -1]])}

        """
        return self.__integral_support_vectors
    def poly(self):
        r"""
        returns a polyhedron in :math:`\mathbb{Q}^{n}` generated by the given affine monoid.
        Notes that even the affine monoid non-normal, this returns only the normal rational cone, i.e., consisting of all lattice points of a rational polyhedron.
        See `Polyhedra <https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/constructor.html>`_ for details to use the object returned.
    
        OUTPUT:

        - ``poly`` -- A ``sage.geometry.polyhedron.parent.Polyhedra_ZZ_ppl_with_category`` object.

        EXAMPLE:

        ::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]]) 
            sage: Q = AffineMonoid(A)
            sage: Q.poly()
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 rays
            sage: B = matrix(ZZ, [[1,1],[0,1]]) 
            sage: R = AffineMonoid(B)
            sage: Q.poly() == R.poly()                                                      
            True
            sage: # They are differ as a polyhedron over ZZ, but the same as an object in SageMath.

        Especially, if you turn on Normaliz option by construction of your affine semigroup, then you can use all methods in normaliz. For example,

        ::

            sage: A = matrix(ZZ,[[1,1],[0,1]])                                              
            sage: Q = AffineMonoid(A, True)                                                 
            sage: Q.poly().hilbert_series([1,0])                                            
            1/(t^2 - 2*t + 1)

        Be careful that Normaliz does not distinguish a non-normal affine semigroup with its saturated one.

        """
        return self.__poly
    def face_lattice(self):
        r"""
        returns a finite lattice of of tuples representing faces ordered by inclusion. Notes that :math:`(-1,)` is the least element for every affine monoid, denoting (-1)-dimesional face. This is not the same as :math:`()`, the 0-dimensional face, which exists only in a pointed affine monoid.
        See `Finite lattice and semilattices <https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/posets/lattices.html>`_ for details to use the object.
    
        OUTPUT:

        - ``face_lattice`` -- A ``sage.combinat.posets.lattices.FiniteLatticePoset_with_category`` object.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])
            sage: Q = AffineMonoid(A) 
            sage: Q.face_lattice()
            Finite lattice containing 5 elements
            sage: list(Q.face_lattice())                                                    
            [(-1,), (), (0,), (1,), (0, 1, 2)]

        """
        return self.__face_lattice
    def mingens(self):
        r"""
        returns the minimal generating set of the given affine monoid as an ``NumPy.ndarray`` object with 2-dimensional shape if ``self`` is pointed affine monoid. If ``self`` is not pointed,
        then return the generating set with warning message.
        
    
        OUTPUT:

        - ``genset`` -- A ``NumPy.ndarray`` object with 2-dimensional shape.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])
            sage: Q = AffineMonoid(A) 
            sage: Q.mingens()
            array([[1, 2],
                   [0, 2]])

        """
        return self.__mingens

    def gens(self):
        r"""
        returns the generating set of the given affine monoid as an ``NumPy.ndarray`` object with 2-dimensional shape. This may not be the minimal generating set.
        In case of the pointed affine monoid, use ``self.mingens()`` to get the minimal generating set.
    
        OUTPUT:

        - ``genset`` -- A ``NumPy.ndarray`` object with 2-dimensional shape.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])
            sage: Q = AffineMonoid(A) 
            sage: Q.gens()
            array([[1, 2, 7],
                   [0, 2, 4]])

        """
        return self.__gens
    def is_empty(self):
        r"""
        returns ``True`` iff it is an empty set as an affine monoid .

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])
            sage: Q = AffineMonoid(A) 
            sage: Q.is_empty()
            False
            sage: R = AffineMonoid(matrix(ZZ,0))
            sage: R.is_empty()
            True
        """
        return self.__is_zeromonoid_bool
    def is_pointed(self):
        r"""
        returns ``True`` iff the given affine monoid is pointed, i.e., no units except 0.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])
            sage: Q = AffineMonoid(A) 
            sage: Q.is_pointed()
            True
            sage: R = AffineMonoid(matrix(ZZ,[[1,2,7,-1],[0,2,4,0]]))
            sage: R.is_pointed()
            False
        """
        return self.__is_pointed_bool
    def save_txt(self):
        r"""
        returns strings containing information about the given affine monoid. One can recover it using a global function ``txt_to_affinemonoid()``.

        EXAMPLE::
                                    
            sage: from stdpairs import AffineMonoid                                         
            sage: from stdpairs import txt_to_affinemonoid                                  
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])                                         
            sage: Q = AffineMonoid(A)                                                       
            sage: Q.save_txt()                                                              
            'Q\n1,2,7|0,2,4\n'
            sage: R=txt_to_affinemonoid(Q.save_txt())                                       
            sage: R == Q                                                                    
            True
        """
        return "Q\n"+_string_interface._np2d_to_string(self.__gens)+"\n"
    def save(self, path_of_file):
        r"""
        saves the given ``AffineMonoid`` object as sobj file. See `Object Persistence <https://doc.sagemath.org/html/en/reference/misc/sage/misc/persist.html>`_ for detail.
        
        INPUT:

        - ``path_of_file`` -- a ``string`` object denoting the path which the object ``self`` will be saved as a binary file.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid                                         
            sage: from pathlib import Path
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])                                         
            sage: Q = AffineMonoid(A)                                                       
            sage: Q == loads(dumps(Q))
            True
        """
        # If path_of_file starts with ~, change it as a Home directory.
        if not isinstance(path_of_file,str):
            raise ValueError("[Error]: The given instance is not a valid path")
        if path_of_file[0] == '~':
            path_of_file = str(Path.home())+path_of_file[1:]
        sobj_save(self, path_of_file)

    def face(self, index):
        r"""
        returns a set of generators of the given face represented by a tuple consisting of column indices of the generating matrix.
        
        INPUT:

        - ``index`` -- A ``tuple`` object containing column indices of the generating matrix which consists a face.

        OUTPUT:

        - A ``NumPy.ndarray`` object representing a face.

        EXAMPLE:

        ::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])                                         
            sage: Q = AffineMonoid(A)
            sage: Q.face((0,))                                                              
            array([[1],
                   [0]])
            sage: Q.face((1,))                                                              
            array([[2],
                   [2]])
            sage: Q.face((0,1,2))                                                           
            array([[1, 2, 7],
                   [0, 2, 4]])
            
        """
        if not isinstance(index,tuple):
            raise ValueError("[Error]:Instance is not a tuple.")
        if index in self.__face_lattice:
            if index != (-1,):
                return self.__gens[:,index]
            else:
                return np.array([], dtype = 'int64')
        else:
            raise ValueError("[Error]:Instance is not a face of the affine monoid")
            
    def index_of_face(self, face_genset):
        r"""
        returns a tuple denoting set of generators of the given face represented by the input ``face_genset``, a matrix consisting of columns representing a face.
        
        INPUT:

        - ``face_genset`` -- A ``NumPy`` 2D-array with integer elements or ``sage.matrix.matrix_integer_dense`` object.

        OUTPUT:

        - A ``tuple`` object representing the generating matrix of a face of the affine semigroup.

        EXAMPLE:

        ::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])                                         
            sage: Q = AffineMonoid(A)
            sage: F = matrix(ZZ,[[1],[0]])
            sage: Q.index_of_face(F) 
            (0,)
        """
        type_mat = type(matrix(ZZ,0))
        if isinstance(face_genset, type_mat):
            face_genset = np.array(face_genset).astype('int64')
        del type_mat
        if len(face_genset.shape) !=2:
            raise ValueError("[Error]: an argument is not 2D numpy array.")
        if face_genset.shape[1]==0:
            return ()
        # Collect all column indices of self.__gens which are also a column of the given face.
        # Return the collected indices as a tuple.
        j=0
        ind=[]
        for i in range(self.__gens.shape[1]):
            if np.array_equal(face_genset[:,(j,)],self.__gens[:,(i,)]):
                j=j+1
                ind.append(i)
            if j == face_genset.shape[1]:
                return tuple(ind)
        raise ValueError("[Error]:An argument's column dimension exceed the generators of ambient monoid.")

    def is_element(self, col_vector): 
        r"""
        Suppose :math:`A` is the generating matrix and :math:`b` is a vector given by the input ``col_vector``. Then
        this function returns a ``sage.matrix.matrix_integer_dense`` object whose row :math:`u` satisfying :math:`A \cdot u^{t}=b`.
        In other words, if the input is an element of the given affine semigroup, it returns non-empty matrix whose rows show how to generated
        the input using the generators of the affine monoid. Otherwise, if the input is not the element of the given affine monoid, it returns 
        the empty matrix.

        INPUT:

        - ``col_vector`` -- A ``NumPy.ndarray`` object with 2-dimensional shape with integer coefficients or ``sage.matrix.matrix_integer_dense`` object with one column. In other words, this argument should be a vector.

        OUTPUT:

        - A ``sage.matrix.matrix_integer_dense`` object whose row :math:`u` is the solution of the equation :math:`A \cdot u^{t}=b`.

        EXAMPLE::

            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])
            sage: Q = AffineMonoid(A) 
            sage: Q.is_element(matrix(ZZ,[[7],[4]]))
            [0 0 1]
            [3 2 0]
            sage: Q.is_element(matrix(ZZ,[[-1],[0]]))
            []
        """
        type_mat = type(matrix(ZZ,0))
        if isinstance(col_vector, type_mat):
            col_vector = np.array(col_vector).astype('int64')
        del type_mat
        if len(col_vector.shape)!=2:
            raise ValueError("[Error]:The given instance is not a 2D-array")
        else:
            if col_vector.shape[1]!=1:
                raise ValueError("[Error]:The given instance is not a column vector")
        return _zsolve_interface._zsolve_exact_solver(self.__gens, col_vector)

    def _intersection_of_pairs(self, monomial_a, face_a, monomial_b, face_b):
        """ 
        Given two proper pairs (``monomial_a``,``face_a``) and (``monomial_b``,``face_b``) find nonnegative vectors :math:`u` and :math:`v` such that

        ::

            monomial_a+ face_a * u = monomial_b+ face_b * v. 
            

        INPUT:

        - ``monomial_a``,``monomial_b``  -- ``NumPy.ndarray`` objects with 2-dimensional shape, integer cofficients and only 1 column or a ``sage.matrix.matrix_integer_dense`` object with 1 column. 
        - ``face_a``, ``face_b`` -- ``tuple`` objects in ``self.face_lattice()``.

        OUTPUT:

        - A ``sage.matrix.matrix_integer_dense`` object whose row is of form :math:`[u;v]`, which gives the solution of the equation above. 
        
        """
        #To have optimization, we intended to use is_element.
        type_mat = type(matrix(ZZ,0))
        if isinstance(monomial_a, type_mat):
            monomial_a = np.array(monomial_a).astype('int64')
        if isinstance(monomial_b, type_mat):
            monomial_b = np.array(monomial_b).astype('int64')
        del type_mat
        if not(isinstance(monomial_a, np.ndarray) and isinstance(monomial_a, np.ndarray)):
            raise ValueError("[Error]:At least one of monomials is not a type of np.ndarray.")
        if (len(monomial_a.shape)!=2) or (len(monomial_b.shape)!=2):
            raise ValueError("[Error]:At least one of monomials is not a 2D array")
        if (monomial_a.shape[1] != 1) or (monomial_b.shape[1] != 1):
            raise ValueError("[Error]:At least one of monomials has more than 1 column")
        if (not face_a in self.__face_lattice.list()) or (not face_b in self.__face_lattice.list()):
            raise ValueError("[Error]:At least one of faces is not in the face of ambient monoid")
        #Generate an affine monoid with pair_a's face with - pair_b's face.
        temp_gens=np.concatenate((self.face(face_a),-self.face(face_b)), axis=1)
        if (temp_gens.size ==0):
            if (monomial_b - monomial_a).any() == True:
                return matrix(ZZ,0)
            else:
                return matrix(ZZ,monomial_a*0)
        solution_vectors = AffineMonoid(temp_gens).is_element(monomial_b - monomial_a)
        #Return the minimal solutions
        return solution_vectors
    def __str__(self):
        return "An affine semigroup whose generating set is \n"+np.array2string(self.__gens)
    def __repr__(self):
        return str(self)
    def __eq__(self, other):
        """ 
        returns ``True`` iff two ``AffineMonoid`` objects are equal, i.e., having the same *minimal* generating matrix.

        
        EXAMPLE::
            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ, [[1,2,7],[0,2,4]])
            sage: Q = AffineMonoid(A) 
            sage: Q.is_element(matrix(ZZ,[[7],[4]]))
            [0 0 1]
            [3 2 0]
            sage: Q.is_element(matrix(ZZ,[[-1],[0]]))
            []
        
        TESTS::
            
            sage: from stdpairs import AffineMonoid
            sage: A = matrix(ZZ,[[0,1,1,0],[0,0,1,1],[1,1,1,1]])
            sage: Q = AffineMonoid(A) 
            sage: G = SymmetricGroup(A.ncols())
            sage: for permutation in G: 
            ....:     if (Q== AffineMonoid(A*permutation)) == False: 
            ....:         raise SyntaxError("[Error]: __eq__ does not work properly.")
            
        """
        if not isinstance(other, type(self)): return NotImplemented
        return _zsolve_interface._check_swap_column_equivalence(self.mingens(), other.mingens())
