"""Common module"""

import numpy as np
from sage.all import matrix
from sage.all import ZZ
from sage.interfaces.four_ti_2 import four_ti_2

import sage.geometry.polyhedron.constructor as const


#### zsolve cases; used in all 
def _zsolve_exact_solver(matrix_b, col_vector, all_solution=False):
    """Find solution u of matrix_b*u=vector. If return_boolean is True, return all results from zsolve."""
    if not isinstance(matrix_b,np.ndarray):
        raise ValueError("[Error]: First instance,", matrix_b ," is not numpy array.")
    elif len(matrix_b.shape) !=2:
        raise ValueError("[Error]: First instance,", matrix_b ," is not 2D numpy array.")
    if not isinstance(col_vector,np.ndarray):
        raise ValueError("[Error]: Second instance,", col_vector ," is not numpy array.")
    elif len(col_vector.shape) !=2:
        raise ValueError("[Error]: First instance,", matrix_b ," is not 2D numpy array.")
    elif col_vector.shape[1] != 1:
        raise ValueError("[Error]: First instance,", matrix_b ," is not a column vector, i.e., its shape[1] != 1")
    # Make a matrix submatrix_b=[A;A]
    submatrix_b=np.concatenate((matrix_b,matrix_b))
    # Transpose a 2D column col_vector into 1D col_vector.
    col_vector=col_vector.flatten()
    # sign for zsolve
    # If submatrix_b is n x m matrix, then sgn is m x 1 matrix whose entries are 1.
    # 1 implies that solution should be nonnegative.
    sgn = [1 for column in submatrix_b.T]
    # relation for zsolve.
    # What we have is [<=; >=], which is needed to solve an integer equation.
    rel=['<'] * matrix_b.shape[0] + ['>']*matrix_b.shape[0]
    # Solve it using zsolve.
    result=four_ti_2.zsolve(submatrix_b.tolist(), rel, col_vector.tolist()+col_vector.tolist(), sgn)
    if all_solution == True:
        return result
    return_sol = result[0]
    inf_sol = result[2]
    del result
    # Notes that result[0]: Inhomogeneous solution Bx = [col_vector;col_vector]
    # result[1]: Homogeneous solution Bx=0
    # result[2]: Affine monoid whose solution set lives in
    if inf_sol: #If solution is infinite,
        raise ValueError("[Error]: Solution space is not finite, which is impossible in this case.")
    return return_sol

def _mingens_for_nonpointed(genset):
    r"""
    This method returns the minimal generators of given generators in cases
    when the genset is not pointed. The result is not unique; it is dependent to orders of columns of genset.

    """
    living_cols = [0]
    col_num = genset.shape[1]
    for idx in range(1,col_num): 
        P = const.Polyhedron(rays=genset[:,living_cols].T.tolist(),base_ring=ZZ) 
        if P.contains(genset[:,idx].T.tolist()): 
            temp=_zsolve_exact_solver(genset[:,living_cols],genset[:,idx][np.newaxis].T) 
            if np.array(temp).size == 0: 
                living_cols.append(idx) 
        else: 
            living_cols.append(idx) 
    return genset[:,living_cols]


def _mingens(genset_of_base_monoid, genset):
    """ This method returns a minimal generators of the given generators of an ideal.
        This is private member; only call by generator.
    """
    #Check that if it is an empty (zero) ideal.
    if (genset.size ==0):
        return np.zeros((genset_of_base_monoid.shape[0],), dtype=int)[np.newaxis].T
    # Sort columns of genset by an ascending order of its norm.
    ind_col=np.argsort(np.linalg.norm(genset,axis=0))
    sorted_gens = (genset.T[ind_col]).T
    # Check that if the ideal is nonproper.
    column_index =0
    parity = True
    while(sorted_gens.size >0 and parity == True):
        if np.linalg.norm(sorted_gens.T[0])==0:
            sorted_gens = sorted_gens[:,1:]
        else:
            parity = False
    if sorted_gens.any()== False: # If sorted gens are 0,
        #Return the zero column only as a generator
        return(np.zeros((genset_of_base_monoid.shape[0],), dtype=int)[np.newaxis].T)
    #Because of for loop updates, we use the transposed one.
    sorted_gens = sorted_gens.T
    idx_for_gens=0
    # To update our loop, we use while instead of for loop.
    while idx_for_gens < sorted_gens.shape[0]-1:
        # Pick a column which is the leftmost part.
        column = sorted_gens[idx_for_gens]
        #Fix all columns which are right of the chosen column in the genset.
        temp_sorted_gens=np.delete(sorted_gens,list(range(idx_for_gens+1)),0)
        #Initialize record for column indices.
        death_note = [];
        for idx, sub_column in enumerate(temp_sorted_gens):
            #Given all columns which are right of the chosen column,
            #Test weather Ax+column = sub_column has a nonnegative solution
            if (_zsolve_exact_solver(genset_of_base_monoid,(sub_column-column)[np.newaxis].T)).nrows()>0:
                #If this has a solution,
                # Record the index of that subcolumn. Since we don't need it for generating the given ideal.
                death_note.append(idx_for_gens+idx+1)
        #Delete all columns which are redundant (i.e., reachable by Ax+fixed column for some x.)
        sorted_gens = np.delete(sorted_gens,death_note,0)
        #Increase index for next choice.
        idx_for_gens=idx_for_gens+1
    return(sorted_gens.T)

def _check_swap_column_equivalence(matrix_a,matrix_b):
    """Two matrices must have the same shape.
    """
    type_mat = type(matrix(ZZ,0))
    if isinstance(matrix_a, type_mat):
        matrix_a = np.array(matrix_a).astype('int64')
    if isinstance(matrix_b, type_mat):
        matrix_b = np.array(matrix_b).astype('int64')
    if not isinstance(matrix_a, np.ndarray):
        raise ValueError("[Error]:1st argument is not numpy.ndarray")
    if not isinstance(matrix_b, np.ndarray):
        raise ValueError("[Error]:2nd argument is not numpy.ndarray")
    if (matrix_a.shape == (0,)) and (matrix_b.shape == (0,)):
        return True
    if (matrix_a.shape != matrix_b.shape):
        raise ValueError("[Error]:Two arguments have difference shapes")
    if len(matrix_a.shape)!=2:
        raise ValueError("[Error]:Two arguments are not numpy 2D-array")
        
    col_list_self = list(map(tuple,matrix_a.T))
    col_list_other = list(map(tuple,matrix_b.T))
    if len(col_list_self) == len(col_list_other):
        if len(set(col_list_self+col_list_other)) == len(set(col_list_self)):
            return True
        else:
            return False
    else:
        return False