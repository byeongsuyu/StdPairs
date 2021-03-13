r"""

Call package:

::

    sage: from stdpairs import *
    #'AffineMonoid', 'MonomialIdeal', 'ProperPair' classes and
    #'prime_ideal','div_pairs','pair_difference','from_macaulay2'
    #'to_macaulay2','txt_to_affinemonoid','txt_to_monomialideal'
    # global functions are imported.

AUTHORS:

- Byeongsu Yu (2021-02-25): initial version.

TESTS:

::

	sage: from stdpairs import *
	sage: import numpy as np
	sage: def print_err(str_var):
	....:     return str_var + "does not work properly."
	sage: print("This is test suits for stdpairs.")
	sage: print("running AffineMonoid and its methods . . .", end =" ")
	sage: A=matrix(ZZ, [[1,1,1],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]).T
	sage: Q=AffineMonoid(A)
	sage: if (matrix(ZZ,Q.mingens()) != A[:,1:5]): raise SyntaxError(print_err(".gens()"))
	sage: if len(Q.index_to_face()) != 11:
	....:     raise SyntaxError(print_err(".index_to_face()"))
	sage: if str(type(Q.poly())) != "<class 'sage.geometry.polyhedron.parent.Polyhedra_ZZ_ppl_with_category.element_class'>":
	....:     raise SyntaxError(print_err(".poly()"))
	sage: if set(Q.face_lattice()) != set(Q.index_to_face().keys()):
	....:     raise SyntaxError(print_err(".face_lattice() or .index_to_face()")) 
	sage: if set(Q.poly().face_lattice()) != set(Q.index_to_face().values()):
	....:     raise SyntaxError(print_err(".poly() or .index_to_face()")) 
	sage: if len(Q.integral_support_vectors()) != 10:
	....:     raise SyntaxError(print_err(".integral_support_vectors()"))
	sage: temp_list = list([np.array(face.ambient_Hrepresentation()).T[1:].T for face in Q.poly().face_lattice()])
	sage: del temp_list[0]
	sage: ttemp_list = list(Q.integral_support_vectors().values())
	sage: if not np.all([np.array_equal(temp_list[i],ttemp_list[i]) for i in range(len(temp_list))]):
	....:     raise SyntaxError(print_err(".integral_support_vectors()"))
	sage: del temp_list,ttemp_list
	sage: if not np.array_equal(A,Q.gens()):
	....:     raise SyntaxError(print_err(".gens()"))
	sage: if (Q.is_empty() != False) or (Q.is_pointed() != True):
	....:     raise SyntaxError(print_err(".is_empty"))
	sage: if not np.all([face_tup== Q.index_of_face(Q.face(face_tup)) for face_tup in list(Q.index_to_face().keys())[1:]]):
	....:     raise SyntaxError(print_err(".face() or .index_of_face()"))
	sage: if not ((not Q.is_element(matrix(ZZ,[[2],[3],[1]]))) and (Q.is_element(matrix(ZZ,[[4],[7],[7]])))):
	....:     raise SyntaxError(print_err(".is_element()"))
	sage: if not ('Q\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\n' == Q.save_txt()):
	....:     raise SyntaxError(print_err(".save_txt()"))
	sage: print(" pass")
	sage: print("running txt_to_affinemonoid() . . .", end ="")
	sage: if (Q!=txt_to_affinemonoid('Q\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\n')): raise SyntaxError(print_err("txt_to_affinemonoid()"))
	sage: print(" pass")
	sage: print("running txt_to_monomialideal() . . .", end ="")
	sage: J = txt_to_monomialideal('I\nQ\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\nis_empty_ideal\n0\ngens\n1,1,2|1,2,2|3,3,3\n__is_std_cover_calculated\n1\n{"(1, 3)": ["0|0|0&(1, 3)"], "(1, 2)": ["0|0|0&(1, 2)"], "(0, 3, 4)": ["0|0|0&(0, 3, 4)"], "(0, 2, 4)": ["0|0|0&(0, 2, 4)"], "(2,)": ["0|1|1&(2,)"]}\n__is_overlap_calculated\n1\n{"(1, 3)": [["0|0|0&(1, 3)"]], "(1, 2)": [["0|0|0&(1, 2)"]], "(0, 3, 4)": [["0|0|0&(0, 3, 4)"]], "(0, 2, 4)": [["0|0|0&(0, 2, 4)"]], "(2,)": [["0|1|1&(2,)"]]}\n__is_max_overlap_calculated\n1\n{"(1, 3)": [["0|0|0&(1, 3)"]], "(1, 2)": [["0|0|0&(1, 2)"]], "(0, 3, 4)": [["0|0|0&(0, 3, 4)"]], "(0, 2, 4)": [["0|0|0&(0, 2, 4)"]], "(2,)": [["0|1|1&(2,)"]]}\n__is_ass_prime_calculated\n1\n{"(1, 3)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n1,1|0,1|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(1, 2)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|1,1|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(0, 3, 4)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|0,0|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(0, 2, 4)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0|0,1|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(2,)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0,1|0,1,1|1,1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n"}\n__is_irr_decom_prime_calculated\n1\n["I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n1,1|0,1|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(1, 3)\\": [\\"0|0|0&(1, 3)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|1,1|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(1, 2)\\": [\\"0|0|0&(1, 2)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|0,0|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(0, 3, 4)\\": [\\"0|0|0&(0, 3, 4)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0|0,1|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(0, 2, 4)\\": [\\"0|0|0&(0, 2, 4)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0,0,1,2|0,1,2,2,2|2,2,2,2,2\\n__is_std_cover_calculated\\n1\\n{\\"(2,)\\": [\\"0|0|1&(2,)\\", \\"1|1|1&(2,)\\", \\"0|0|0&(2,)\\", \\"0|1|1&(2,)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n"]')
	sage: print(" pass")
	sage: print("running MonomialIdeal and its methods . . .", end ="")
	sage: I = MonomialIdeal(matrix(ZZ,[[1,1,3],[1,2,3],[2,2,3]]).T, Q)
	sage: if (I==J) != True: raise SyntaxError(print_err(".__eq__()"))
	sage: if (I.ambient_monoid() == J.ambient_monoid()) != True: raise SyntaxError(print_err(".ambient_monoid()"))
	sage: if (np.array_equal(I.gens(), J.gens())) != True: raise SyntaxError(print_err(".gens()"))
	sage: if (I.associated_primes() == J.associated_primes()) != True: raise SyntaxError(print_err(".associated_primes()"))
	sage: if (I.irreducible_decomposition() == J.irreducible_decomposition()) != True: raise SyntaxError(print_err(".irreducible_decomposition()"))
	sage: if (I.is_element(J.gens()[:,[0]]) == J.is_element(I.gens()[:,[0]])) != True: raise SyntaxError(print_err(".is_element()"))
	sage: if (I.is_empty() == J.is_empty()) != True: raise SyntaxError(print_err(".is_empty()"))
	sage: if (I.is_irreducible() == J.is_irreducible()) != True: raise SyntaxError(print_err(".is_irreducible()"))
	sage: if (I.is_primary() == J.is_primary()) != True: raise SyntaxError(print_err(".is_primary()"))
	sage: if (I.is_prime() == J.is_prime()) != True: raise SyntaxError(print_err(".is_prime()"))
	sage: if (I.radical() == J.radical()) != True: raise SyntaxError(print_err(".radical()"))
	sage: if (I.is_radical() == J.is_radical()) != True: raise SyntaxError(print_err(".is_radical()"))
	sage: if (I.is_principal() == J.is_principal()) != True: raise SyntaxError(print_err(".is_principal()"))
	sage: if (I.is_standard_monomial(J.gens()[:,[0]]) == J.is_standard_monomial(I.gens()[:,[0]])) != True: raise SyntaxError(print_err(".is_standard_monomial()"))
	sage: if (I.maximal_overlap_classes() == J.maximal_overlap_classes()) != True: raise SyntaxError(print_err(".is_standard_monomial()"))
	sage: if ([I.multiplicity(item) for item in J.associated_primes()] == [J.multiplicity(item) for item in I.associated_primes()]) != True: raise SyntaxError(print_err(".multiplicity()"))
	sage: if (I.overlap_classes() == J.overlap_classes()) != True: raise SyntaxError(print_err(".overlap_classes()"))
	sage: if (I.standard_cover() == J.standard_cover()) != True: raise SyntaxError(print_err(".standard_cover()"))
	sage: H=J.associated_primes()[(1,3)].intersect(J.associated_primes()[(1,2)])
	sage: T=txt_to_monomialideal('I\nQ\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\nis_empty_ideal\n0\ngens\n1|1|1\n__is_std_cover_calculated\n1\n{"(1, 3)": ["0|0|0&(1, 3)"], "(1, 2)": ["0|0|0&(1, 2)"]}\n__is_overlap_calculated\n0\n__is_max_overlap_calculated\n0\n__is_ass_prime_calculated\n0\n__is_irr_decom_prime_calculated\n0\n')
	sage: if (T == H) != True: raise SyntaxError(print_err(".intersect()"))
	sage: if ( I== I+J ) != True: raise SyntaxError(print_err(".__add__()"))
	sage: if ((I*J) == MonomialIdeal(matrix(ZZ, [[2, 2, 3, 2, 3, 4],[2, 3, 3, 4, 4, 4],[6, 6, 6, 6, 6, 6]]), Q)) != True: raise SyntaxError(print_err(".__mul__()"))
	sage: print("pass")
	sage: print("running prime_ideal() . . .", end ="")
	sage: if prime_ideal((), Q) != MonomialIdeal(Q.gens(),Q): raise SyntaxError(print_err("prime_ideal()"))
	sage: print(" pass")
	sage: print("running ProperPair class and its methods . . .", end ="")
	sage: P = list(I.standard_cover().values())[0][0]
	sage: if (np.array_equal(P.monomial(),0*P.ambient_ideal().gens()[:,[0]])) != True: raise SyntaxError(print_err(".monomial()"))
	sage: if (P.face() == (1,3)) != True: raise SyntaxError(print_err(".__eq__()"))
	sage: if (P.ambient_ideal() == I) != True: raise SyntaxError(print_err(".ambient_ideal()"))
	sage: if (P.is_maximal() == True) != True: raise SyntaxError(print_err(".is_maximal()"))
	sage: if (P.is_element(matrix(ZZ,[[1],[1],[1]])) == False) != True: raise SyntaxError(print_err(".is_element()"))
	sage: print(" pass")
	sage: print("running div_pairs() . . .", end ="")
	sage: P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,4), I)            
	sage: T = ProperPair(matrix(ZZ,[[1], [1], [1]]), (3,), I)
	sage: if (div_pairs(P,T) != False): raise SyntaxError(print_err("div_pairs()"))
	sage: TT = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,4), I)
	sage: if (np.array(div_pairs(P,TT)).size != 7): raise SyntaxError(print_err("div_pairs()"))
	sage: print(" pass")
	sage: print("running macaulay2() . . .", end="")
	sage: try: 
	....:     aaa=macaulay2('3/5 + 7/11') 
	....:     del aaa
	....: except: 
	....:     raise SyntaxError("Interface for Macaulay2 does not work properly.")
	sage: print(" pass")
	sage: print("running from_macaulay2() . . .", end ="")
	sage: R = macaulay2.eval('ZZ[x,y,z]')
	sage: test = macaulay2.eval("loadPackage Normaliz")
	sage: test = macaulay2.eval('S=createMonomialSubalgebra {x^5*y, y*z^2, z^3}')
	sage: if from_macaulay2("S").save_txt() != 'Q\n5,0,0|1,1,0|0,2,3\n': raise SyntaxError(print_err("from_macaulay2()"))
	sage: print(" pass")
	sage: print("running to_macaulay2() . . .", end ="")
	sage: D=to_macaulay2(J)
	sage: if (D['AffineSemigroupRing'] != macaulay2.eval('R')) or (D['MonomialIdeal'] != macaulay2.eval('I')) or (D['StandardCover'] != macaulay2.eval('SC')): raise SyntaxError(print_err("to_macaulay2()"))
	sage: print(" pass")
	sage: print("StdPair package test was done! Everything works well.") 




"""

#*****************************************************************************
#       Copyright (C) 2021 Byeongsu Yu <byeongsu.yu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; version 3 of
#  the License. 
#                  http://www.gnu.org/licenses/
#*****************************************************************************

__all__ = ['AffineMonoid', 'MonomialIdeal', 'ProperPair', 'prime_ideal','div_pairs','pair_difference','from_macaulay2','to_macaulay2','txt_to_affinemonoid','txt_to_monomialideal']


from .affinemonoid import AffineMonoid
from .monomialideal import MonomialIdeal, prime_ideal
from .properpair import ProperPair, div_pairs
from ._stdpairs import pair_difference
from .macaulay2_interface import from_macaulay2,to_macaulay2
from .txt_to_object import txt_to_affinemonoid, txt_to_monomialideal




    