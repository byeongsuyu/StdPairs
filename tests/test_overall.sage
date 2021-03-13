from stdpairs import *
import numpy as np
def print_err(str_var):
    return str_var + "does not work properly."
print("This is test suits for stdpairs.")
print("running AffineMonoid and its methods . . .", end =" ")
A=matrix(ZZ, [[1,1,1],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]).T
Q=AffineMonoid(A)
if (matrix(ZZ,Q.mingens()) != A[:,1:5]): raise SyntaxError(print_err(".gens()"))
if len(Q.index_to_face()) != 11:
    raise SyntaxError(print_err(".index_to_face()"))
if str(type(Q.poly())) != "<class 'sage.geometry.polyhedron.parent.Polyhedra_ZZ_ppl_with_category.element_class'>":
    raise SyntaxError(print_err(".poly()"))
if set(Q.face_lattice()) != set(Q.index_to_face().keys()):
    raise SyntaxError(print_err(".face_lattice() or .index_to_face()")) 
if set(Q.poly().face_lattice()) != set(Q.index_to_face().values()):
    raise SyntaxError(print_err(".poly() or .index_to_face()")) 
if len(Q.integral_support_vectors()) != 10:
    raise SyntaxError(print_err(".integral_support_vectors()"))
temp_list = list([np.array(face.ambient_Hrepresentation()).T[1:].T for face in Q.poly().face_lattice()])
del temp_list[0]
ttemp_list = list(Q.integral_support_vectors().values())
if not np.all([np.array_equal(temp_list[i],ttemp_list[i]) for i in range(len(temp_list))]):
    raise SyntaxError(print_err(".integral_support_vectors()"))
del temp_list,ttemp_list
if not np.array_equal(A,Q.gens()):
    raise SyntaxError(print_err(".gens()"))
if (Q.is_empty() != False) or (Q.is_pointed() != True):
    raise SyntaxError(print_err(".is_empty"))
if not np.all([face_tup== Q.index_of_face(Q.face(face_tup)) for face_tup in list(Q.index_to_face().keys())[1:]]):
    raise SyntaxError(print_err(".face() or .index_of_face()"))
if not ((not Q.is_element(matrix(ZZ,[[2],[3],[1]]))) and (Q.is_element(matrix(ZZ,[[4],[7],[7]])))):
    raise SyntaxError(print_err(".is_element()"))
if not ('Q\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\n' == Q.save_txt()):
    raise SyntaxError(print_err(".save_txt()"))
print(" pass")
print("running txt_to_affinemonoid() . . .", end ="")
if (Q!=txt_to_affinemonoid('Q\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\n')): raise SyntaxError(print_err("txt_to_affinemonoid()"))
print(" pass")
print("running txt_to_monomialideal() . . .", end ="")
J = txt_to_monomialideal('I\nQ\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\nis_empty_ideal\n0\ngens\n1,1,2|1,2,2|3,3,3\n__is_std_cover_calculated\n1\n{"(1, 3)": ["0|0|0&(1, 3)"], "(1, 2)": ["0|0|0&(1, 2)"], "(0, 3, 4)": ["0|0|0&(0, 3, 4)"], "(0, 2, 4)": ["0|0|0&(0, 2, 4)"], "(2,)": ["0|1|1&(2,)"]}\n__is_overlap_calculated\n1\n{"(1, 3)": [["0|0|0&(1, 3)"]], "(1, 2)": [["0|0|0&(1, 2)"]], "(0, 3, 4)": [["0|0|0&(0, 3, 4)"]], "(0, 2, 4)": [["0|0|0&(0, 2, 4)"]], "(2,)": [["0|1|1&(2,)"]]}\n__is_max_overlap_calculated\n1\n{"(1, 3)": [["0|0|0&(1, 3)"]], "(1, 2)": [["0|0|0&(1, 2)"]], "(0, 3, 4)": [["0|0|0&(0, 3, 4)"]], "(0, 2, 4)": [["0|0|0&(0, 2, 4)"]], "(2,)": [["0|1|1&(2,)"]]}\n__is_ass_prime_calculated\n1\n{"(1, 3)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n1,1|0,1|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(1, 2)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|1,1|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(0, 3, 4)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|0,0|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(0, 2, 4)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0|0,1|1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "(2,)": "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0,1|0,1,1|1,1,1\\n__is_std_cover_calculated\\n0\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n"}\n__is_irr_decom_prime_calculated\n1\n["I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n1,1|0,1|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(1, 3)\\": [\\"0|0|0&(1, 3)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|1,1|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(1, 2)\\": [\\"0|0|0&(1, 2)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,1|0,0|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(0, 3, 4)\\": [\\"0|0|0&(0, 3, 4)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0|0,1|1,1\\n__is_std_cover_calculated\\n1\\n{\\"(0, 2, 4)\\": [\\"0|0|0&(0, 2, 4)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n", "I\\nQ\\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\\nis_empty_ideal\\n0\\ngens\\n0,0,0,1,2|0,1,2,2,2|2,2,2,2,2\\n__is_std_cover_calculated\\n1\\n{\\"(2,)\\": [\\"0|0|1&(2,)\\", \\"1|1|1&(2,)\\", \\"0|0|0&(2,)\\", \\"0|1|1&(2,)\\"]}\\n__is_overlap_calculated\\n0\\n__is_max_overlap_calculated\\n0\\n__is_ass_prime_calculated\\n0\\n__is_irr_decom_prime_calculated\\n0\\n"]')
print(" pass")
print("running MonomialIdeal and its methods . . .", end ="")
I = MonomialIdeal(matrix(ZZ,[[1,1,3],[1,2,3],[2,2,3]]).T, Q)
if (I==J) != True: raise SyntaxError(print_err(".__eq__()"))
if (I.ambient_monoid() == J.ambient_monoid()) != True: raise SyntaxError(print_err(".ambient_monoid()"))
if (np.array_equal(I.gens(), J.gens())) != True: raise SyntaxError(print_err(".gens()"))
if (I.associated_primes() == J.associated_primes()) != True: raise SyntaxError(print_err(".associated_primes()"))
if (I.irreducible_decomposition() == J.irreducible_decomposition()) != True: raise SyntaxError(print_err(".irreducible_decomposition()"))
if (I.is_element(J.gens()[:,[0]]) == J.is_element(I.gens()[:,[0]])) != True: raise SyntaxError(print_err(".is_element()"))
if (I.is_empty() == J.is_empty()) != True: raise SyntaxError(print_err(".is_empty()"))
if (I.is_irreducible() == J.is_irreducible()) != True: raise SyntaxError(print_err(".is_irreducible()"))
if (I.is_primary() == J.is_primary()) != True: raise SyntaxError(print_err(".is_primary()"))
if (I.is_prime() == J.is_prime()) != True: raise SyntaxError(print_err(".is_prime()"))
if (I.radical() == J.radical()) != True: raise SyntaxError(print_err(".radical()"))
if (I.is_radical() == J.is_radical()) != True: raise SyntaxError(print_err(".is_radical()"))
if (I.is_principal() == J.is_principal()) != True: raise SyntaxError(print_err(".is_principal()"))
if (I.is_standard_monomial(J.gens()[:,[0]]) == J.is_standard_monomial(I.gens()[:,[0]])) != True: raise SyntaxError(print_err(".is_standard_monomial()"))
if (I.maximal_overlap_classes() == J.maximal_overlap_classes()) != True: raise SyntaxError(print_err(".is_standard_monomial()"))
if ([I.multiplicity(item) for item in J.associated_primes()] == [J.multiplicity(item) for item in I.associated_primes()]) != True: raise SyntaxError(print_err(".multiplicity()"))
if (I.overlap_classes() == J.overlap_classes()) != True: raise SyntaxError(print_err(".overlap_classes()"))
if (I.standard_cover() == J.standard_cover()) != True: raise SyntaxError(print_err(".standard_cover()"))
H=J.associated_primes()[(1,3)].intersect(J.associated_primes()[(1,2)])
T=txt_to_monomialideal('I\nQ\n1,0,1,0,1|1,0,0,1,1|1,1,1,1,1\nis_empty_ideal\n0\ngens\n1|1|1\n__is_std_cover_calculated\n1\n{"(1, 3)": ["0|0|0&(1, 3)"], "(1, 2)": ["0|0|0&(1, 2)"]}\n__is_overlap_calculated\n0\n__is_max_overlap_calculated\n0\n__is_ass_prime_calculated\n0\n__is_irr_decom_prime_calculated\n0\n')
if (T == H) != True: raise SyntaxError(print_err(".intersect()"))
if ( I== I+J ) != True: raise SyntaxError(print_err(".__add__()"))
if ((I*J) == MonomialIdeal(matrix(ZZ, [[2, 2, 3, 2, 3, 4],[2, 3, 3, 4, 4, 4],[6, 6, 6, 6, 6, 6]]), Q)) != True: raise SyntaxError(print_err(".__mul__()"))
print("pass")
print("running prime_ideal() . . .", end ="")
if prime_ideal((), Q) != MonomialIdeal(Q.gens(),Q): raise SyntaxError(print_err("prime_ideal()"))
print(" pass")
print("running ProperPair class and its methods . . .", end ="")
P = list(I.standard_cover().values())[0][0]
if (np.array_equal(P.monomial(),0*P.ambient_ideal().gens()[:,[0]])) != True: raise SyntaxError(print_err(".monomial()"))
if (P.face() == (1,3)) != True: raise SyntaxError(print_err(".__eq__()"))
if (P.ambient_ideal() == I) != True: raise SyntaxError(print_err(".ambient_ideal()"))
if (P.is_maximal() == True) != True: raise SyntaxError(print_err(".is_maximal()"))
if (P.is_element(matrix(ZZ,[[1],[1],[1]])) == False) != True: raise SyntaxError(print_err(".is_element()"))
print(" pass")
print("running div_pairs() . . .", end ="")
P = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,4), I)            
T = ProperPair(matrix(ZZ,[[1], [1], [1]]), (3,), I)
if (div_pairs(P,T) != False): raise SyntaxError(print_err("div_pairs()"))
TT = ProperPair(matrix(ZZ,[[1], [1], [1]]), (0,4), I)
if (np.array(div_pairs(P,TT)).size != 7): raise SyntaxError(print_err("div_pairs()"))
print(" pass")
print("running macaulay2() . . .", end="")
try: 
    aaa=macaulay2('3/5 + 7/11') 
    del aaa
except: 
    raise SyntaxError("Interface for Macaulay2 does not work properly.")
print(" pass")
print("running from_macaulay2() . . .", end ="")
R = macaulay2.eval('ZZ[x,y,z]')
test = macaulay2.eval("loadPackage Normaliz")
test = macaulay2.eval('S=createMonomialSubalgebra {x^5*y, y*z^2, z^3}')
if from_macaulay2("S").save_txt() != 'Q\n5,0,0|1,1,0|0,2,3\n': raise SyntaxError(print_err("from_macaulay2()"))
print(" pass")
print("running to_macaulay2() . . .", end ="")
D=to_macaulay2(J)
if (D['AffineSemigroupRing'] != macaulay2.eval('R')) or (D['MonomialIdeal'] != macaulay2.eval('I')) or (D['StandardCover'] != macaulay2.eval('SC')): raise SyntaxError(print_err("to_macaulay2()"))
print(" pass")
print("StdPair package test was done! Everything works well.")











