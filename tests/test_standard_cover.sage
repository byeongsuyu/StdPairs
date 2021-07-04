print("test MonomialIdeal.standard_cover() . . . .", end="")
from stdpairs import AffineMonoid, MonomialIdeal
import numpy as np
import ast
numvar = 3
numgen = 3
A = np.identity(numvar,dtype="int64")
B = np.random.randint(3, size=(numvar, numgen))
Q = AffineMonoid(A)
I = MonomialIdeal(B,Q)
# Calculate it by Macaulay2
# Set a polynomial ring first.
temp=macaulay2.eval('R=ZZ[vars (0..'+str(numvar-1)+')]')
gens_ideal='I=monomialIdeal(';
for row in I.gens().T:
    eq='' 
    for idx in range(numvar): 
        eq=eq+'R_'+str(idx)+'^'+str(row[idx])+'*' 
    eq = eq[:-1] 
    gens_ideal=gens_ideal + eq+',' 
gens_ideal=gens_ideal[:-1] + ')'
temp=macaulay2.eval(gens_ideal)
temp=macaulay2.eval('L=standardPairs I')
result=list(macaulay2('T=apply(L, i -> {(exponents i_0)_0, apply(i_1, j -> index j)})')) 
for item in result: temp=item[1].sort();
faces = list(set([tuple(item[1]) for item in result])) 
faces = [ast.literal_eval(str(face)) for face in faces]
cover_from_mac2={}
for face in faces: 
    cover_from_mac2[face]=[] 
    for item in result: 
        if face == tuple(item[1]): 
            cover_from_mac2[face].append(ast.literal_eval(str(tuple(item[0]))))
# Calculate standard pairs by stdpairs package
cover_from_sage ={}
for face,list_of_pairs in I.standard_cover().items(): 
    cover_from_sage[face] = [tuple(pair.monomial().T.tolist()[0]) for pair in list_of_pairs]
# Check equality:
faces=list(cover_from_mac2.keys())+list(cover_from_sage.keys())
faces=list(set([str(item) for item in faces]))
faces=[ast.literal_eval(item) for item in faces]
for face in faces:
    if set([str(item) for item in cover_from_sage[face]]) != set([str(item) for item in cover_from_mac2[face]]):
        raise SyntaxError("Standard pairs from Macaulay2 are not equal to those from stdpairs package.")
print(" pass")