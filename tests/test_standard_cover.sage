sage: print("test MonomialIdeal.standard_cover() . . . .", end="")
sage: from stdpairs import AffineMonoid, MonomialIdeal
sage: import numpy as np
sage: import ast
sage: numvar = 3
sage: numgen = 3
sage: A = np.identity(numvar,dtype="int64")
sage: B = np.random.randint(3, size=(numvar, numgen))
sage: Q = AffineMonoid(A)
sage: I = MonomialIdeal(B,Q)
sage: # Calculate it by Macaulay2
sage: # Set a polynomial ring first.
sage: temp=macaulay2.eval('R=ZZ[vars (0..'+str(numvar-1)+')]')
sage: gens_ideal='I=monomialIdeal(';
sage: for row in I.gens().T:
....:     eq='' 
....:     for idx in range(numvar): 
....:         eq=eq+'R_'+str(idx)+'^'+str(row[idx])+'*' 
....:     eq = eq[:-1] 
....:     gens_ideal=gens_ideal + eq+',' 
sage: gens_ideal=gens_ideal[:-1] + ')'
sage: temp=macaulay2.eval(gens_ideal)
sage: temp=macaulay2.eval('L=standardPairs I')
sage: result=list(macaulay2('T=apply(L, i -> {(exponents i_0)_0, apply(i_1, j -> index j)})')) 
sage: for item in result: temp=item[1].sort();
sage: faces = list(set([tuple(item[1]) for item in result])) 
sage: faces = [ast.literal_eval(str(face)) for face in faces]
sage: cover_from_mac2={}
sage: for face in faces: 
....:     cover_from_mac2[face]=[] 
....:     for item in result: 
....:         if face == tuple(item[1]): 
....:             cover_from_mac2[face].append(ast.literal_eval(str(tuple(item[0]))))
sage: # Calculate standard pairs by stdpairs package
sage: cover_from_sage ={}
sage: for face,list_of_pairs in I.standard_cover().items(): 
....:     cover_from_sage[face] = [tuple(pair.monomial().T.tolist()[0]) for pair in list_of_pairs]
sage: # Check equality:
sage: faces=list(cover_from_mac2.keys())+list(cover_from_sage.keys())
sage: faces=list(set([str(item) for item in faces]))
sage: faces=[ast.literal_eval(item) for item in faces]
sage: for face in faces:
....:     if set([str(item) for item in cover_from_sage[face]]) != set([str(item) for item in cover_from_mac2[face]]):
....:         raise SyntaxError("Standard pairs from Macaulay2 are not equal to those from stdpairs package.")
sage: print(" pass")