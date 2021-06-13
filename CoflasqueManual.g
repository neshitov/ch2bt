Read("CaratReader.g");

InvSubspace := function(hh, n)
  local num_gens, eqs, inv;
  num_gens := Length(hh);
  eqs :=  MatrixByBlockMatrix(
                  BlockMatrix(List([1 .. Length(hh)],
                                   i -> [i, 1, hh[i] - IdentityMat(n, Integers)]),
                              num_gens,
                              1)
  );
  inv := NullspaceMat(TransposedMat(eqs));
  inv := TransposedMat(inv);
  return inv;
end;


gens := Carat(6, 6137, 18);
G := Group(gens);
subgens := SmallGeneratingSet(Representative(ConjugacyClassesSubgroups(G)[2]));

#G := Group(gens);
#Print(Order(G));
#Print("nr subgroups\n");
#Print(Length(ConjugacyClassesSubgroups(G)));
