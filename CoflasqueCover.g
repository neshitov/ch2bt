###############################################################################
# Wrapper around algorithm for coflabby resolution by A. Hoshi and A. Yamasaki
###############################################################################

Read("./RatProbAlgTori/FlabbyResolution.gap");

CoflasqueCover := function(gens)
  local g, subgroups, proj;
  g := Group(List(gens, x -> TransposedMat(x)),
             IdentityMat(DimensionsMat(gens[1])[1]));
  subgroups :=List(ConjugacyClassesSubgroups2(g),Representative);
  proj:=FindCoflabbyResolutionBase(g, subgroups);
  proj:=TransposedMat(proj);

  return proj;
end;


LeftInverse := function(A)
  local g;
  g := NullspaceIntMat(A);
  g := Concatenation(TransposedMat(A), g)^-1;
  g := TransposedMat(g);
  return g{[ 1 .. DimensionsMat(A)[2]]};
end;


CoflasqueResolution := function(pr, gens)
  local cr, inverse_injection, t_p, gen;
  cr := rec();
  cr.surjection := pr;
  cr.injection := NullspaceIntMat(TransposedMat(pr));
  cr.injection := TransposedMat(cr.injection);
  inverse_injection := LeftInverse(cr.injection);
  t_p := TransposedMat(pr);
  cr.actionP := List(gens,
                     x -> TransposedMat(PermutationMat(Permutation(TransposedMat(x), t_p), DimensionsMat(pr)[2]))
                    );
  cr.actionC :=List(cr.actionP, x -> inverse_injection * x * cr.injection);

  return cr;
end;
