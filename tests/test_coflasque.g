##############################################
# Tests for coflasque resolutions
#############################################

Read("./RatProbAlgTori/FlabbyResolution.gap");

CheckResolution := function(gens, res)
  local diff, i;
  if Length(gens) <> Length(res.actionP) then
    Error("Wrong number of P gens");
  fi;
  if Length(gens) <> Length(res.actionC) then
    Error("Wrong number of C gens");
  fi;

  for i in [1..Length(gens)] do
    diff:= gens[i] * res.surjection - res.surjection * res.actionP[i];
    if diff <> NullMat(DimensionsMat(diff)[1], DimensionsMat(diff)[2]) then
      Error("Projection equivariance Error");
    fi;
  od;

  for i in [1..Length(gens)] do
    diff:= res.actionP[i] * res.injection - res.injection * res.actionC[i];
    if diff <> NullMat(DimensionsMat(diff)[1], DimensionsMat(diff)[2]) then
      Error("Injection equivariance Error");
    fi;
  od;

  if not IsCoflabby(TransposedMatrixGroup(Group(res.actionC))) then
    Error("Kernel is not coflasque");
  fi;
  Print("pass");
end;

#gens := Carat(5, 697, 4);
#cr := CoflasqueResolution(CoflasqueCover(gens), gens);
#CheckResolution(gens, cr);
