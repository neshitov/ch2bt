# This code tests the SparsReduction algorithm

Read("SparseReduction.gap");

IsREF :=function(sm)
  local pivot, i;
  pivot := -1;
  for i in [1..Nrows(sm)] do
    if sm!.indices[i] <> [] then
      if sm!.indices[i][1] <= pivot then
        return false;
      else
        pivot:=sm!.indices[i][1];
      fi;
    else
      pivot := Ncols(sm) + 1;
    fi;
  od;
  return true;
end;


for i in [1..10] do
  M := RandomMat(10, 5, Integers);
  Ms := SparseMatrix(M);
  tr := ReduceMatrix(Ms);
  tr_dense := ConvertSparseMatrixToMatrix(tr.transform);

  if AbsInt(Determinant(tr_dense)) <> 1 then
    Error("det");
  fi;

  if not IsREF(tr.transform * SparseMatrix(M)) then
    Error(M);
  else
    Print("Pass\n");
  fi;
od;
