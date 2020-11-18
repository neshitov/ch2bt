# This code tests the SparsReduction algorithm

Read("SparseReduction.gap");

for i in [1..10] do
  M := RandomMat(10, 5, Integers);
  Ms := SparseMatrix(M);
  tr := ReduceMatrix(Ms);

  if tr * SparseMatrix(M) <> Ms then
    Error(M);
  else
    Print("Pass\n");
  fi;
od;
