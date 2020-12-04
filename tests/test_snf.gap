# Function to test Smith Normal Form

CheckAnswer := function(tr, A)
  local nrows, ncols, power, i, j;
  nrows := DimensionsMat(A)[1];
  ncols := DimensionsMat(A)[2];
  power := PrimePowersInt(Characteristic(Ring(A[1][1])))[2];

  for i in [1 .. nrows ] do
    for j in [1 .. ncols ] do
      if j <> i and tr.SNF[i][j] <> Zero(Integers mod 2^power) then
        Error("SNF not diagonal");
      fi;
    od;
  od;

  for i in [ 1 .. tr.rank ] do
    if tr.SNF[i][i] = Zero(Integers mod 2^power) then
      Error( " rank incorrect " );
    fi;
  od;

  for i in [ tr.rank + 1 .. Minimum(nrows, ncols) ] do
    if tr.SNF[i][i] <> Zero(Integers mod 2^power) then
      Error( " rank incorrect " );
    fi;
  od;

  if (tr.row_t * A * tr.col_t) <> tr.SNF then
    Error( "Transformation is not equal to SNF " );
  fi;

  if (tr.row_t * tr.row_t_inverse) <> IdentityMat(nrows) * One(Integers mod 2^power) then
    Error( "Row transform inverse fail " );
  fi;
  Print("pass");

end;

Read("../LocalSmithNormalForm.gap");
x := RandomMat(20, 15) * One(Integers mod 512);;
tr := SNFTransform(x: num_threads:=4);;
CheckAnswer(tr, x);
