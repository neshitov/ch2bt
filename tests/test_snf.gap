# Function to test Smith Normal Form

CheckAnswer := function(tr, A, power)
  local nrows, ncols, i, j, one;
  one := One(Integers mod 2^power);
  nrows := DimensionsMat(A)[1];
  ncols := DimensionsMat(A)[2];

  for i in [1 .. nrows ] do
    for j in [1 .. ncols ] do
      if j <> i and tr.SNF[i][j] * one <> Zero(Integers mod 2^power) then
        Error("SNF not diagonal");
      fi;
    od;
  od;

  for i in [ 1 .. tr.rank ] do
    if tr.SNF[i][i] * one = Zero(Integers mod 2^power) then
      Error( " rank incorrect " );
    fi;
  od;

  for i in [ tr.rank + 1 .. Minimum(nrows, ncols) ] do
    if tr.SNF[i][i] * one <> Zero(Integers mod 2^power) then
      Error( " rank incorrect " );
    fi;
  od;

  if ((tr.row_t * one)
      * (A * one)
      * (tr.col_t * one)) <> tr.SNF * one then
    Error( "Transformation is not equal to SNF " );
  fi;

  if ((tr.row_t * one)
      * (tr.row_t_inverse * one)) <> IdentityMat(nrows) * one then
    Error( "Row transform inverse fail " );
  fi;
  Print("\npass");

end;

Read("../LocalSmithNormalForm.gap");
x := RandomMat(100, 90);;
tr := SNFTransform(x, 10: num_threads:=4);;
CheckAnswer(tr, x, 10);
