GeneralTwoByTwoTransform:= function(A, i, j, a, b, c, d)
  # Performs the row operation
  # R[i] <- aR[i] + bR[j]
  # R[j] <- cR[i] + dR[j]
  # where [[a,b],[c,d]] is unimodular
  local indices_i, entries_i, indices_j, entries_j;
  indices_i := ShallowCopy(IndicesOfSparseMatrix(A)[i]);
  entries_i := ShallowCopy(EntriesOfSparseMatrix(A)[i]);
  indices_j := ShallowCopy(IndicesOfSparseMatrix(A)[j]);
  entries_j := ShallowCopy(EntriesOfSparseMatrix(A)[j]);

  AddRow(indices_i, (a - 1) * entries_i, IndicesOfSparseMatrix(A)[i], EntriesOfSparseMatrix(A)[i]);
  AddRow(indices_j, b * entries_j, IndicesOfSparseMatrix(A)[i], EntriesOfSparseMatrix(A)[i]);
  AddRow(indices_i, c * entries_i, IndicesOfSparseMatrix(A)[j], EntriesOfSparseMatrix(A)[j]);
  AddRow(indices_j, (d - 1) * entries_j, IndicesOfSparseMatrix(A)[j], EntriesOfSparseMatrix(A)[j]);

end;

#ClearRowTail:= function(A, i)
#  # Clear all entries in the row except the last first one.
#  A!.indices[i] := [A!.indices[i][1]];
#  A!.entries[i] := [A!.entries[i][1]];
#end;

#SwitchElements:= function(A, i, j)
#  local temp;
#  if i<>j then
#    temp := A[i];
#    A[i] := A[j];
#    A[j] := temp;
#  fi;
#end;

#SwitchRows:= function(A, i, j)
#  local temp_entries, temp_indices;
#  if i<>j then
#    temp_entries := EntriesOfSparseMatrix(A)[i];
#    temp_indices := IndicesOfSparseMatrix(A)[i];
#    A!.entries[i] := EntriesOfSparseMatrix(A)[j];
#    A!.indices[i] := IndicesOfSparseMatrix(A)[j];
#    A!.entries[j] := temp_entries;
#    A!.indices[j] := temp_indices;
#  fi;
#end;

#ReduceColumnGcd:= function(A, col, pivot_num, tr)
#ReduceColumnGcd:= function(A, column_num, pivot_num, tr)
  #Reduces column if it has gcd entry.
#  local counter1, counter2, entries_num, gcd, d, col, filter;

#  col := GetColumn(A, column_num);
#  filter := PositionsProperty(col.indices, x-> x>=pivot_num);
#  col.indices := col.indices{filter};
#  col.entries := col.entries{filter};
#  entries_num := Length(col.indices);

#  if entries_num = 0 then
#    return rec(finished:=true, pivot_num:=pivot_num);
#  fi;
#  gcd := Gcd(col.entries);
#  for counter1 in [1..entries_num] do
#    if AbsInt(col.entries[counter1]) = gcd then # If there is gcd in the column
#      for counter2 in [1..entries_num] do
#        if counter2 <> counter1 then
#          d:= QuoInt(col.entries[counter2], col.entries[counter1]);
#          AddRowMultiple(A, col.indices[counter2], col.indices[counter1], - d);
#          AddRowMultiple(tr, col.indices[counter2], col.indices[counter1], - d);
#        fi;
#      od;
#      SwitchRows(A, pivot_num, col.indices[counter1]);
#      SwitchRows(tr, pivot_num, col.indices[counter1]);
#      Print("pivot num ");
#      Print(pivot_num);
#      ClearRowTail(A, pivot_num);  # To save space;
#      return rec(finished:=true, pivot_num:=pivot_num+1);
#    fi;
#  od;
#  return rec(finished:=false, pivot_num:=pivot_num);
#end;

RandomPair := function(N)
  local a, b;
  a := Random([1..N-1]);
  b := Random([a+1..N]);
  return [a,b];
end;


FixTwoRows := function(A, column_num, i, j, tr)
  # Reduce rows i and j to its gcd in column col_num
  local r1, r2, gcd_, g, a, b, c, d;

  r1 := GetEntry(A, i, column_num);
  r2 := GetEntry(A, j, column_num);
  gcd_ := Gcdex(r1, r2);
  g := gcd_.gcd;
  a := gcd_.coeff1;
  b := gcd_.coeff2;  # ar1+br2=g
  c := QuoInt(-r2, g);
  d := QuoInt(r1, g); # ad-bc=1
  GeneralTwoByTwoTransform(A, i, j, a, b, c, d);
  GeneralTwoByTwoTransform(tr, i, j, a, b, c, d);
end;


ReduceColumn := function(A, column_num, pivot_num, tr)
  local col, result, i, j, i_num, j_num, k, pair, gcd, gcd_, best_gcd, col_filter,
        filter, r1, r2, a, b, c, d, best_pair, entries_num, g;

  #result := ReduceColumnGcd(A, col, pivot_num, tr);
  result := ReduceColumnGcd(A, column_num, pivot_num, tr);
  while not result.finished do
    # Randomly choose pair of entries with small gcd
    col := GetColumn(A, column_num);
    filter := PositionsProperty(col.indices, x-> x>=pivot_num);
    col.indices := col.indices{filter};
    col.entries := col.entries{filter};
    filter := PositionsProperty(col.entries, x-> x<>0);
    col.indices := col.indices{filter};
    col.entries := col.entries{filter};
    entries_num := Length(col.indices);

    if Length(col.indices) > 1 then
      best_pair := RandomPair(Length(col.indices));
      best_gcd := GcdInt(col.entries[best_pair[1]], col.entries[best_pair[2]]);
      for k in [1..3] do
        pair := RandomPair(entries_num);
        gcd := GcdInt(col.entries[pair[1]], col.entries[pair[2]]);
        if gcd < best_gcd then
          best_pair := pair;
          best_gcd := gcd;
        fi;
      od;

      i_num := best_pair[1];
      j_num := best_pair[2];
      i := col.indices[i_num];
      j := col.indices[j_num];
      FixTwoRows(A, column_num, i, j, tr);
    fi;
    result := ReduceColumnGcd(A, column_num, pivot_num, tr);
  od;
  return result.pivot_num;
end;

#NonzeroInColumns :=function(m)
#  local mt;
#  mt := TransposedSparseMat(m);
#  return List([1..Ncols(m)], i-> Length(mt!.indices[i]));
#end;

ReduceMatrix := function(A)
  # Reduces SparseMatrix A to upper triangular in-place.
  # Returns the sparse row transform matrix
  local col_num, pivot_num, tr;
  pivot_num := 1;
  tr := SparseIdentityMatrix(Nrows(A), Integers);
  for col_num in [1..Ncols(A)] do
    #Print("\ncol_num \n");
    #Print(col_num);
    #Print("\n");
    #Display(A);

    pivot_num := ReduceColumn(A, col_num, pivot_num, tr);
    Print(col_num);
    Print(" ");
    Print(pivot_num);
    Print("\n");
  od;
  return rec(transform:= tr,
             rank:=pivot_num - 1);
end;

ElementaryMatrix := function(n, i, j, c, ring)
  # Elementary matrix M:
  # _*M performs column operation C[j] <- c*C[i] + C[j]
  # M*_ performs row operation R[i] <- R[i] + c*R[j]
  local res;
  res := SparseIdentityMatrix(n, ring);
  SetEntry(res, i, j, c);
  return res;
end;

MatProduct := function(a,b)
  local dim_a, dim_b;
  dim_a := DimensionsMat(a);
  dim_b := DimensionsMat(b);
  if dim_a[2] = dim_b[1] then
    return a * b;
  else
    return Error(Concatenation([dim_a, dim_b]));
  fi;
end;
