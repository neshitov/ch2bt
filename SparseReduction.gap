LoadPackage("Gauss");

GetColumn:= function(A, i)
  local sm;
  sm:= TransposedSparseMat(CertainColumns(A, [i]));
  return rec(indices:=sm!.indices[1],
             entries:=sm!.entries[1]);
end;

AddRowMultiple:= function(A, i, j, c)
  # Performs row operation R[i] <- R[i] + c * R[j].
  local mat;
  AddRow(A!.indices[j], c * A!.entries[j], A!.indices[i], A!.entries[i]);
end;

GeneralTwoByTwoTransform:= function(A, i, j, a, b, c, d)
  # Performs the row operation
  # R[i] <- aR[i] + bR[j]
  # R[j] <- cR[i] + dR[j]
  # where [[a,b],[c,d]] is unimodular
  local temp_i;
  temp_i := ShallowCopy(A!.entries[i]);
  AddRow(A!.indices[i], (a - 1) * A!.entries[i], A!.indices[i], A!.entries[i]);
  AddRow(A!.indices[j], b * A!.entries[j], A!.indices[i], A!.entries[i]);
  AddRow(A!.indices[j], (d - 1) * A!.entries[j], A!.indices[j], A!.entries[j]);
  AddRow(A!.indices[i], c * temp_i, A!.indices[j], A!.entries[j]);
end;

#SwitchElements:= function(A, i, j)
#  local temp;
#  if i<>j then
#    temp := A[i];
#    A[i] := A[j];
#    A[j] := temp;
#  fi;
#end;

SwitchRows:= function(A, i, j)
  local temp_entries, temp_indices;
  if i<>j then
    temp_entries := A!.entries[i];
    temp_indices := A!.indices[i];
    A!.entries[i] := A!.entries[j];
    A!.indices[i] := A!.indices[j];
    A!.entries[j] := temp_entries;
    A!.indices[j] := temp_indices;
  fi;
end;

#ReduceColumnGcd:= function(A, col, pivot_num, tr)
ReduceColumnGcd:= function(A, column_num, pivot_num, tr)
  #local counter1, counter2, entries_num, gcd, d;
  local counter1, counter2, entries_num, gcd, d, col, filter;

  col := GetColumn(A, column_num);
  filter := PositionsProperty(col.indices, x-> x>=pivot_num);
  col.indices := col.indices{filter};
  col.entries := col.entries{filter};
  entries_num := Length(col.indices);

  if entries_num = 0 then
    return rec(finished:=true, pivot_num:=pivot_num);
  fi;
  gcd := Gcd(col.entries);
  for counter1 in [1..entries_num] do
    if AbsInt(col.entries[counter1]) = gcd then # If there is gcd in the column
      for counter2 in [1..entries_num] do
        if counter2 <> counter1 then
          d:= QuoInt(col.entries[counter2], col.entries[counter1]);
          AddRowMultiple(A, col.indices[counter2], col.indices[counter1], - d);
          AddRowMultiple(tr, col.indices[counter2], col.indices[counter1], - d);
        fi;
      od;
      SwitchRows(A, pivot_num, col.indices[counter1]);
      SwitchRows(tr, pivot_num, col.indices[counter1]);
      return rec(finished:=true, pivot_num:=pivot_num+1);
    fi;
  od;
  return rec(finished:=false, pivot_num:=pivot_num);
end;

#BinaryGcd := function(A)
#  #~5 times faster than general gcd
#  local n;
#  n:=Length(A);
#  if n=0 then
#    return 0;
#  elif n=1 then
#    return A[1];
#  elif n=2 then
#    return GcdInt(A[1], A[2]);
#  fi;
#  return GcdInt(BinaryGcd(A{[1..QuoInt(n,2)]}),
#                BinaryGcd(A{[QuoInt(n,2)+1, n]}));
#end;

RandomPair := function(N)
  local a, b;
  a := Random([1..N-1]);
  b := Random([a+1..N]);
  return [a,b];
end;


ReduceColumn := function(A, column_num, pivot_num, tr)
  local col, result, i, j, k, pair, gcd, best_gcd, col_filter,
        filter, r1, r2, a, b, c, d, best_pair, entries_num, g;
  col := GetColumn(A, column_num);
  filter := PositionsProperty(col.indices, x-> x>=pivot_num);
  col.indices := col.indices{filter};
  col.entries := col.entries{filter};
  entries_num := Length(col.indices);

  #result := ReduceColumnGcd(A, col, pivot_num, tr);
  result := ReduceColumnGcd(A, column_num, pivot_num, tr);
  while not result.finished do
    # Randomly choose pair of entries with small gcd
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

      i := best_pair[1];
      j := best_pair[2];
      r1 := col.entries[i];
      r2 := col.entries[j];
      gcd := Gcdex(r1, r2);
      g := gcd.gcd;
      a := gcd.coef1;
      b := gcd.coef2;  # ar1+br2=g
      c := QuoInt(-r2, g);
      d := QuoInt(r1, g); # ad-bc=1
      GeneralTwoByTwoTransform(A, i, j, a, b, c, d);
      GeneralTwoByTwoTransform(tr, i, j, a, b, c, d);
      col_filter := PositionsProperty(col.entries, x-> x<>0);
      col.entries := col.entries{col_filter};
      col.indices := col.indices{col_filter};
    fi;
    result := ReduceColumnGcd(A, col, pivot_num, tr);
    Print("\nresult finished\n");
    Print(result.finished);
  od;
  return result.pivot_num;
end;

ReduceMatrix := function(A)
  local col_num, pivot_num, tr;
  pivot_num := 1;
  tr := SparseIdentityMatrix(Nrows(A), Integers);
  for col_num in [1..Ncols(A)] do
    Print("\ncol_num \n");
    Print(col_num);
    Print("\n");
    Display(A);

    pivot_num := ReduceColumn(A, col_num, pivot_num, tr);
  od;
  return tr;
end;
