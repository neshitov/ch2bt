LoadPackage("Gauss");

GetColumn:= function(A, i)
  local sm;
  sm:= TransposedSparseMat(CertainColumns(A, [i]));
  return rec(indices:=IndicesOfSparseMatrix(sm)[1],
             entries:=EntriesOfSparseMatrix(sm)[1]);
end;


GetRow:= function(A, i)
  local sm;
  sm:=CertainRows(A, [i]);
  return rec(indices:=IndicesOfSparseMatrix(sm)[1],
             entries:=EntriesOfSparseMatrix(sm)[1]);
end;


AddRowMultiple:= function(A, i, j, c)
  # Performs row operation R[i] <- R[i] + c * R[j].
  local row, k;
  row := GetRow(A, j);
  for k in [1..Length(row.indices)] do
    AddToEntry(A, i, row.indices[k], c * row.entries[k]);
  od;
  #AddRow(IndicesOfSparseMatrix(A)[j], c * EntriesOfSparseMatrix(A)[j],
  #       IndicesOfSparseMatrix(A)[i], EntriesOfSparseMatrix(A)[i]);
end;


AddColumnMultiple:= function(A, i, j, c)
  # Performs column operation C[i] <- C[i] + c * C[j]
  local col, k;
  col := GetColumn(A, j);
  for k in [1..Length(col.indices)] do
    AddToEntry(A, col.indices[k], i, c * col.entries[k]);
  od;
end;


PermuteRows := function(A, perm)
  A!.indices := Permuted(A!.indices, perm);
  A!.entries := Permuted(A!.entries, perm);
end;


SwitchRows:= function(A, i, j)
  local temp_entries, temp_indices;
  if i<>j then
    temp_entries := EntriesOfSparseMatrix(A)[i];
    temp_indices := IndicesOfSparseMatrix(A)[i];
    A!.entries[i] := EntriesOfSparseMatrix(A)[j];
    A!.indices[i] := IndicesOfSparseMatrix(A)[j];
    A!.entries[j] := temp_entries;
    A!.indices[j] := temp_indices;
  fi;
end;


ClearRow := function(H, col_t, i, p)
  # Records column operations to clear row i with pivot at position (i,i)
  local pivot, gcd_pval, j, c, tmp_indices, tmp_entries, tmp_index;

  pivot := GetEntry(H, i, i);
  gcd_pval := Minimum(List(EntriesOfSparseMatrix(H)[i], x -> PValuation(Int(x), p)));

  if PValuation(Int(pivot), p) <> gcd_pval then
    j := PositionProperty(EntriesOfSparseMatrix(H)[i], x -> PValuation(Int(x), p) = gcd_pval);
    tmp_index := IndicesOfSparseMatrix(H)[i][j];
    AddColumnMultiple(H, i, tmp_index, One(RingOfDefinition(H)));
    AddColumnMultiple(col_t, i, tmp_index, One(RingOfDefinition(H)));
  fi;

  tmp_indices := ShallowCopy(IndicesOfSparseMatrix(H)[i]);
  tmp_entries := ShallowCopy(EntriesOfSparseMatrix(H)[i]);
  pivot := GetEntry(H, i, i);

  for j in [2..Length(EntriesOfSparseMatrix(H)[i])] do
    AddColumnMultiple(H, tmp_indices[j], i, - tmp_entries[j] / pivot);
    AddColumnMultiple(col_t, tmp_indices[j], i, - tmp_entries[j] / pivot);
  od;
end;


ClearColumn := function(H, row_t, row_t_inverse, i, p)
  # Records column operations to clear row i with pivot at position (i,i)
  local pivot, gcd_pval, j, c, tmp_indices, tmp_entries, column;

  pivot := GetEntry(H, i, i);
  column := GetColumn(H, i);
  gcd_pval := Minimum(List(column.entries, x -> PValuation(Int(x), p)));
  if PValuation(Int(pivot), p) <> gcd_pval then
    j := PositionProperty(column.entries, x -> PValuation(Int(x), p) = gcd_pval);
    AddRowMultiple(H, i, column.indices[j], One(RingOfDefinition(H)));
    AddRowMultiple(row_t, i, column.indices[j], One(RingOfDefinition(H)));
    AddRowMultiple(row_t_inverse, i, column.indices[j], - One(RingOfDefinition(H)));
  fi;

  column := GetColumn(H, i);
  pivot := GetEntry(H, i, i);
  for j in [2..Length(column.entries)] do
    AddRowMultiple(H, column.indices[j], i, - column.entries[j] / pivot);
    AddRowMultiple(row_t, column.indices[j], i, - column.entries[j] / pivot);
    AddRowMultiple(row_t_inverse, column.indices[j], i, column.entries[j] / pivot);
  od;
end;


MultiplyRowByUnit := function(H, i, c)
  # Multiples i-th row by c in H and
  H!.entries[i] := c * H!.entries[i];
end;


IsPivot := function(H, i)
  if GetEntry(H, i, i) = Zero(RingOfDefinition(H)) then
    return false;
  else
    return (Length(H!.indices[i]) = 1 and Length(GetColumn(H, i).indices) = 1);
  fi;
end;


PermuteSparseVector := function(indices, entries, perm, perm_left, perm_right)
  local i;
  if Length(indices) > 0 then
    if not (indices[1] > perm_right or indices[Length(indices)] < perm_left) then
      for i in [1..Length(indices)] do
        indices[i] := indices[i] ^ perm;
      od;
      SortParallel(indices, entries);
    fi;
  fi;
end;


PermuteColumns := function(A, perm)
  local i, perm_left, perm_right;
  perm_left := SmallestMovedPoint(perm);
  perm_right := LargestMovedPoint(perm);
  for i in [1..Nrows(A)] do
    PermuteSparseVector(A!.indices[i], A!.entries[i], perm, perm_left, perm_right);
  od;
end;


FixPivot := function(H, col_t, row_t, row_t_inverse, i, p)
  local next_row, pivot, unit;
  if GetEntry(H, i, i) = Zero(RingOfDefinition(H)) then
    next_row := PositionProperty([i + 1..Nrows(H)], j -> Length(H!.indices[j]) > 0);
    if next_row = fail then
      return false;
    fi;
    next_row := next_row + i;
    SwitchRows(H, i, next_row);
    SwitchRows(row_t, i, next_row);
    SwitchRows(row_t_inverse, i, next_row);
  fi;

  while not IsPivot(H, i) do
    ClearRow(H, col_t, i, p);  # Will make pivot entry nonzero if it was zero
    ClearColumn(H, row_t, row_t_inverse, i, p);
  od;
  pivot := GetEntry(H, i, i);
  #if pivot = Zero(RingOfDefinition(H)) then
  #  Error("Zero pivot");
  #fi;
  unit := (p ^ PValuation(Int(pivot), p)) / pivot;
  MultiplyRowByUnit(H, i, unit);
  MultiplyRowByUnit(row_t, i, unit);
  MultiplyRowByUnit(row_t_inverse, i, One(RingOfDefinition(H)) / unit);
  return true;
end;


ReduceModp := function(A, p)
  # Reduces matrix Z/2^n -> Z/2
  return List([1..DimensionsMat(A)[1]], i->List(A[i], x -> Int(x) * One(Z(2))));
end;


SNFTransform := function(A, p)
  local row_t, row_t_inverse, col_t, col_t_inverse, i, perm;
  row_t := SparseIdentityMatrix(Nrows(A), RingOfDefinition(A));
  row_t_inverse := SparseIdentityMatrix(Nrows(A), RingOfDefinition(A));
  col_t := SparseIdentityMatrix(Ncols(A), RingOfDefinition(A));

  i := 1;
  while FixPivot(A, col_t, row_t, row_t_inverse, i, p) do
    Print(i);
    Print(" done\n");
    i := i + 1;
  od;

  perm := SortingPerm(List([1..i], x-> PValuation(Int(GetEntry(A, x, x)), p)));
  PermuteRows(A, perm);
  PermuteRows(row_t, perm);
  PermuteColumns(row_t_inverse, Inverse(perm));

  PermuteColumns(A, perm);
  PermuteColumns(col_t, perm);

  return rec(rank := i - 1,
             row_transform := row_t,
             row_transform_inverse := row_t_inverse,
             column_transform := col_t);
end;
