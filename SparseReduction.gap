LoadPackage("Gauss");

GetColumn:= function(A, i)
  local sm;
  #if i > Ncols(A) then
  #  Error(i);
  #fi;IndicesOfSparseMatrix
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


HNFTransform := function(A)
  # Returns a matrix M such that M*A=H in Hermite normal form
  local tr, rank, transform, inv_transform, sm, ring;
  ring := RingOfDefinition(A);
  # Hermite form of sparse matrix sm
  tr := EchelonMatTransformation(A);
  rank := Nrows(tr.vectors);
  transform := UnionOfRows(tr.coeffs, tr.relations);
  return rec(transform:=transform,
             H:=transform * A,
             pivots:=tr.heads);
end;


ColPermutation := function(_li)
  local current_col, perm, x;
  current_col := 1;
  perm := ();
  for x in _li do
    if x <> current_col then
      perm := perm * (current_col, x);
    fi;
    current_col := current_col + 1;
  od;
  return perm;
end;

ClearRow := function(H, col_t, i, p)
  # Records column operations to clear row i with pivot at position (i,i)
  local pivot, gcd_pval, j, c, tmp_indices, tmp_entries;

  pivot := GetEntry(H, i, i);
  gcd_pval := Minimum(List(EntriesOfSparseMatrix(H)[i], x -> PValuation(Int(x), p)));

  if PValuation(Int(pivot), p) <> gcd_pval then
    j := PositionProperty(EntriesOfSparseMatrix(H)[i], x -> PValuation(Int(x), p) = gcd_pval);
    AddColumnMultiple(H, i, IndicesOfSparseMatrix(H)[i][j], One(RingOfDefinition(H)));
    AddColumnMultiple(col_t, i, IndicesOfSparseMatrix(H)[i][j], One(RingOfDefinition(H)));
  fi;

  tmp_indices := ShallowCopy(IndicesOfSparseMatrix(H)[i]);
  tmp_entries := ShallowCopy(EntriesOfSparseMatrix(H)[i]);
  pivot := GetEntry(H, i, i);
  for j in [2..Length(EntriesOfSparseMatrix(H)[i])] do
    AddColumnMultiple(H, tmp_indices[j], i, - tmp_entries[j] / pivot);
    AddColumnMultiple(col_t, tmp_indices[j], i, - tmp_entries[j] / pivot);
  od;
end;

ClearColumn := function(H, row_t, i, p)
  # Records column operations to clear row i with pivot at position (i,i)
  local pivot, gcd_pval, j, c, tmp_indices, tmp_entries, column;

  pivot := GetEntry(H, i, i);
  column := GetColumn(H, i);
  gcd_pval := Minimum(List(column.entries, x -> PValuation(Int(x), p)));
  if PValuation(Int(pivot), p) <> gcd_pval then
    j := PositionProperty(column.entries, x -> PValuation(Int(x), p) = gcd_pval);
    AddRowMultiple(H, i, column.indices[j], One(RingOfDefinition(H)));
    AddRowMultiple(row_t, i, column.indices[j], One(RingOfDefinition(H)));
  fi;

  column := GetColumn(H, i);
  pivot := GetEntry(H, i, i);
  for j in [2..Length(column.entries)] do
    AddRowMultiple(H, column.indices[j], i, - column.entries[j] / pivot);
    AddRowMultiple(row_t, column.indices[j], i, - column.entries[j] / pivot);
  od;
end;

IsPivot := function(H, i)
  if GetEntry(H, i, i) = Zero(RingOfDefinition(H)) then
    return false;
  else
    return (Length(H!.indices[i]) = 1 and Length(GetColumn(H, i).indices) = 1);
  fi;
end;

IsProblem := function(M)
  local x, y;
  for x in EntriesOfSparseMatrix(M) do
    for y in x do
      if y = Zero(RingOfDefinition(M)) then
        return true;
      fi;
    od;
  od;
  return false;
end;

FixPivot := function(H, col_t, row_t, i, p)
  local next_row;
  #if IsProblem(H) then
  #  Error("Sparse fuckup");
  #fi;
  #while not IsPivot(H, i) do
  #  Print("pivot valuation ");
  #  Print(PValuation(Int(GetEntry(H, i, i)), p));
  #  Print("\n");
  #  ClearRow(H, col_t, i, p);
  #  ClearColumn(H, row_t, i, p);
  #od;
  if GetEntry(H, i, i) = Zero(RingOfDefinition(H)) then
    next_row := PositionProperty([i + 1..Nrows(H)], j -> Length(H!.indices[j]) > 0);
    if next_row = fail then
      return fail;
    fi;
    SwitchRows(H, i, next_row);
    SwitchRows(row_t, i, next_row);
  fi;

  while not IsPivot(H, i) do
    ClearRow(H, col_t, i, p);  # Will make pivot entry nonzero if it was zero
    ClearColumn(H, row_t, i, p);
  od;
end;

SNFTransform := function(A, p)
  local row_t, col_t, i;
  row_t := SparseIdentityMatrix(Nrows(A), RingOfDefinition(A));
  col_t := SparseIdentityMatrix(Ncols(A), RingOfDefinition(A));

  for i in [1 .. Maximum(Nrows(A), Ncols(A))] do
    Print(i);
    Print("\n");
    FixPivot(A, col_t, row_t, i, p);
  od;
end;

PreSNFTransform := function(A, p)
  # Returns Smith Normal form and transformation matrices over Z/p^N
  local hnf, pivot_columns, col_t, row_t, pm, H, rank,
        col, pvals, i, j, gcd_pval, c, pivot, coef, tmp_indices, tmp_entries, cr;
  hnf := HNFTransform(A);
  H := hnf.H;
  row_t := hnf.transform;
  Print("HNF done\n");
  pivot_columns := PositionsProperty([1..Ncols(A)], i -> hnf.pivots[i] > 0);
  rank := Length(pivot_columns);
  pm := PermutationMat(ColPermutation(pivot_columns), Ncols(A));
  pm := SparseMatrix(pm * One(RingOfDefinition(A)));
  col_t := pm;
  # Move all pivots to the left
  H := H * pm;
#  return H;
#end;
#  #i := 1;
  for i in [1..rank] do
    Print("do pivot ");
    Print(i);
    Print("\n");
    FixPivot(H, col_t, row_t, i, p);
    if ((Length(H!.indices[i]) <> 1
        or Length(GetColumn(H, i).indices) <> 1)
        or GetEntry(H, i, i) = 0) then
          Error("FixPivot fail");
    fi;
  od;

  #return rec(SNF:=H,
  #           row_transform:=row_t,
  #           col_transform:=col_t);
  return H;
end;
