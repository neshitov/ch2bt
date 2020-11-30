LoadPackage("Gauss");

#GetColumn:= function(A, i)
#  local sm;
#  sm:= TransposedSparseMat(CertainColumns(A, [i]));
#  return rec(indices:=sm!.indices[1],
#             entries:=sm!.entries[1]);
#end;


GetRow:= function(A, i)
  return rec(indices:=A!.indices[i],
             entries:=A!.entries[i]);
end;


AddRowMultiple:= function(A, i, j, c)
  # Performs row operation R[i] <- R[i] + c * R[j].
  local m;
  m := MultRow(A!.indices[j], A!.entries[j], c);
  AddRow(m.indices, m.entries,
         A!.indices[i], A!.entries[i]);
end;


AddColumnMultiple:= function(A, i, j, c)
  # Performs column operation C[i] <- C[i] + c * C[j]
  local k, x;
  for k in [ 1 .. Nrows(A) ] do
    x := GetEntry(A, k, j);
    AddToEntry(A, k, i, c * x);
  od;
end;


PermuteRows := function(A, perm)
  A!.indices := Permuted(A!.indices, perm);
  A!.entries := Permuted(A!.entries, perm);
end;


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


SwitchColumns := function(A, i, j)
  local t, nr_i, nr_j, tmp, insert_place, entry;
  for t in [ 1 .. Nrows(A) ] do
    nr_i := PositionSet(A!.indices[t], i);
    nr_j := PositionSet(A!.indices[t], j);

    if nr_i <> fail and nr_j <> fail then
      tmp := A!.entries[t][nr_i];
      A!.entries[t][nr_i] := A!.entries[t][nr_j];
      A!.entries[t][nr_j] := tmp;

    elif nr_i <> fail and nr_j = fail then
      Remove(A!.indices[t], nr_i);
      entry := Remove(A!.entries[t], nr_i);
      insert_place := PositionSorted(A!.indices[t], j);
      Add(A!.indices[t], j, insert_place);
      Add(A!.entries[t], entry, insert_place);

    elif nr_i = fail and nr_j <> fail then
      Remove(A!.indices[t], nr_j);
      entry := Remove(A!.entries[t], nr_j);
      insert_place := PositionSorted(A!.indices[t], i);
      Add(A!.indices[t], i, insert_place);
      Add(A!.entries[t], entry, insert_place);
    fi;
  od;
end;


#ClearRow := function(H, col_t, i, p)
#  # Records column operations to clear row i with pivot at position (i,i)
#  local pivot, gcd_pval, j, c, tmp_indices, tmp_entries, tmp_index;

#  pivot := GetEntry(H, i, i);
#  gcd_pval := Minimum(List(H!.entries[i], x -> PValuation(Int(x), p)));

#  if PValuation(Int(pivot), p) <> gcd_pval then
#    j := PositionProperty(H!.entries[i], x -> PValuation(Int(x), p) = gcd_pval);
#    AddColumnMultiple(H, i, tmp_index, One(RingOfDefinition(H)));
#    AddColumnMultiple(col_t, i, tmp_index, One(RingOfDefinition(H)));
#  fi;

#  tmp_entries := ShallowCopy(H!.indices[i]);
#  pivot := GetEntry(H, i, i);

#  for j in [2..Length(H!.entries[i])] do
#    AddColumnMultiple(H, tmp_indices[j], i, - tmp_entries[j] / pivot);
#    AddColumnMultiple(col_t, tmp_indices[j], i, - tmp_entries[j] / pivot);
#  od;
#end;


#ClearColumn := function(H, row_t, row_t_inverse, i, p)
  # Records column operations to clear row i with pivot at position (i,i)
#  local pivot, gcd_pval, j, c, tmp_indices, tmp_entries, column;

#  pivot := GetEntry(H, i, i);
#  column := GetColumn(H, i);
#  gcd_pval := Minimum(List(column.entries, x -> PValuation(Int(x), p)));
#  if PValuation(Int(pivot), p) <> gcd_pval then
#    j := PositionProperty(column.entries, x -> PValuation(Int(x), p) = gcd_pval);
#    AddRowMultiple(H, i, column.indices[j], One(RingOfDefinition(H)));
#    AddRowMultiple(row_t, i, column.indices[j], One(RingOfDefinition(H)));
#    AddRowMultiple(row_t_inverse, i, column.indices[j], - One(RingOfDefinition(H)));
#  fi;

#  column := GetColumn(H, i);
#  pivot := GetEntry(H, i, i);
#  for j in [2..Length(column.entries)] do
#    AddRowMultiple(H, column.indices[j], i, - column.entries[j] / pivot);
#    AddRowMultiple(row_t, column.indices[j], i, - column.entries[j] / pivot);
#    AddRowMultiple(row_t_inverse, column.indices[j], i, column.entries[j] / pivot);
#  od;
#end;


MultiplyRowByUnit := function(H, i, c)
  # Multiples i-th row by c in H and
  H!.entries[i] := c * H!.entries[i];
end;


MultiplyColumnByUnit := function(H, i, c)
  # Multiples i-th column by c in H and
  local nr, k;
  for k in [ 1 .. Nrows(H) ] do
    nr := PositionSet( H!.indices[k], i );
    if nr <> fail then
      H!.entries[k][nr] := H!.entries[k][nr] * c;
    fi;
  od;
end;

#IsPivot := function(H, i)
#  if GetEntry(H, i, i) = Zero(RingOfDefinition(H)) then
#    return false;
#  else
#    return (Length(H!.indices[i]) = 1 and Length(GetColumn(H, i).indices) = 1);
#  fi;
#end;


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
  for i in [ 1..Nrows(A) ] do
    PermuteSparseVector(A!.indices[i], A!.entries[i], perm, perm_left, perm_right);
  od;
end;


#FixPivot := function(H, col_t, row_t, row_t_inverse, i, p)
#  local next_row, pivot, unit;
#  if GetEntry(H, i, i) = Zero(RingOfDefinition(H)) then
#    next_row := PositionProperty([i + 1..Nrows(H)], j -> Length(H!.indices[j]) > 0);
#    if next_row = fail then
#      return false;
#    fi;
#    next_row := next_row + i;
#    SwitchRows(H, i, next_row);
#    SwitchRows(row_t, i, next_row);
#    SwitchRows(row_t_inverse, i, next_row);
#  fi;

#  while not IsPivot(H, i) do
#    ClearRow(H, col_t, i, p);  # Will make pivot entry nonzero if it was zero
#    ClearColumn(H, row_t, row_t_inverse, i, p);
#    Print("Do Clear row and col\n");
#  od;
#  pivot := GetEntry(H, i, i);
  #if pivot = Zero(RingOfDefinition(H)) then
  #  Error("Zero pivot");
  #fi;
#  unit := (p ^ PValuation(Int(pivot), p)) / pivot;
#  MultiplyRowByUnit(H, i, unit);
#  MultiplyRowByUnit(row_t, i, unit);
#  MultiplyRowByUnit(row_t_inverse, i, One(RingOfDefinition(H)) / unit);
#  return true;
#end;


ReduceModp := function(A, p)
  # Reduces matrix Z/2^n -> Z/2
  return List([1..DimensionsMat(A)[1]], i->List(A[i], x -> Int(x) * One(Z(2))));
end;

HNFTransform := function(A)
  # Returns a matrix M such that M*A=H in Hermite normal form
  local tr, rank, transform, transform_inverse, sm, ring;
  ring := RingOfDefinition(A);
  # Hermite form of sparse matrix sm
  tr := EchelonMatTransformation(A);
  rank := Nrows(tr.vectors);
  transform := UnionOfRows(tr.coeffs, tr.relations);
  #Print("HNF done");

  tr := EchelonMatTransformation(transform);
  transform_inverse := UnionOfRows(tr.coeffs, tr.relations);
  #Print("HNF tr_inverse done");

  return rec(transform:=transform,
             transform_inverse:=transform_inverse,
             H:=transform * A,
             pivots:=tr.heads);
end;

SparseInverse := function(A)
  return HNFTransform(A).transform;
end;


FindSmallestColumnEntry := function(A, t, j_t, p)
  # Finds entries with minimal valuations in row and column containing (t, j_t)
  local pivot_p_val, entry_p_val,
        min_p_val_in_column, min_p_val_in_column_index,
        i, nr;

  pivot_p_val := PValuation( Int( A!.entries[t][1] ), p );
  min_p_val_in_column := pivot_p_val;
  min_p_val_in_column_index := t; # row nr with minimal p val in column

  # find smallest p valuation in column
  for i in [ t + 1 .. Nrows(A)] do
    entry_p_val := PValuation( Int( GetEntry(A, i, j_t) ), p );
    if entry_p_val < min_p_val_in_column then
      min_p_val_in_column := entry_p_val;
      min_p_val_in_column_index := i;
    fi;
  od;

  return rec(min_p_val_in_column := min_p_val_in_column,
             min_p_val_in_column_index := min_p_val_in_column_index);
end;


FindSmallestRowColumnEntries := function(A, t, j_t, p)
  # Finds entries with minimal valuations in row and column containing (t, j_t)
  local pivot_p_val, entry_p_val,
        min_p_val_in_column, min_p_val_in_column_index,
        min_p_val_in_row, min_p_val_in_row_index,
        i, nr;

  pivot_p_val := PValuation( Int( A!.entries[t][1] ), p );
  min_p_val_in_column := pivot_p_val;
  min_p_val_in_column_index := t; # row nr with minimal p val in column
  min_p_val_in_row := pivot_p_val;
  min_p_val_in_row_index := j_t; # column nr with minimal p val in row

  # find smallest p valuation in column
  for i in [ t + 1 .. Nrows(A)] do
    entry_p_val := PValuation( Int( GetEntry(A, i, j_t) ), p );
    if entry_p_val < min_p_val_in_column then
      min_p_val_in_column := entry_p_val;
      min_p_val_in_column_index := i;
    fi;
  od;

  # find smallest p valuation in row
  for nr in [ 2 .. Length( A!.indices[t] ) ] do
    entry_p_val := PValuation( Int( A!.entries[t][nr] ), p );
    if entry_p_val < min_p_val_in_row then
      min_p_val_in_row := entry_p_val;
      min_p_val_in_row_index := A!.indices[t][nr];
    fi;
  od;

  return rec(min_p_val_in_column := min_p_val_in_column,
             min_p_val_in_column_index := min_p_val_in_column_index,
             min_p_val_in_row := min_p_val_in_row,
             min_p_val_in_row_index := min_p_val_in_row_index);
end;


#SNFTransformDestructive := function(A, row_t, row_t_inverse, col_t)
SNFTransformDestructive := function(A, row_t, row_t_inverse, col_t)
  # Transforms given matrix A into its Smith Normal Form
  # adding necessary row operations to row_t, their inverses to row_t_inverse
  # and column operations to col_t
  local p, char, heads,
        j, # column number
        i,
        k,
        pivot_i,
        pivot_j,
        nrows,
        t, i_t, j_t, min_v,
        coef,
        sm_rec,
        min_p_val_in_column, min_p_val_in_column_index,
        min_p_val_in_row, min_p_val_in_row_index,
        pivot_p_val,
        keep_going,
        tmp,
        rank_count,
        unit,
        tmp_indices,
        pivots,
        perm;

  #A := CopyMat(mat);
  rank_count := 1;
  nrows := Nrows(A);
  #row_t := SparseIdentityMatrix(nrows, A!.ring);
  #row_t_inverse := SparseIdentityMatrix(nrows, A!.ring);
  #col_t := SparseIdentityMatrix(Ncols(A), A!.ring);
  char := Characteristic( A!.ring );
  pivots := [];

  if char = 0 or Length( PrimePowersInt( char ) ) > 2 then
        Error( "only Z / p^n supported" );
  fi;
  p := PrimePowersInt( char )[1];

  # Do Hermite normal form reduction first

  for t in [ 1 .. Nrows(A) ] do

    # Choose j_t such that A[>=t,j_t] is the leftmost nonzero column in A[>=t,:]
    # and A[i_t, j_t] has the smallest p-valuation in the column A[:,j_t]
    j_t := infinity;
    min_v := infinity;
    i_t := 0;
    for i in [ t .. nrows ] do
      if Length(A!.indices[i]) > 0 then
        if A!.indices[i][1] < j_t then
          j_t := A!.indices[i][1];
          i_t := i;
          min_v := PValuation( Int( A!.entries[i][1] ), p );
        elif A!.indices[i][1] = j_t and PValuation( Int( A!.entries[i][1] ), p ) < min_v then
          i_t := i;
          min_v := PValuation( Int( A!.entries[i][1] ), p );
        fi;
      fi;
    od;

    if j_t = infinity then
      break;
    fi;

    # Move pivot to position to (t, j_t)
    if i_t <> t then
      SwitchRows(A, i_t, t);
      SwitchRows(row_t, i_t, t);
      SwitchColumns(row_t_inverse, i_t, t);
    fi;

    # get rid of all entries below the pivot in column j_t
    for i in [ t + 1 .. nrows] do
      if Length(A!.indices[i]) > 0 and A!.indices[i][1] = j_t then
        coef := - A!.entries[i][1] / A!.entries[t][1];
        AddRowMultiple( A, i, t, coef );
        AddRowMultiple( row_t, i, t, coef );
        AddColumnMultiple( row_t_inverse, t, i, - coef );
      fi;
    od;


    Print("hnf: do ");
    Print(t);
    Print("\n");
  od;


  for t in [ 1 .. Nrows(A) ] do

    # Choose j_t such that A[>=t,j_t] is the leftmost nonzero column in A[>=t,:]
    # and A[i_t, j_t] has the smallest p-valuation in the column A[:,j_t]
    j_t := infinity;
    min_v := infinity;
    i_t := 0;
    for i in [ t .. nrows ] do
      if Length(A!.indices[i]) > 0 then
        if A!.indices[i][1] < j_t then
          j_t := A!.indices[i][1];
          i_t := i;
          min_v := PValuation( Int( A!.entries[i][1] ), p );
        elif A!.indices[i][1] = j_t and PValuation( Int( A!.entries[i][1] ), p ) < min_v then
          i_t := i;
          min_v := PValuation( Int( A!.entries[i][1] ), p );
        fi;
      fi;
    od;

    if j_t = infinity then
      break;
    fi;
    Add( pivots, j_t );

    # Move pivot to position to (t, j_t)
    if i_t <> t then
      SwitchRows(A, i_t, t);
      SwitchRows(row_t, i_t, t);
      SwitchColumns(row_t_inverse, i_t, t);
    fi;

    # Make pivot divide all elements in its row and column
    keep_going := true;
    while keep_going do
      pivot_p_val := PValuation( Int( A!.entries[t][1] ), p );
      sm_rec := FindSmallestRowColumnEntries(A, t, j_t, p);
      min_p_val_in_column := sm_rec.min_p_val_in_column;
      min_p_val_in_column_index := sm_rec.min_p_val_in_column_index;
      min_p_val_in_row := sm_rec.min_p_val_in_row;
      min_p_val_in_row_index := sm_rec.min_p_val_in_row_index;

      if pivot_p_val = min_p_val_in_column and pivot_p_val = min_p_val_in_row then
        keep_going := false;
      elif min_p_val_in_column <= min_p_val_in_row then
        AddRowMultiple( A, t, min_p_val_in_column_index, One( A!.ring ) );
        AddRowMultiple( row_t, t, min_p_val_in_column_index, One( A!.ring ) );
        AddColumnMultiple( row_t_inverse, min_p_val_in_column_index, t, - One( A!.ring ) );
      else
        AddColumnMultiple( A, j_t, min_p_val_in_row_index, One( A!.ring ) );
        AddColumnMultiple( col_t, j_t, min_p_val_in_row_index, One( A!.ring ) );
      fi;
    od;

    # Make pivot power of p
    pivot_p_val := PValuation( Int( A!.entries[t][1] ), p );
    unit := ( p ^ pivot_p_val * One( A!.ring ) ) / A!.entries[t][1];
    MultiplyRowByUnit( A, t, unit );
    MultiplyRowByUnit( row_t, t, unit );
    MultiplyColumnByUnit( row_t_inverse, t, One( A!.ring ) / unit );

    # get rid of all entries in column j_t
    for i in [ t + 1 .. nrows] do
      if Length(A!.indices[i]) > 0 and A!.indices[i][1] = j_t then
        coef := - A!.entries[i][1] / A!.entries[t][1];
        AddRowMultiple( A, i, t, coef );
        AddRowMultiple( row_t, i, t, coef );
        AddColumnMultiple( row_t_inverse, t, i, - coef );
      fi;
    od;

    # get rid of all entries in row t
    j := 2;
    tmp_indices := ShallowCopy(A!.indices[t]);
    for i in tmp_indices do
      if i > j_t then
        coef := - GetEntry( A, t, i ) / GetEntry( A, t, j_t );
        AddColumnMultiple( A, i, j_t, coef );
        AddColumnMultiple( col_t, i, j_t, coef );
      fi;
    od;

  Print(t);
  Print(" done\n");
  rank_count := rank_count + 1;
  od;

  # move pivots to the main diagonal
  for t in [ 1 .. Length( pivots ) ] do
    if t <> pivots[t] then
      SwitchColumns( A, t, pivots[t] );
      SwitchColumns( col_t, t, pivots[t] );
    fi;
  od;

  # sort pivots

  perm := SortingPerm( List( [ 1 .. Length( pivots ) ],
                       k-> PValuation( Int( A!.entries[k][1] ), p) ) );
  #Print(perm);
  PermuteRows(A, perm);
  PermuteRows(row_t, perm);
  PermuteColumns(row_t_inverse, perm);

  PermuteColumns(A, perm);
  PermuteColumns(col_t, perm);

  return rank_count - 1;
end;


SNFTransform := function(mat)
  # Returns Smith Normal Form transform with row_transform, its inverse,
  # and column transform
  # Uses Hermite Normal form first to speed up computations.

  local hnf, col_t, row_t, row_t_inverse, rank, A;
  #hnf := HNFTransform(mat);
  #Print("HNF done\n");
  col_t := SparseIdentityMatrix(Ncols(mat), mat!.ring);
  row_t := SparseIdentityMatrix(Nrows(mat), mat!.ring);
  row_t_inverse := SparseIdentityMatrix(Nrows(mat), mat!.ring);
  #row_t := hnf.transform;
  #row_t_inverse := hnf.transform_inverse;
  #A := hnf.H;
  A := CopyMat( mat );
  rank := SNFTransformDestructive(A, row_t, row_t_inverse, col_t);
  return rec(SNF := A,
             rank := rank,
             row_t := row_t,
             row_t_inverse := row_t_inverse,
             col_t := col_t);
end;


PrintMatrix := function(file, mat)
  local nrows, ncols, power, i, j;
  if IsSparseMatrix(mat) then
    power := PrimePowersInt(Characteristic(mat!.ring))[2];
    nrows := Nrows(mat);
    ncols := Ncols(mat);
    PrintTo(file, nrows, " ", ncols, " ", power, "\n");
    for i in [1 .. nrows] do
      for j in [1 .. ncols] do
        AppendTo(file, Int(GetEntry(mat, i, j)));
        AppendTo(file, " ");
      od;
      AppendTo(file, "\n");
    od;
  else
    power := PrimePowersInt(Characteristic(Ring(mat[1][1])))[2];
    nrows := DimensionsMat(mat)[1];
    ncols := DimensionsMat(mat)[2];
    PrintTo(file, nrows, " ", ncols, " ", power, "\n");
    for i in [1 .. nrows] do
      for j in [1 .. ncols] do
        AppendTo(file, Int(mat[i][j]));
        AppendTo(file, " ");
      od;
      AppendTo(file, "\n");
    od;
  fi;
end;

ReadCResult := function(file)
  local tr, power;
  tr := ReadAsFunction(file)();
  power := tr.power;
  tr.SNF :=  tr.SNF * One(Integers mod 2^power);
  tr.row_t :=  tr.row_t * One(Integers mod 2^power);
  tr.row_t_inverse :=  tr.row_t_inverse * One(Integers mod 2^power);
  tr.col_t :=  tr.col_t * One(Integers mod 2^power);
  return tr;
end;

SNFTransform_c := function(mat)
  local tmp_dir, tmp_input_file, tmp_output_file, exec_file;

  tmp_dir := DirectoryTemporary();
  tmp_input_file := Filename(tmp_dir, "c_input");
  tmp_output_file := Filename(tmp_dir, "c_input_out");
  exec_file := Filename(DirectoryCurrent(), "c_code/a.out");
  PrintMatrix(tmp_input_file, mat);
  Exec(exec_file, tmp_input_file);
  return ReadCResult(tmp_output_file);
end;

CheckCAnswer := function(tr, A)
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

CheckAnswer := function(tr, A)
  local i;

  for i in [ 1 .. Nrows(A) ] do
    if tr.SNF!.indices[i] <> [] and tr.SNF!.indices[i] <> [i] then
      Error( " SNF not diagonal " );
    fi;
  od;

  for i in [ 1 .. tr.rank ] do
    if tr.SNF!.indices[i] <> [i] then
      Error( " rank incorrect " );
    fi;
  od;

  for i in [ tr.rank + 1 .. Nrows(A) ] do
    if tr.SNF!.indices[i] <> [] then
      Error( " rank incorrect " );
    fi;
  od;

  if (tr.row_t * A * tr.col_t) <> tr.SNF then
    Error( "Transformation is not equal to SNF " );
  fi;

  if (tr.row_t * tr.row_t_inverse) <> SparseIdentityMatrix(Nrows(A), A!.ring) then
    Error( "Row transform inverse fail " );
  fi;
  Print("pass");
end;
