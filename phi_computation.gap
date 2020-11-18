LoadPackage("Gauss");

Read("SparseReduction.gap");
carat_folder := "/home/alexander/ch2bt/carat_tables";
files := [Concatenation(carat_folder, "/cryst1.txt"),
          Concatenation(carat_folder, "/cryst2.txt"),
          Concatenation(carat_folder, "/cryst3.txt"),
          Concatenation(carat_folder, "/cryst4.txt"),
          Concatenation(carat_folder, "/cryst5.txt"),
          Concatenation(carat_folder, "/cryst6.txt"),
#          Concatenation(carat_folder, "/H1cryst.txt"),
#          Concatenation(carat_folder, "/caratchpol.txt"),
#          Concatenation(carat_folder, "/carat2crystcat.txt"),
#          Concatenation(carat_folder, "/caratnumber.gap"),
          "/home/alexander/ch2bt/FlabbyResolution.gap"
];
for file in files do
  Read(file);
od;
cryst := [cryst1, cryst2, cryst3, cryst4, cryst5, cryst6];
#Print(cryst1);

Carat := function(d, n, i)
  return cryst[d][n][i];
end;

Dual := function(mat_gens)
  # Return dual lattice to the lattice with action of group generated by mat_gens.
  local result, m;
  for m in mat_gens do
    Add(result, TransposedMat(Inverse(m)));
  od;
  return result;
end;

Mod2FastInverse := function(m)
  # Return inverse of a finite orded matrix
  local inverse, count;
  count := 0;
  inverse := Z(2) * IdentityMat(DimensionsMat(m)[1]);
  while (inverse * m <> Z(2) * IdentityMat(DimensionsMat(m)[1])) or (count > 100) do
    inverse := inverse * m;
    count := count + 1;
  od;
  if inverse * m <> Z(2) * IdentityMat(DimensionsMat(m)[1]) then
    Error("Fast inverse fail");
  fi;
  return inverse;
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

Mod2Image := function(m)
  # Return the image of the matrix map m
  local pivots;
  pivots := PositionsProperty(EchelonMat(m).heads, i-> i>0);
  return m{[1..DimensionsMat(m)[1]]}{pivots};
end;

LeftInverse := function(m)
  # Return a matrix a such that a*m = I if such a exists
  local emt, inv;
  emt := EchelonMatTransformation(m);
  inv := TransposedMat(emt.vectors) * emt.coeffs;
  if inv * m <> Z(2) * IdentityMat(DimensionsMat(m)[2]) then
    Error("Left inverse failure");
  fi;
  return inv;
end;

FlasqueResolution := function(g)
  # Flasque resolution (all matrices act on the left)
  local fr;
  fr:= FlabbyResolution(TransposedMatrixGroup(g));
  return rec(injection:=TransposedMat(fr.injection),
             surjection:=TransposedMat(fr.surjection),
             actionP:=TransposedMatrixGroup(fr.actionP),
             actionF:=TransposedMatrixGroup(fr.actionF));
end;

CoflasqueResolution := function(gens)
  # Coflasque resolution (all matrices act on the left)
  local d, fr;
  d := DimensionsMat(gens[1])[1];
  fr:= FlabbyResolution(Group(gens, IdentityMat(d)));
  return rec(injection:=fr.surjection,
             surjection:=fr.injection,
             actionP:=GeneratorsOfGroup(fr.actionP),
             actionC:=GeneratorsOfGroup(fr.actionF));
end;

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
end;

MatrixRelations := function(m)
  #  Get relations between matrices in the given list.
  local G, rel, rels, result;
  result := [];
  G := Group(m);
  G := Image(IsomorphismFpGroupByGenerators(G, m));
  rels := RelatorsOfFpGroup(G);
  for rel in rels do
    Add(result, TietzeWordAbstractWord(rel));
  od;
  return result;
end;

CocycleRelation := function(gens, d, w)
  # Create map corresponding to the relation w.
  local n_gens, proj, result, action, proj_w_inverse, i, proj_, action_, suka;
  #Print(" Cocycle relation start\n");
  n_gens := Length(gens);
  #result := SparseMatrix(Z(2) * NullMat(d, d * n_gens));
  #result := SparseZeroMatrix(d, d * n_gens);

  #Print("functions init start\n");
  action := function(i)
    if i > 0 then
      return gens[i];
    else
      return Mod2FastInverse(gens[-i]);
    fi;
  end;

  proj := function(i)
    local res;
    res := BlockMatrix([[1, i, IdentityMat(d)]], 1, n_gens);
    res := MatrixByBlockMatrix(res);
    #res := SparseMatrix(res);
    return res;
  end;

  proj_w_inverse := function(i)
    if i > 0 then
      return proj(i);
    else
      return - action(i) * proj(-i);
    fi;
  end;

  #Print("functions init end\n");
  result := proj_w_inverse(w[Length(w)]);
  for i in Reversed(w{[1..Length(w)-1]}) do
    proj_ := proj_w_inverse(i);
    action_ := action(i);
    suka:= action_*result;
    #Print("suka\n");
    #Print(suka);
    #Print("\nsuka\n");
    result:= proj_ + suka;
    #Print("addition done\n");
  od;
  #Print(" Cocycle relation end\n");
  return result;
end;

CocycleRelations := function(gens, d, relations)
  #local w;
  #w := MatrixRelations(gens);
  #return BlockMatrix(List([1..Length(relations)], i -> [i, 1, CocycleRelation(gens, d, relations[i])]),
  #                   Length(relations), 1);
  local cocycle_relations, relation;
  cocycle_relations := [];
  for relation in relations do
    cocycle_relations := Concatenation(cocycle_relations, CocycleRelation(gens, d, relation));
  od;
  return cocycle_relations;
end;

Mod2Z1 := function(gens, relations)
  # Construct the inclusion Z^1(G, M) -> M^gens G generated by matrices over Z/2.
  local c_rel, result, d;
  #sparse_gens := List(gens, g-> SparseMatrix(g));
  #Print("Mod2Z1: start construct relations\n");
  d := DimensionsMat(gens[1])[1];

  #c_rel := CocycleRelations(sparse_gens, d, relations);
  c_rel := CocycleRelations(gens, d, relations);
  #Print("Mod2Z1: relations constructed\n");
  result := NullspaceMat(TransposedMat(c_rel));
  #Print("Mod2Z1: Nullspace computed\n");
  return TransposedMat(result);
end;

Mod2B1 := function(gens, z1_incl)
  # Construct the inclusion B^1(G,M) -> Z^1(G,M) when gens are over Z/2.
  local result, d, b1_incl;
  d := DimensionsMat(gens[1])[1];
  b1_incl := BlockMatrix(List([1..Length(gens)], i -> [i, 1, gens[i] - IdentityMat(d)]),
                         Length(gens), 1); # Inclusion B^1(G,M) -> M^gens
  b1_incl := Mod2Image(b1_incl);
  return LeftInverse(z1_incl) * b1_incl;
end;


Testsuka := function(gens)
  local b, i, d;
  d := DimensionsMat(gens[1])[1];
  b := SparseMatrix(gens[1]);
  for i in [2..Length(gens)] do
    b := UnionOfRows(b, SparseMatrix(gens[i]));
  od;
  return b;
end;

B1 := function(gens)
  local b, i, d;
  d := DimensionsMat(gens[1])[1];
  b := SparseMatrix(gens[1] - IdentityMat(d));
  for i in [2..Length(gens)] do
    b := UnionOfRows(b, SparseMatrix(gens[i] - IdentityMat(d)));
  od;
  return b;
end;

Z1 := function(gens)
  # Construct the inclusion Z^1(G, M) -> M^gens G generated by integral matrices gens.
  local gen, b, smith_form, s, p, _rank, d, d_tr, tr, inverse_tr, rank, tr_rec;
  d := DimensionsMat(gens[1])[1];
  b := [];
  if gens=[] then
      return [];
  else
      for gen in gens do
          b := Concatenation(b, gen - IdentityMat(d));
      od;
      #Print("Z1: Start triangulization");

      #d_tr := TriangulizedIntegerMatTransform(b);
      #inverse_tr:= Inverse(d_tr.rowtrans);
      #rank:= d_tr.rank;

      b := SparseMatrix(b);
      tr_rec := ReduceMatrix(b);
      tr := ConvertSparseMatrixToMatrix(tr_rec.transform);
      inverse_tr := Inverse(tr);
      rank := tr_rec.rank;
      #d_tr :


      #Print("Z1: Triangulization done");
      #smith_form := SmithNormalFormIntegerMatTransforms(b);
      #s := smith_form.normal;
      #p := smith_form.rowtrans; # psq = smith_form
      #_rank := Rank(s);
      #return Inverse(p) * DiagonalMat(List([1.._rank], x-> 1));
      return inverse_tr * DiagonalMat(List([1..rank], x-> 1));
  fi;
end;

ExtSquareIterator :=function(d)
  # iterates over basis of exterior square in order 12 13 .. 1d, 23 24.. 2d,
  local result, i, j;
  result := [];
  for i in [1..d-1] do
    for j in [i+1..d] do
      Add(result, [i, j]);
    od;
  od;
  return result;
end;


ExteriorProduct :=function(u, v)
  local d, product, i, j, mi, pair, iterator;
  d := Length(u);
  product := [];
  iterator := ExtSquareIterator(d);
  for pair in iterator do
    i := pair[1];
    j := pair[2];
    mi := [[u[i], v[i]],
           [u[j], v[j]]];
    Add(product, Determinant(mi));
  od;
  return product;
end;

ExteriorSquare := function(gens)
  local square_gens, square_gen, column, i, j, m, m_columns, d, pair, iterator;
  square_gens := [];
  if Length(gens) > 0 then
    d := DimensionsMat(gens[1])[1];
    iterator := ExtSquareIterator(d);
    for m in gens do
      square_gen := [];
      m_columns := TransposedMat(m);
      for pair in iterator do
        i := pair[1];
        j := pair[2];
        column := ExteriorProduct(m_columns[i], m_columns[j]);
        Add(square_gen, column);
      od;
      square_gen := TransposedMat(square_gen);
      Add(square_gens, square_gen);
    od;
    return square_gens;
  else
    return [];
  fi;
end;

MapK := function(res)
  # For a resolution N -> P -> M given by matrix m constructs the matrix of the map
  # Lambda^2 N -> M/2
  local m_columns, dim_N, dim_P, result, column, i, j, k, iterator, pair;
  result := [];
  m_columns := TransposedMat(res.injection);
  dim_P := Length(res.injection);
  dim_N := Length(m_columns);
  iterator := ExtSquareIterator(dim_N);

  for pair in iterator do
    i := pair[1];
    j := pair[2];
    column := [];
    for k in [1..dim_P] do
      Add(column, res.injection[k][i] * res.injection[k][j]);
    od;
    Add(result, column);
  od;
  result := TransposedMat(result);
  result := Z(2) * result;
  result := (Z(2) * res.surjection) * result;
  return result;
end;

MatConcat := function(m_list)
  local nrows;
  nrows := DimensionsMat(m_list[1])[1];
  return List([1..nrows], i -> Concatenation(List(m_list, x -> x[i])));
end;

ComputePhi := function(gens, cr)
  local ZP, ZM, PM, ZPZM, BMZM, ZM_rank, HP_im, result, BM_rank, ZL2N, L2N,
        L2N2, mapk, ZL2N2, HL2N2_im, HL2N_im, relations, _start, _end;

  relations := MatrixRelations(gens);
  result := rec();

  _start := NanosecondsSinceEpoch();
  ZP := Mod2Z1(Z(2) * cr.actionP, relations); # Z1(P/2) -> P/2^gens
  ZM := Mod2Z1(Z(2) * gens, relations); # Z1(M/2) -> M/2^gens
  PM := BlockMatrix(List([1..Length(gens)], i -> [i, i, Z(2) * cr.surjection]),
                    Length(gens), Length(gens));
  ZPZM := LeftInverse(ZM) * PM * ZP; # Z1(P/2) -> Z1(M/2)
  ZPZM := Mod2Image(ZPZM); # image Z1(P/2) -> Z1(M/2)
  BMZM := Mod2B1(gens, ZM); # B1(M/2) -> Z1(M/2)
  ZM_rank := DimensionsMat(BMZM)[1];
  #if ZM_rank <> DimensionsMat(ZM)[2] then
  #  Error("ZM rank");
  #fi;
  BM_rank := RankMat(BMZM);
  result.H1M2_rank := ZM_rank - BM_rank;
  #HP_im := BlockMatrix([[1, 1, BMZM], [1, 2, ZPZM]], 1, 2);
  HP_im := MatConcat([BMZM, ZPZM]);
  result.im_H1P2_rank := RankMat(HP_im) - BM_rank;
  _end := NanosecondsSinceEpoch();
  Print("\n Mod 2 calculations done in\n");
  Print(Float((_end-_start)/(10^9)));
  Print("\n");

  _start := NanosecondsSinceEpoch();
  L2N := ExteriorSquare(cr.actionC);
  #L2N2 := Z(2) * L2N;
  mapk := MapK(cr);
  mapk := BlockMatrix(List([1..Length(gens)], i -> [i, i, mapk]),
                           Length(gens), Length(gens)); # L2n^gens -> M/2^gens
  mapk := LeftInverse(ZM) * mapk;
  _end := NanosecondsSinceEpoch();
  Print("\n Map k computed in\n");
  Print(Float((_end-_start)/(10^9)));
  Print("\n");
  #Print("ZL2N2 started\n");
  #ZL2N2 := Mod2Z1(L2N2);
  #Print("ZL2N2 computed\n");
  #HL2N2_im := BlockMatrix([[1, 1, BMZM], [1, 2, ZPZM], [1, 3, mapk * ZL2N2]], 1, 3);

  #Print("mapk ");
  #Print(DimensionsMat(mapk));
  #Print("\n BMZM ");
  #Print(DimensionsMat(BMZM));
  #Print("\n ZPZM ");
  #Print(DimensionsMat(ZPZM));

  #Print("\nZL2N2 started\n");
  #ZL2N2 := Mod2Z1(L2N2, relations);
  #Print("ZL2N2 computed\n");
  #HL2N2_im := MatConcat([BMZM, ZPZM, (mapk * ZL2N2)]);

  #result.im_H1L2N2_plus_H1P2_rank := RankMat(HL2N2_im) - BM_rank;

  #Print("\nZL2N started\n");
  _start := NanosecondsSinceEpoch();
  ZL2N := Z1(L2N);
  _end := NanosecondsSinceEpoch();
  Print("\n Z1 computed in\n");
  Print(Float((_end-_start)/(10^9)));
  Print("\n");

  #Print("ZL2N computed\n");
  HL2N_im := MatConcat([BMZM, ZPZM, (mapk * (Z(2) * ZL2N))]);

  result.im_H1L2N_plus_H1P2_rank := RankMat(HL2N_im) - BM_rank;
  #result.im_H1L2N2_plus_H1P2_rank := RankMat(HL2N2_im) - BM_rank;
  Print("\nResult:\n");
  return result;
end;


#m := [[1,3],[4,5]];
#Print(m);
#Print(MapK(m));
#G := CaratMatGroupZClass(5,697,4);
#cr := CoflasqueResolution(G);

#CheckEquivariance(G, cr);
