################################################
# Code to compute group Phi(G, M)
# written by Alexander Neshitov
################################################

LoadPackage("Gauss");
LoadPackage("LocalSNF");


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


Mod2FastInverse := function(m)
  # Return inverse of a finite order matrix
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


Mod2Image := function(m)
  # Return the image of the matrix map m
  local pivots;
  pivots := PositionsProperty(EchelonMat(m).heads, i-> i>0);
  return m{[1..DimensionsMat(m)[1]]}{pivots};
end;


CocycleRelation := function(gens, d, w)
  # Create map corresponding to the relation w.
  local n_gens, proj, result, action, proj_w_inverse, i, proj_, action_;
  n_gens := Length(gens);
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
    res := Z(2) * MatrixByBlockMatrix(res);
    return res;
  end;

  proj_w_inverse := function(i)
    if i > 0 then
      return proj(i);
    else
      return - action(i) * proj(-i);
    fi;
  end;

  result := proj_w_inverse( w[Length(w)] );
  for i in Reversed( w{[ 1.. Length(w) - 1 ]} ) do
    proj_ := proj_w_inverse(i);
    action_ := action(i);
    result:= proj_ + action_*result;
  od;
  return result;
end;


CocycleRelations := function(gens, d, relations)
  local cocycle_relations, relation;
  cocycle_relations := [];
  for relation in relations do
    cocycle_relations := Concatenation(cocycle_relations, CocycleRelation(gens, d, relation));
  od;
  return cocycle_relations;
end;


Z1Mod2 := function(gens, relations)
  # Construct the inclusion Z^1(G, M/2) -> M^gens/2
  local c_rel, result, d;
  d := DimensionsMat(gens[1])[1];
  c_rel := CocycleRelations(gens, d, relations);
  result := NullspaceMat(TransposedMat(c_rel));
  return TransposedMat(result);
end;


B1Mod2 := function(gens)
  # Construct the inclusion B^1(G,M/2) -> M^gens/2
  local result, d, b1_incl;
  d := DimensionsMat(gens[1])[1];
  b1_incl := BlockMatrix(List([1..Length(gens)], i -> [i, 1, Z(2) * ( gens[i] - IdentityMat(d) )]),
                         Length(gens), 1); # Inclusion B^1(G,M) -> M^gens
  b1_incl := Mod2Image(b1_incl);
  return b1_incl;
end;


B1 := function(gens)
  # Constructs the matrix of conoundary map \delta M -> M^gens
  local b, i, d;
  d := DimensionsMat(gens[1])[1];
  b := BlockMatrix(List([1..Length(gens)], i -> [i, 1, ( gens[i] - IdentityMat(d) )]),
                   Length(gens), 1); # Inclusion B^1(G,M) -> M^gens
  b := MatrixByBlockMatrix(b);
  return b;
end;


Z1Part := function(gens, power)
  # Finds generators of Z1 that do not belong to B1 and returns them modulo 2
  local b, vec, num_threads;

  num_threads := ValueOption("num_threads");
  if num_threads = fail then
    num_threads := 4;
  fi;

  b := B1(gens) mod ( 2 ^ (power + 1) );
  vec := SaturationVectors(b, power + 1: num_threads:=num_threads);
  vec := vec * One(Integers mod 2);
  return vec;
end;




ExtSquareIterator := function(d)
  # iterates over basis of exterior square in order 12 13 .. 1d, 23 24.. 2d, ...
  local result, i, j;
  result := [];
  for i in [1..d-1] do
    for j in [i+1..d] do
      Add(result, [i, j]);
    od;
  od;
  return result;
end;


ExteriorProduct := function(u, v)
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

MapK := function(res, num_gens)
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
    for k in [ 1 .. dim_P ] do
      Add(column, res.injection[k][i] * res.injection[k][j]);
    od;
    Add(result, column);
  od;
  result := TransposedMat(result);
  result := Z(2) * result;
  result := (Z(2) * res.surjection) * result;
  return result;
end;

MapKgens := function(res, num_gens)
  # For a resolution N -> P -> M given by matrix m constructs the matrix of the map
  # Lambda^2 N^gens -> M/2^gens
  local result;
  result := MapK(res, num_gens);
  result := BlockMatrix( List( [ 1.. num_gens ], i -> [ i, i, result ] ),
                         num_gens, num_gens ); # L2n^gens -> M/2^gens
  return result;
end;


MatConcat := function(m_list)
  local nrows;
  nrows := DimensionsMat(m_list[1])[1];
  return List([1..nrows], i -> Concatenation(List(m_list, x -> x[i])));
end;


ComputePhi := function(gens, cr)
  local relations, group_order, power, result, ZP2, ZM2, PM, BM2, _start, _end, L, L2,
        ZL, ZL2, mapk, num_threads;

  num_threads := ValueOption("num_threads");
  if num_threads = fail then
    num_threads := 4;
  fi;

  relations := MatrixRelations(gens);
  group_order := Order( Group( gens ) );
  power := PValuation(group_order, 2);

  if group_order <> 2 ^ power then
      Error("Not a 2 group");
  fi;

  result := rec();

  ZP2 := Z1Mod2(Z(2) * cr.actionP, relations); # Z^1(P/2) -> P/2^gens
  ZM2 := Z1Mod2(Z(2) * gens, relations); # Z^1(M/2) -> M/2^gens
  PM := BlockMatrix(List([1..Length(gens)], i -> [i, i, Z(2) * cr.surjection]),
                    Length(gens), Length(gens)); # P^gens -> M^gens;

  ZP2 := PM * ZP2; # Z1(P/2) -> M^gens/2
  BM2 := B1Mod2(gens); # B1(M/2) -> M^gens/2

  result.HM2_rank := RankMat( ZM2 ) - RankMat( BM2 );
  result.im_HP2_rank := RankMat( MatConcat( [BM2, ZP2] ) ) - RankMat( BM2 );

  if result.im_HP2_rank = result.HM2_rank then
    result.Phi_rank := 0;
    return result;
  fi;

  L := ExteriorSquare(cr.actionC); # Lambda^2(N)
  mapk := MapKgens(cr, Length(gens)); # L2N^gens -> M^gens/2

  ZL := Z1Part(L, power : num_threads:=num_threads);
  ZL := mapk * ZL;
  result.im_HL_rank := RankMat( MatConcat( [BM2, ZL] ) ) - RankMat(BM2);
  result.Phi_rank := RankMat( MatConcat( [BM2, ZL, ZP2] ) ) - RankMat( MatConcat( [BM2, ZP2] ) );

  return result;

end;
