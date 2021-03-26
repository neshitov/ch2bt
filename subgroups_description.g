Read("./dim_6_logs/carat_6_6129_4_filtered_subgroups_nontrivial_carat_ids.txt");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("ComputePhi.g");
Read("RatProbAlgTori/caratnumber.gap");

results := [];

for i in [1 .. Length(nontrivial_carat_ids)] do

  r := nontrivial_carat_ids[i];
  gens := Carat(r.Carat_id[1], r.Carat_id[2], r.Carat_id[3]);
  G := Group(gens);
  fpg := Image(IsomorphismFpGroupByGenerators(G, gens));

  result := rec();
  result.carat_id := r.Carat_id;
  result.phi_rank := r.Phi_rank;
  result.order := Order(G);
  result.generators := GeneratorsOfGroup(fpg);
  result.relators := RelatorsOfFpGroup(fpg);

  Add(results, result);
od;
out_path := "./dim_6_logs/carat_6_6129_4_filtered_subgroups_nontrivial_description.txt";
PrintTo(out_path, results);

gens := Carat(6,2377,10);
pr := CoflasqueCover(gens);
cr := CoflasqueResolution(pr, gens);

Print(ComputePhi(gens, cr));

Read("./dim_6_logs/carat_6_6129_4_filtered_subgroups_nontrivial.txt");
Read("dim_6_filtered_subgroups/carat_6_6129_2_filtered_subgroups.txt");

total := Length(carat_6_6129_2_filtered_subgroups);
out := [];

Print("\n");
for k in [1 .. total] do
  sub_gens := carat_6_6129_2_filtered_subgroups[k];
  G := Group(sub_gens);
  q_class := CaratQClass(G);
  Print(k);
  Print("\r");
  if q_class = [6, 2377] then
    Print("\n");
    Print(k);
    Print("\n");
    Print(q_class);
    Print("\n");
    Print("\n");
  fi;
od;
