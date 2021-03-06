Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");

Read("./dim_6_logs/carat_6_6129_4_filtered_subgroups_nontrivial.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_4_filtered_subgroups.txt");
Read("RatProbAlgTori/caratnumber.gap");

out_path := "./dim_6_logs/carat_6_6129_4_filtered_subgroups_nontrivial_carat_ids.txt";

total_nontrivial := Length(carat_6_6129_4_filtered_subgroups_nontrivial);
out := [];

for k in [1 .. total_nontrivial] do
  subgroup_nr := carat_6_6129_4_filtered_subgroups_nontrivial[k].subgroup_num;
  sub_gens := carat_6_6129_4_filtered_subgroups[subgroup_nr];
  G := Group(sub_gens);
  q_class := CaratQClass(G);


  N := Carat_nr_Z_classes(q_class[1], q_class[2]);
  Print("\n", N, " Z classes \n");

  for i in [1 .. N] do
    gens := Carat(q_class[1], q_class[2], i);
    Print("\r\n", [i, k], " out of ", [N, total_nontrivial], "\n");
    # construct a coflasque resolution:
    pr := CoflasqueCover(gens);
    cr := CoflasqueResolution(pr, gens);

    #Print(gens);
    #Print(B1Mod2(gens));

    # compute the group Phi(G, M):

    result := ComputePhi(gens, cr: num_threads:=4);
    result.Carat_id := [q_class[1], q_class[2], i];
    if result.Phi_rank <> 0 then
      Print(result);
      Print("\n");
      Add(out, result);
    fi;
  od;
  PrintTo(out_path, out);
od;
