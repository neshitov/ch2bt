#######################################################################
# Code to compute Phi(G,M) for all filtered 2-subgroups G of GL(6,Z)
#######################################################################

Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("config.g");
Read("./dim_6/carat_ids.txt");


out_path := "./dim_6/result_copy.txt";

PrintTo(out_path, "result := [");

for i in [1 .. Length(filtered_subgroups_ids)] do
  id := filtered_subgroups_ids[i];
  gens := Carat(id[1], id[2], id[3]);
  Reset(GlobalMersenneTwister);
  Reset(GlobalRandomSource);

  pr := CoflasqueCover(gens);
  cr := CoflasqueResolution(pr, gens);

  result := rec();
  result.carat_id := id;
  result.coflasque_cover_dim := DimensionsMat(cr.actionC[1])[1];

  res := ComputePhi(gens, cr: num_threads:=n_threads);
  result.result := res;
  AppendTo(out_path, result, ",");
od;
AppendTo(out_path, "];");
