Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("./dim_6/carat_ids.txt");


out_path := "./dim_6/filtered_subgroups_result_small_rank.txt";

PrintTo(out_path, "result := [");

for i in [1 .. Length(filtered_subgroups_ids)] do

  id := filtered_subgroups_ids[i];
  gens := Carat(id[1], id[2], id[3]);
  pr := CoflasqueCover(gens);
  Print(DimensionsMat(pr), "\n");
  cr := CoflasqueResolution(pr, gens);
  Print(DimensionsMat(cr.actionC[1]),"\n");

  result := rec();
  result.carat_id := id;
  result.coflasque_cover_dim := DimensionsMat(cr.actionC[1])[1];

  if result.coflasque_cover_dim <= 85 then
    res := ComputePhi(gens, cr: num_threads:=20);
    result.result := res;
  else
    result.result := "NOT COMPUTED";
  fi;

  AppendTo(out_path, result, ",");
  Print("\n", i, " done of ", Length(filtered_subgroups_ids), "\n" );
od;
AppendTo(out_path, "];");
