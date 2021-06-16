Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("./dim_6/carat_ids.txt");

do_calculation := function()
  local i, id, gens, cr, pr, res;
  PrintTo("./dim_6/filtered_subgroups_result.txt", "result := [");
  for i in [1 .. Length(filtered_subgroups_ids)] do
    id := filtered_subgroups_ids[i];
    gens := Carat(id[1], id[2], id[3]);
    pr := CoflasqueCover(gens);
    cr := CoflasqueResolution(pr, gens);
    res := ComputePhi(gens, cr);
    res.carat_id := id;
    AppendTo("./dim_6/filtered_subgroups_result.txt", res, ",");
    Print("\n", i, " done of ", Length(filtered_subgroups_ids), "\n" );
  od;
  AppendTo("./dim_6/filtered_subgroups_result.txt", "];");
end;

do_calculation();
