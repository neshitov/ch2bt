Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("./dim_5/carat_ids.txt");

do_calculation := function()
  local id, gens, cr, pr, res;
  PrintTo("./dim_5/result.txt", "result := [");
  for id in filtered_subgroups_ids do
    gens := Carat(id[1], id[2], id[3]);
    pr := CoflasqueCover(gens);
    cr := CoflasqueResolution(pr, gens);
    res := ComputePhi(gens, cr);
    res.carat_id := id;
    AppendTo("./dim_5/result.txt", res, ",");
  od;
  AppendTo("./dim_5/result.txt", "];");
end;
