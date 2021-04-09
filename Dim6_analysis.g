Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("./dim_6/carat_ids.txt");


Print(Length(filtered_subgroups_ids), "\n");

i := 1;

id := filtered_subgroups_ids[i];
gens := Carat(id[1], id[2], id[3]);
pr := CoflasqueCover(gens);
Print(DimensionsMat(pr), "\n");
cr := CoflasqueResolution(pr, gens);
Print(DimensionsMat(cr.actionC[1]),"\n");

result := rec();
result.carat_id := id;
result.coflasque_cover_dim := DimensionsMat(cr.actionC[1])[1];

if result.coflasque_cover_dim <= 40 then
  res := ComputePhi(gens, cr);
  result.result := res;
else
  result.result := "NOT COMPUTED";
fi;

Print(result);
