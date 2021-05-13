Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
result_path := "./dim_6/result.txt";
Read(result_path);

for i in [1 .. Length(result)] do
  if result[i].result = "NOT COMPUTED" then
    id := result[i].carat_id;
    Print("Doing ", id);
    Print("\n");
    gens := Carat(id[1], id[2], id[3]);
    pr := CoflasqueCover(gens);
    cr := CoflasqueResolution(pr, gens);
    res := ComputePhi(gens, cr: num_threads:=4);
    result[i].result := res;
    PrintTo(result_path, "result := ", result, ";");
  fi;
od;
