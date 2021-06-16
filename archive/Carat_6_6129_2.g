Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("./dim_6_filtered_subgroups/carat_6_6129_2_filtered_subgroups.txt");

output_path := "./dim_6_logs/carat_6_6129_2_filtered_subgroups_result.txt";

#AppendTo(output_path, "carat_6_6129_2_filtered_results:=[\n");

gens := Carat(6, 6129, 2);
pr := CoflasqueCover(gens);

N := Length(carat_6_6129_2_filtered_subgroups);
Print(N);

for i in [N+1 .. N] do
  sub_gens := carat_6_6129_2_filtered_subgroups[i];
  cr := CoflasqueResolution(pr, sub_gens);
  res := ComputePhi(sub_gens, cr: num_threads:=20);
  res.subgroup_num := i;
  AppendTo(output_path, res);
  AppendTo(output_path, ", ");
  Print( StringFormatted("\n{} / {} done\n", i, N) );
od;
AppendTo(output_path, "\n];");