Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("./dim_5_filtered_subgroups/carat_5_697_4_filtered_subgroups.txt");

output_path := "./dim_5_logs/carat_5_697_4_filtered_subgroups_result.txt";

#AppendTo(output_path, "carat_6_6129_4_filtered_results:=[\n");

gens := Carat(5, 697, 4);
pr := CoflasqueCover(gens);

for i in [1 .. Length(carat_5_697_4_filtered_subgroups)] do
  sub_gens := carat_5_697_4_filtered_subgroups[i];
  cr := CoflasqueResolution(pr, sub_gens);
  res := ComputePhi(sub_gens, cr: num_threads:=4);
  res.subgroup_num := i;
  AppendTo(output_path, res);
  AppendTo(output_path, ", ");
  Print( StringFormatted("\n{} / {} done\n", i, Length(carat_5_697_4_filtered_subgroups)) );
od;
AppendTo(output_path, "\n];");
