Read("phi_computation.gap");
Read("/home/alexander/ch2bt/carat_tables/carat_6_6129_4_filtered_subgroups.txt");

output_path := "/home/alexander/ch2bt/dim_6_logs/carat_6_6129_4_filtered_subgroups_result.txt";

AppendTo(output_path, "cryst6_filtered_results:=[\n");

gens := Carat(6, 6129, 4);
#gens := Carat(4, 10, 1);
pr := CoflasqueCover(gens);

for i in [10000 .. Length(carat_6_6129_4_filtered_subgroups)] do
  sub_gens := carat_6_6129_4_filtered_subgroups[i];
  cr := ResolutionFromProjection(pr, sub_gens);
  res := ComputePhi(sub_gens, cr);
  res.subgroup_num := i;
  AppendTo(output_path, res);
  AppendTo(output_path, ", ");
  Print( StringFormatted("\n{} / {} done\n", i, Length(carat_6_6129_4_filtered_subgroups)) );
od;
AppendTo(output_path, "\n];");
