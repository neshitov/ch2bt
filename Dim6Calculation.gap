Read("phi_computation.gap");
input_path := "/home/alexander/ch2bt/carat_tables/cryst6_filtered_ids.txt";
Read(input_path);
output_path := "/home/alexander/ch2bt/dim_6_logs/cryst6_filtered_phi_result.txt";

PrintTo(output_path, "cryst6_filtered_results:=[\n");
for id_num in [1..Length(cryst6_filtered_ids)] do
  id := cryst6_filtered_ids[id_num];
  gens := cryst[id[1]][id[2]][id[3]];
  cr := CoflasqueResolution(gens);
  result := ComputePhi(gens, cr: num_threads:=6);
  result.carat_id := id;
  Print( StringFormatted("\n{} / {} done\n", id_num, Length(cryst6_filtered_ids)) );

  AppendTo(output_path, result);
  AppendTo(output_path, ",");
  AppendTo(output_path, "\n");
od;
AppendTo(output_path, "\n];");
