Read("phi_computation.gap");
input_path := "/home/alexander/ch2bt/carat_tables/cryst6_filtered.txt";
Read(input_path);
output_path := "/home/alexander/ch2bt/carat_tables/cryst6_phi_result.txt";

for nr_q_cl in [ nr ] do
  for i in [ 1 .. Length(cryst6_filtered[nr_q_cl]) ] do
    gens := cryst6_filtered[nr_q_cl][i];
    result := ComputePhi(gens, CoflasqueResolution(gens));
    result.group_id := [nr_q_cl, i];
    Print(result);

    output_path := "/home/alexander/ch2bt/carat_tables/cryst6_phi_result.txt";
    AppendTo(output_path, result);
    AppendTo(output_path, ",");
    AppendTo(output_path, "\n");

  od;
od;
