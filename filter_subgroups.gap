Read("filter_tools.gap");
input_path := "/home/alexander/ch2bt/carat_tables/cryst6.txt";
Read(input_path);

ids_to_keep := [];
#for q_cl in cryst6 do
for q_cl_num in [1 .. Length(cryst6)] do
#for q_cl_num in [9 .. 9] do
  for z_cl_num in [1 .. Length(cryst6[q_cl_num])] do
  #for z_cl_num in [1 .. 6] do
    #Print(z_cl_num);
    if KeepGroup(6, cryst6[q_cl_num][z_cl_num]) then
      Add(ids_to_keep, [6, q_cl_num, z_cl_num]);
    fi;
  od;
  Print(Concatenation("\r", StringFormatted("Q_cl  {} / {} done", q_cl_num, Length(cryst6))));
od;

output_path := "/home/alexander/ch2bt/carat_tables/cryst6_filtered_ids.txt";
PrintTo(output_path, "cryst6_filtered_ids:=");
AppendTo(output_path, ids_to_keep);
AppendTo(output_path, ";");
