Read("filter_tools.gap")

input_path := "/home/alexander/toric_ZG_module/ZG_module/carat_tables/cryst6.txt";
Read(input_path);

to_keep := [];
for q_cl in cryst6 do
  z_classes_to_keep := [];
  for z_cl in q_cl do
    if KeepGroup(6, z_cl) then
      Add(z_classes_to_keep, z_cl);
    fi;
  od;
  if Length(z_classes_to_keep) > 0 then
    Add(to_keep, z_classes_to_keep);
  fi;
od;

output_path := "/home/alexander/toric_ZG_module/ZG_module/carat_tables/cryst6_filtered.txt";
PrintTo(output_path, "cryst6_filtered:=");
AppendTo(output_path, to_keep);
AppendTo(output_path, ";");
