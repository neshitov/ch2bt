Read("carat_6_6129_4_filtered_subgroups_nontrivial_carat_ids.txt");
Print(Length(nontrivial_carat_ids));
Print("\n");
unique_carat_ids := [];
unique_records := [];
for i in [1 .. Length(nontrivial_carat_ids)] do
  if not nontrivial_carat_ids[i].Carat_id in unique_carat_ids then
    Add(unique_carat_ids, nontrivial_carat_ids[i].Carat_id);
    Add(unique_records, nontrivial_carat_ids[i]);
  fi;
od;

Print(unique_carat_ids);

out_path := "carat_6_6129_4_filtered_subgroups_nontrivial_carat_ids_unique.txt";
PrintTo(out_path, "nontrivial_carat_ids:= ");
AppendTo(out_path, unique_records);
AppendTo(out_path, ";");
