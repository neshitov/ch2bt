Read("./dim_6_logs/carat_6_6129_4_filtered_subgroups_result.txt");
out_path := "./dim_6_logs/carat_6_6129_4_filtered_subgroups_nontrivial.txt";
#non_trivial_list := [];
#for i in [1 .. Length(carat_6_6129_4_filtered_subgroups_result)] do
#  if carat_6_6129_4_filtered_subgroups_result[i].Phi_rank <> 0 then
#    Add(non_trivial_list, carat_6_6129_4_filtered_subgroups_result[i]);
#  fi;
#od;
#PrintTo(out_path, non_trivial_list);
Read(out_path);
