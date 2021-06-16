Read("./dim_6_logs/carat_6_6129_2_filtered_subgroups_result.txt");
Print(Length(carat_6_6129_2_filtered_results));

nums := [];
for x in carat_6_6129_2_filtered_results do
  Add(nums, x.subgroup_num);
od;
suka := [];
for i in [1 .. Length(nums)] do
  if not nums[i] in suka then
    Add(suka, nums[i]);
  else
    Print(nums[i]);
    Print("\n");
  fi;
od;
Print(Length(suka));
