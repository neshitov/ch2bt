#######################################################################################
# Looking for (6,2377,10)
#######################################################################################

Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");

Read("dim_6_filtered_subgroups/carat_6_6129_4_filtered_subgroups.txt");
Read("dim_6_filtered_subgroups/carat_6_6129_2_filtered_subgroups.txt");
Read("./RatProbAlgTori/caratnumber.gap");


total := Length(carat_6_6129_4_filtered_subgroups);

for k in [1 .. total] do
  sub_gens := carat_6_6129_4_filtered_subgroups[k];
  G := Group(sub_gens);
  q_class := CaratQClass(G);
  if (q_class[1] = 6) and (q_class[2] = 2377) then
    Print("\n 6_6129_4 filtered_subgroup nr = ", k, "\n");
    Print("\n", q_class, "\n");
  fi;
od;

total := Length(carat_6_6129_2_filtered_subgroups);
out := [];

for k in [1 .. total] do
  sub_gens := carat_6_6129_2_filtered_subgroups[k];
  G := Group(sub_gens);
  q_class := CaratQClass(G);
  if (q_class[1] = 6) and (q_class[2] = 2377) then
    Print("\n 6_6129_2 filtered_subgroup nr = ", k, "\n");
    Print("\n", q_class, "\n");
  fi;
od;
