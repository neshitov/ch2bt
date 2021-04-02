Read("CoflasqueCover.g");
Read("CaratReader.g");
Read("ComputePhi.g");
Read("tests/test_coflasque.g");
Read("./dim_6_filtered_subgroups/carat_6_6129_4_filtered_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_2_filtered_subgroups.txt");
Read("./RatProbAlgTori/caratnumber.gap");

gens := Carat(6, 6129, 2);
pr := CoflasqueCover(gens);
cr := CoflasqueResolution(pr, gens);


suspicious := [409];

Print("\n6 6129 2:\n");

for i in suspicious do
  Print("\n", i, "\n");
  sub_gens := carat_6_6129_2_filtered_subgroups[i];
  Print("\nCarat number=", CaratZClass(Group(sub_gens)), "\n");
  big_cr := CoflasqueResolution(pr, sub_gens);
  CheckResolution(sub_gens, big_cr);
  Print("\nPhi for big resolution:\n");
  Print(ComputePhi(sub_gens, big_cr));

  small_pr := CoflasqueCover(sub_gens);
  small_cr := CoflasqueResolution(small_pr, sub_gens);
  CheckResolution(sub_gens, small_cr);
  Print("\nPhi for small resolution:\n");
  Print(ComputePhi(sub_gens, small_cr));

od;


gens := Carat(6, 6129, 4);
pr := CoflasqueCover(gens);
cr := CoflasqueResolution(pr, gens);


suspicious := [46, 140, 182, 183, 299, 413];
#suspicious := [46];

Print("\n6 6129 4:\n");

for i in suspicious do
  Print("\n", i, "\n");
  sub_gens := carat_6_6129_4_filtered_subgroups[i];
  Print("\nCarat number=", CaratQClass(Group(sub_gens)), "\n");
  big_cr := CoflasqueResolution(pr, sub_gens);
  CheckResolution(sub_gens, big_cr);
  Print("\nPhi for big resolution:\n");
  Print(ComputePhi(sub_gens, big_cr));

  small_pr := CoflasqueCover(sub_gens);
  small_cr := CoflasqueResolution(small_pr, sub_gens);
  CheckResolution(sub_gens, small_cr);
  Print("\nPhi for small resolution:\n");
  Print(ComputePhi(sub_gens, small_cr));

od;
