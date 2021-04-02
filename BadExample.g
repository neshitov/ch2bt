Read("CoflasqueCover.g");
Read("CaratReader.g");
Read("ComputePhi.g");
Read("FilterTools.g");
Read("tests/test_coflasque.g");
Read("./dim_6_filtered_subgroups/carat_6_6129_4_filtered_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_2_filtered_subgroups.txt");
Read("./RatProbAlgTori/caratnumber.gap");

for i in [1 .. 4] do
  gens := Carat(6, 6129, i);
  Print([6,6129,i], " keep: ", KeepGroup(6, gens), "\n");
od;

for i in [1 .. 10] do
  A := RandomUnimodularMat(6);
  Print("\n random mat:\n");
  Print(A);
  Print("\n\n Carat = ", [6, 2377, 10], "\n");
  gens := Carat(6, 2377, 10);
  Print("Keep: ", KeepGroup(6, gens), "\n");
  pr := CoflasqueCover(gens);
  cr := CoflasqueResolution(pr, gens);
  Print(ComputePhi(gens, cr));

  gens_twisted := List(gens, x-> A*x*Inverse(A));
  Print("\nKeep_twisted: ", KeepGroup(6, gens_twisted), "\n");
  pr := CoflasqueCover(gens_twisted);
  cr := CoflasqueResolution(pr, gens_twisted);
  Print(ComputePhi(gens_twisted, cr));

  #Print("\nCarat number=", CaratQClass(Group(gens)), "\n");


od;
