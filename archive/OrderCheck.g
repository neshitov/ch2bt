Read("CoflasqueCover.g");
Read("CaratReader.g");
Read("ComputePhi.g");
Read("tests/test_coflasque.g");

Read("CaratReader.g");
Read("FilterTools.g");
#Read("./dim_6_filtered_subgroups/carat_6_6129_4_filtered_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_1_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_2_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_3_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_4_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_4_fp_subgroups.txt");
Read("./dim_6_filtered_subgroups/carat_6_6129_2_fp_subgroups.txt");
#Read("./dim_6_filtered_subgroups/carat_6_6129_2_filtered_subgroups.txt");
Read("./carat_tables/caratchpol.txt");
Read("./RatProbAlgTori/caratnumber.gap");

characters := function(G)
  local chars;
  chars := List(ConjugacyClasses(G), x -> Trace(Representative(x)));
  Sort(chars);
  return chars;
end;

write_6_6129_subgroups := function(i, output_path)
  local gens, hh, h, gens_h;
  gens := Carat(6, 6129, i);
  hh := ConjugacyClassesSubgroups(Group(gens));;


  PrintTo(output_path, "carat_6_6129_");
  AppendTo(output_path, i);
  AppendTo(output_path, "_subgroups:=[ ");

  for h in hh do
    gens_h := SmallGeneratingSet(Representative(h));
    AppendTo(output_path, gens_h);
    AppendTo(output_path, ", ");
  od;

  AppendTo(output_path, " ];");
end;

#write_6_6129_subgroups(1, "./dim_6_filtered_subgroups/carat_6_6129_1_subgroups.txt");
#write_6_6129_subgroups(3, "./dim_6_filtered_subgroups/carat_6_6129_3_subgroups.txt");

carat_subgroups := [carat_6_6129_1_subgroups, carat_6_6129_2_subgroups,
                   carat_6_6129_3_subgroups, carat_6_6129_4_subgroups];

look_up := function(nr)
  local gens, i, sub_gens, target_h1, subgrp, grp, pr, cr;
  gens := Carat(6, 2377, 10);
  grp := Group(gens);
  Print([6,2377,10], "\n");
  target_h1 := Filtered(H1(Group(gens)), i-> i<>1);
  Print("H1=", target_h1, "\n");
  Print("H0=", H0(Group(gens)), "\n");
  Print("Hminus1=", Hminus1(Group(gens)), "\n");
  Print(characters(grp),"\n");

  pr := CoflasqueCover(gens);
  cr := CoflasqueResolution(pr, gens);
  Print(ComputePhi(gens, cr));
  Print("\n");

  for i in [2 .. Length(carat_subgroups[nr])] do
    sub_gens := carat_subgroups[nr][i];
    subgrp := Group(sub_gens);
      #if CaratQClass(Group(sub_gens))[2] = 2377 then
    if (Order(subgrp) = 8)
        and (Filtered(H1(subgrp), x -> x<>1) = target_h1)
        and (characters(subgrp) = characters(grp))
        #and CharPolyList(subgrp) = CharPolyList(Group(gens)) then
        then
      Print("(6,6129, ", nr, ") subgroup ", i, "\n");
      Print("Carat Q class = ", CaratQClass(Group(sub_gens)), "\n");
      Print(characters(subgrp),"\n");
      Print("H1=", Filtered(H1(subgrp), i -> i<>1), "\n");
      Print("H0=", H0(subgrp), "\n");
      Print("Hminus1=", Hminus1(subgrp), "\n");

      pr := CoflasqueCover(sub_gens);
      cr := CoflasqueResolution(pr, sub_gens);
      Print(ComputePhi(sub_gens, cr));
      Print("\n\n");

    fi;
    #nr := nr + 1;
  od;

end;

for i in [1 .. 4] do
  look_up(i);
od;
#gens := Carat(6, 6129, 1);
#Print(Order(Group(gens)));
