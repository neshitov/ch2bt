
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
#Read("./dim_6_filtered_subgroups/carat_6_6129_2_filtered_subgroups.txt");


Read("./RatProbAlgTori/caratnumber.gap");
out_path := "./dim_6_filtered_subgroups/carat_6_6129_2_fp_subgroups.txt";

#gens := Carat(6,2377,10);
gens := Carat(6,6129,2);
save_fp_subgroups := function(gens, out_path)
  local G, Gfp, iso, inv_iso, subgroups_fp, subgroup_from_fp;
  G := Group(gens);
  iso := IsomorphismFpGroupByGenerators(G, gens);
  Gfp := Image(iso);
  subgroups_fp := List(ConjugacyClassesSubgroups(Gfp), Representative);
  subgroups_from_fp := List(subgroups_fp, grp -> List(GeneratorsOfGroup(grp),
                                                      x -> PreImagesRepresentative(iso, x)));
  PrintTo(out_path, "carat_6_6129_2_fp_subgroups:=");
  AppendTo(out_path, subgroups_from_fp);
  AppendTo(out_path, ";");
end;

#save_fp_subgroups(gens, out_path);
#Read(out_path);
#Print(carat_6_6129_4_fp_subgroups);
