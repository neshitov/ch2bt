#################################################################
# List all subgroups of CARAT 6,6129,4 that are
# not decomposable
# not sign permutations
# have order > 4
#################################################################

Read("FilterTools.g");
Read("CaratReader.g");

gens := Carat(6, 6129, 4);;
hh := ConjugacyClassesSubgroups(Group(gens));;

output_path :="./dim_6_filtered_subgroups/carat_6_6129_4_filtered_subgroups.txt";
PrintTo(output_path, "carat_6_6129_4_filtered_subgroups:=[ ");

for h in hh do
  gens_h := SmallGeneratingSet(Representative(h));
  if KeepGroup(6, gens_h) then
    AppendTo(output_path, gens_h);
    AppendTo(output_path, ", ");
  fi;
od;


AppendTo(output_path, " ];");
