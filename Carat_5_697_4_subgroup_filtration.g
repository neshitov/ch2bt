#################################################################
# List all subgroups of CARAT 6,6129,4 that are
# not decomposable
# not sign permutations
# have order > 4
#################################################################

Read("FilterTools.g");
Read("CaratReader.g");

gens := Carat(5, 697, 4);;
hh := ConjugacyClassesSubgroups(Group(gens));;

output_path :="./dim_5_filtered_subgroups/carat_5_697_4_filtered_subgroups.txt";
PrintTo(output_path, "carat_5_697_4_filtered_subgroups:=[ ");

for h in hh do
  gens_h := SmallGeneratingSet(Representative(h));
  if KeepGroup(5, gens_h) then
    AppendTo(output_path, gens_h);
    AppendTo(output_path, ", ");
  fi;
od;


AppendTo(output_path, " ];");
