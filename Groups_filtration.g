Read("CaratReader.g");
Read("FilterTools.g");

filter_groups := function(d, out_path)
  local i, j, gens, filtered_subgroups;
  filtered_subgroups := [];
  for i in [1 .. Carat_nr_Q_classes(d)] do
    for j in [1 .. Carat_nr_Z_classes(d, i)] do
      gens := Carat(d, i, j);
      if KeepGroup(d, gens) then
        Add(filtered_subgroups, [d, i, j]);
      fi;
    od;
  od;
  PrintTo(out_path, "filtered_subgroups_ids := ");
  AppendTo(out_path, filtered_subgroups);
  AppendTo(out_path, ";");
end;

filter_groups(6, "./dim_6/carat_ids.txt");
