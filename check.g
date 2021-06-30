#######################################################################
# Code to compute Phi(G,M) for all filtered 2-subgroups G of GL(6,Z)
#######################################################################

Read("ComputePhi.g");
Read("CaratReader.g");
Read("CoflasqueCover.g");
Read("config.g");

check_dimension := function(carat_ids_path, result_path)
  local i, id, gens, pr, cr, result, res;

  Read(carat_ids_path);
  PrintTo(result_path, "result := [");

  for i in [1 .. Length(filtered_subgroups_ids)] do
  #for i in [1 .. 2] do
    id := filtered_subgroups_ids[i];

    #if id <> [3, 14, 2] then
    if true then
      gens := Carat(id[1], id[2], id[3]);
      Reset(GlobalMersenneTwister);
      Reset(GlobalRandomSource);

      pr := CoflasqueCover(gens);
      cr := CoflasqueResolution(pr, gens);

      result := rec();
      result.carat_id := id;
      result.coflasque_resolution_dim := DimensionsMat(cr.actionC[1])[1];

      if result.coflasque_resolution_dim > 1 then
        res := ComputePhi(gens, cr: num_threads:=n_threads);
      else
        res := rec(HM2_rank := -1, Phi_rank := 0,
                   im_HL_rank := 0, im_HP2_rank :=-1);
      fi;
      result.result := res;
      AppendTo(result_path, result, ",");
    fi;
  od;
  AppendTo(result_path, "];");

end;

#check_dimension("./dim_3/carat_ids.txt", "./dim_3/result.txt");
#check_dimension("./dim_4/carat_ids.txt", "./dim_4/result.txt");
check_dimension("./dim_5/carat_ids.txt", "./dim_5/result.txt");
#check_dimension("./dim_6/carat_ids.txt", "./dim_6/result.txt");
