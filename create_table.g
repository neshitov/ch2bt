Read("./dim_6/result.txt");
Read("CaratReader.g");

table_path := "./dim_6/result_table.tex";

for i in [1 .. Length(result)] do
	if result[i].result.Phi_rank <> 0 then
		carat_id:=result[i].carat_id;
		gens:=Carat(carat_id[1], carat_id[2], carat_id[3]);
		G:=Group(gens);
		fpg:=Image(IsomorphismFpGroupByGenerators(G, gens));
		AppendTo(table_path, result[i].carat_id);
		AppendTo(table_path, " & ");
		AppendTo(table_path, Order(G));
		AppendTo(table_path, " & ");
		tietze_words := List(RelatorsOfFpGroup(fpg), r->TietzeWordAbstractWord(r));
		AppendTo(table_path, tietze_words);
		AppendTo(table_path, " & ");
		AppendTo(table_path, result[i].result.Phi_rank);
		AppendTo(table_path, "\\\\ \n");
	fi;
od;
