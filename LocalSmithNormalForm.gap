PKG_DIR := "/home/alexander/ch2bt";

PrintMatrix := function(file, mat)
  local nrows, ncols, power, i, j;

  #if IsSparseMatrix(mat) then
  #  power := PrimePowersInt(Characteristic(mat!.ring))[2];
  #  nrows := Nrows(mat);
  #  ncols := Ncols(mat);
  #  PrintTo(file, nrows, " ", ncols, " ", power, "\n");
  #  for i in [1 .. nrows] do
  #    for j in [1 .. ncols] do
  #      AppendTo(file, Int(GetEntry(mat, i, j)));
  #      AppendTo(file, " ");
  #    od;
  #    AppendTo(file, "\n");
  #  od;

  #else
  power := PrimePowersInt(Characteristic(Ring(mat[1][1])))[2];
  nrows := DimensionsMat(mat)[1];
  ncols := DimensionsMat(mat)[2];
  PrintTo(file, nrows, " ", ncols, " ", power, "\n");
  for i in [1 .. nrows] do
    for j in [1 .. ncols] do
      AppendTo(file, Int(mat[i][j]));
      AppendTo(file, " ");
    od;
    AppendTo(file, "\n");
  od;
  #fi;
end;


ReadCResult := function(file)
  local tr, power;

  tr := ReadAsFunction(file)();
  power := tr.power;
  tr.SNF :=  tr.SNF * One(Integers mod 2^power);
  tr.row_t :=  tr.row_t * One(Integers mod 2^power);
  tr.row_t_inverse :=  tr.row_t_inverse * One(Integers mod 2^power);
  tr.col_t :=  tr.col_t * One(Integers mod 2^power);
  return tr;

end;

SNFTransform := function(mat)
  local tmp_dir, tmp_input_file, tmp_output_file, exec_file, num_threads, cmd;

  tmp_dir := DirectoryTemporary();
  tmp_input_file := Filename(tmp_dir, "c_input");
  tmp_output_file := Filename(tmp_dir, "c_input_out");
  num_threads := ValueOption("num_threads");
  if num_threads = fail then
    num_threads := 4;
  fi;
  cmd := StringFormatted("src/snf -num_threads {}", num_threads);
  exec_file := Filename(Directory(PKG_DIR), cmd);

  PrintMatrix(tmp_input_file, mat);
  Exec(exec_file, tmp_input_file);
  return ReadCResult(tmp_output_file);

end;
