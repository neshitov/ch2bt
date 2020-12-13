# Usage tr:=SNFTransform(A : num_threads:=2);;

#PKG_DIR := "/home/alexander/ch2bt";
PKG_DIR := DirectoryCurrent();


SNFTransform := function(mat, power)
  local tmp_dir, tmp_input_file, tmp_output_file,
        exec_file, num_threads, cmd, nrows, ncols, s, path, a_str, gg, tr;
  nrows:= DimensionsMat(mat)[1];
  ncols:= DimensionsMat(mat)[2];
  num_threads := ValueOption("num_threads");
  if num_threads = fail then
    num_threads := 4;
  fi;
  cmd := StringFormatted("src/gap_snf");
  exec_file := Filename(Directory(PKG_DIR), cmd);
  path := DirectoryCurrent();;

  s := InputOutputLocalProcess(path, exec_file, ["-num_threads", String(num_threads)]);
  a_str := String(mat);;
  RemoveCharacters(a_str, "[,]");;

  WriteLine(s, Concatenation(StringFormatted("{} {} {}", nrows, ncols, power), a_str));

  gg := ReadLine(s);

  while gg<>"result\n" do
    RemoveCharacters(gg, "\n");
    Print(Concatenation("\r", gg));;
    gg := ReadLine(s);

  od;
  tr:=ReadAsFunction(s)();

  return tr;

end;

SaturationVectors := function(mat, power)
  local nrows, ncols, num_threads, vectors, cmd, exec_file, path, s, a_str, gg;

  nrows:= DimensionsMat(mat)[1];
  ncols:= DimensionsMat(mat)[2];
  num_threads := ValueOption("num_threads");
  if num_threads = fail then
    num_threads := 4;
  fi;
  cmd := StringFormatted("src/gap_saturation_vectors");
  exec_file := Filename(Directory(PKG_DIR), cmd);
  path := DirectoryCurrent();;

  s := InputOutputLocalProcess(path, exec_file, ["-num_threads", String(num_threads)]);
  a_str := String(mat);;
  RemoveCharacters(a_str, "[,]");;

  WriteLine(s, Concatenation(StringFormatted("{} {} {}", nrows, ncols, power), a_str));

  gg := ReadLine(s);

  while gg<>"result\n" do
    RemoveCharacters(gg, "\n");
    Print(Concatenation("\r", gg));;
    gg := ReadLine(s);

  od;
  vectors := ReadAsFunction(s)();
  CloseStream(s);
  return vectors;

end;
