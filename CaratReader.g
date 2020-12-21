#########################################################
# Function to read CARAT tables
#########################################################

carat_folder := "./carat_tables";
files := [
  Concatenation(carat_folder, "/cryst1.txt"),
  Concatenation(carat_folder, "/cryst2.txt"),
  Concatenation(carat_folder, "/cryst3.txt"),
  Concatenation(carat_folder, "/cryst4.txt"),
  Concatenation(carat_folder, "/cryst5.txt"),
  Concatenation(carat_folder, "/cryst6.txt"),
];

for file in files do
  Read(file);
od;
cryst := [cryst1, cryst2, cryst3, cryst4, cryst5, cryst6];


Carat := function(d, n, i)
  return cryst[d][n][i];
end;