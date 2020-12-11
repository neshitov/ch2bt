// usage gap_snf -num_threads n input_file or usage snf input_file
// compile gcc -fopenmp gap_snf.c snf.c -o gap_snf

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <stdbool.h>

#include "snf.h"

int NUM_THREADS = 4;



int main(int argc, char** argv){

  int** x;
  int** row_t;
  int** row_t_inverse;
  int** col_t;
  //FILE * fp;
  char* output_name;
  char* input_name;
  int nrows;
  int ncols;
  int power;
  int num_threads;

  if (strcmp(argv[1], "-num_threads")==0){
    num_threads = atoi(argv[2]);
  } else{
    num_threads = NUM_THREADS;
  }

  int** a = read_matrix(&nrows, &ncols, &power);

  row_t = get_identity_matrix(nrows);
  row_t_inverse = get_identity_matrix(nrows);
  col_t = get_identity_matrix(ncols);

  int rank = snf_transform(a, row_t, row_t_inverse, col_t, nrows, ncols, power, num_threads);

  printf("result\n");
  printf("local tr;\n");
  printf("tr:=rec(\n");
  printf("SNF:=");
  write_matrix(a, nrows, ncols);
  printf(",\n");

  printf("rank:= %d,\n", rank);
  printf("power:= %d,\n", power);

  printf("row_t:=");
  write_matrix(row_t, nrows, nrows);
  printf(",\n");

  printf("row_t_inverse:=");
  write_matrix(row_t_inverse, nrows, nrows);
  printf(",\n");

  printf("col_t:=");
  write_matrix(col_t, ncols, ncols);
  printf(");\n return tr;");
}
