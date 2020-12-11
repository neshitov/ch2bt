// usage gap_saturation_vectors -num_threads n
// compile gcc -fopenmp gap_saturation_vectors.c snf.c -o gap_saturation_vectors

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <stdbool.h>

#include "snf.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

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

  // output all columns of row_t_inverse corresponding to non-invertible pivots
  printf("result\n");
  printf("local sat_vectors;\n");
  printf("sat_vectors:=");

  printf("[ ");
  for(int i = 0; i < nrows; ++i){
    printf("[ ");
    for(int j = 0; j < MIN(nrows, ncols); ++j){
      if (a[j][j] > 1){
        printf("%d, ", row_t_inverse[i][j]);
      }
    }
    printf("], \n");
  }
  printf(" ];\n");
  printf("return sat_vectors;");
}
