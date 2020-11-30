#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
#include <stdbool.h>
#include <omp.h>
//#include <iostream>
//using namespace std;

int** allocate_memory(int nrows, int ncols){
  int** a;
  a = malloc(nrows * sizeof(int*));
  for(int i = 0; i < nrows; ++i){
    a[i] = malloc(ncols * sizeof(int));
  }
  return a;
}

int reduce(int a, int power){
  // a mod order when order is power of 2
  return a & ((1 << power) - 1);
}

void print_matrix(int** a, int nrows, int ncols){
  for(int i = 0; i < nrows; ++i){
      for(int j = 0; j < ncols; ++j){
        printf("%d ", a[i][j]);
      }
      printf("\n");
  }
}

void switch_rows(int** mat, int i, int j){
  int* temp;
  if (i!=j){
    temp = mat[i];
    mat[i] = mat[j];
    mat[j] = temp;
  }
}

void permute_rows(int** mat, int* p_original, int perm_len){
  int next;
  int temp;
  int p[perm_len];
  for (int i = 0; i < perm_len; ++i){
    p[i] = p_original[i];
  }

  for (int i = 0; i < perm_len; ++i){
    next = i;
    while (p[next] >= 0){
      switch_rows(mat, i, p[next]);
      temp = p[next];
      p[next] -= perm_len;
      next = temp;
    }
  }
}

void switch_columns(int** mat, int nrows, int i, int j){
  int temp;
  if (i!=j){
    for (int k = 0; k < nrows; k++){
      temp = mat[k][i];
      mat[k][i] = mat[k][j];
      mat[k][j] = temp;
    }
  }
}

void permute_columns(int** mat, int* p_original, int nrows, int perm_len){
  int next;
  int temp;
  int p[perm_len];

  for (int i = 0; i < perm_len; ++i){
    p[i] = p_original[i];
  }

  for (int i = 0; i < perm_len; ++i){
    next = i;
    while (p[next] >= 0){
      switch_columns(mat, nrows, i, p[next]);
      temp = p[next];
      p[next] -= perm_len;
      next = temp;
    }
  }
}


int two_val(int n, int power){
  if (n == 0){
    return power;
  } else {
    return __builtin_ctz(n);
  }
}

int inverse(int a, int power){
  // return 1/a mod order where order is power of 2
  int order = 1 << power;
  int s = two_val(a - 1, power);
  int t = a / (1 << s);
  int u = 2 - a;
  int amone = a - 1;
  for (int i = 1; i < power / s; i <<= 1){
    amone *= amone;
    amone &= (order - 1);
    u *= (amone + 1);
    u &= (order - 1);
  }
  return u;
}


int divide(int a, int b, int power){
  // return  a / b
  int val_a = two_val(a, power);
  int val_b = two_val(b, power);
  int b_odd = b / (1 << val_b);
  int a_odd = a / (1 << val_a);
  return reduce(a_odd * inverse(b_odd, power) * (1 << (val_a - val_b)), power);
}

void multiply_row_by_unit(int** mat, int i, int unit, int ncols, int power){
  for (int k = 0; k < ncols; ++k){
    mat[i][k] = reduce(mat[i][k] * unit, power);
  }
}

void multiply_col_by_unit(int** mat, int i, int unit, int nrows, int power){
  for (int k = 0; k < nrows; ++k){
    mat[k][i] = reduce(mat[k][i] * unit, power);
  }
}

void add_row_mult(int** mat, int i, int j, int c, int ncols, int power){
  #pragma omp parallel num_threads(4)
  {
  #pragma omp for

  for (int k = 0; k < ncols; ++k){
    mat[i][k] = reduce(mat[i][k] + c * mat[j][k], power);
  }
}
}

void add_col_mult(int** mat, int i, int j, int c, int nrows, int power){
  #pragma omp parallel num_threads(4)
  {
  #pragma omp for

  for (int k = 0; k < nrows; ++k){
    mat[k][i] = reduce(mat[k][i] + c * mat[k][j], power);
  }
}
}

void find_smallest_row_col_entries(
  int** mat,
  int nrows,
  int ncols,
  int t,
  int* min_p_val_in_column,
  int* min_p_val_in_column_index,
  int* min_p_val_in_row,
  int* min_p_val_in_row_index,
  int power
){

  // Finds entries with minimal valuations in row and column containing (t, j_t)
  int entry_p_val;
  int pivot_p_val = two_val(mat[t][t], power);
  *min_p_val_in_column = pivot_p_val;
  *min_p_val_in_column_index = t;
  *min_p_val_in_row = pivot_p_val;
  *min_p_val_in_row_index = t;

  // find smallest 2 valuation in column
  for (int i = t + 1; i < nrows; i++){
    entry_p_val = two_val(mat[i][t], power);
    if (entry_p_val < *min_p_val_in_column){
      *min_p_val_in_column = entry_p_val;
      *min_p_val_in_column_index = i;
    }
  }

  // find smallest 2 valuation in row
  for (int j = t + 1; j < ncols; j++){
    entry_p_val = two_val(mat[t][j], power);
    if (entry_p_val < *min_p_val_in_row){
      *min_p_val_in_row = entry_p_val;
      *min_p_val_in_row_index = j;
    }
  }
}


int snf_transform(
  int** mat, int** row_t,
  int** row_t_inverse, int** col_t,
  int nrows, int ncols,
  int power
){


  int rank_count = 1;
  int min_p_val_in_column;
  int min_p_val_in_column_index;
  int min_p_val_in_row;
  int min_p_val_in_row_index;
  int pivot_p_val;
  int coef;


  for (int t = 0; t < nrows; ++t){

    // Choose j_t such that A[>=t,j_t] is the leftmost nonzero column in A[>=t,:]
    int j_t = ncols + 1;
    int min_v = power;
    int i_t = 0;

    for (int i = t; i < nrows; i++){
      for (int j = t; j < ncols; j++){
        if (mat[i][j] != 0) {
          if (j == t){
            i_t = i;
            j_t = j;
            goto end;
          }
          if (j < j_t){
            j_t = j;
            i_t = i;
          }
          break;
        }
      }
    }
    end:

    // if A[>=t,:] is zero, exit loop
    if (j_t == ncols + 1){
      break;
    }

    // Move pivot to position to (t, t)
    if (i_t != t){
      switch_rows(mat, i_t, t);
      switch_rows(row_t, i_t, t);
      switch_columns(row_t_inverse, nrows, i_t, t);
    }
    if (j_t != t){
      switch_columns(mat, nrows, j_t, t);
      switch_columns(col_t, ncols, j_t, t);
    }


    // Make pivot divide all elements in its row and column
    bool keep_going = true;
    while (keep_going){
      pivot_p_val = two_val(mat[t][t], power);
      find_smallest_row_col_entries(mat, nrows, ncols, t,
                                    &min_p_val_in_column, &min_p_val_in_column_index,
                                    &min_p_val_in_row, &min_p_val_in_row_index,
                                    power);

      if ((pivot_p_val == min_p_val_in_column)
          && (pivot_p_val == min_p_val_in_row)){
            keep_going = false;
      }
      else if (min_p_val_in_column <= min_p_val_in_row){
        add_row_mult(mat, t, min_p_val_in_column_index, 1, ncols, power);
        add_row_mult(row_t, t, min_p_val_in_column_index, 1, nrows, power);
        add_col_mult(row_t_inverse, min_p_val_in_column_index, t, -1, nrows, power);
      }
      else {
        add_col_mult(mat, t, min_p_val_in_row_index, 1, nrows, power);
        add_col_mult(col_t, t, min_p_val_in_row_index, 1, ncols, power);
      }
    }


    // get rid of all entries in column t

    //#pragma omp parallel num_threads(4)
    //{
    //#pragma omp for
      for (int i = t + 1; i < nrows; ++i){
        if (mat[i][t] != 0) {
          coef = - divide(mat[i][t], mat[t][t], power);
          //printf("R%d <- R%d + %d * R%d\n", i, i, coef, t);
          add_row_mult(mat, i, t, coef, ncols, power);
          add_row_mult(row_t, i, t, coef, nrows, power);
          add_col_mult(row_t_inverse, t, i, - coef, nrows, power);
        }
      }
    //}
    //printf("\nsnf:\n");
    //print_matrix(mat, nrows, ncols);
    //printf("\nrow_t:\n");
    //print_matrix(row_t, nrows, nrows);
    //printf("\ncol_t:\n");
    //print_matrix(col_t, ncols, ncols);


    // get rid of all entries in row t

    //#pragma omp parallel num_threads(4)
    //{
    //#pragma omp for
      for (int j = t + 1; j < ncols; ++j){
        if (mat[t][j] != 0) {
          coef = - divide(mat[t][j], mat[t][t], power);
          add_col_mult(mat, j, t, coef, nrows, power);
          add_col_mult(col_t, j, t, coef, ncols, power);
        }
      }
    //}

    // Make pivot power of p
    pivot_p_val = two_val(mat[t][t], power);
    int unit = divide((1 << pivot_p_val), mat[t][t], power);
    mat[t][t] = (1 << pivot_p_val);
    multiply_row_by_unit(row_t, t, unit, nrows, power);
    multiply_col_by_unit(row_t_inverse, t, divide(1, unit, power), nrows, power);

  rank_count++;
  printf("c_code: %d done\n", t);
  //printf("\nsnf:\n");
  //print_matrix(mat, nrows, ncols);
  //printf("\nrow_t:\n");
  //print_matrix(row_t, nrows, nrows);
  //printf("\ncol_t:\n");
  //print_matrix(col_t, ncols, ncols);
  }
  rank_count--;

  // Rearrange pivots in divisible ascending order
  int index[rank_count];
  for (int i = 0; i < rank_count; ++i){index[i] = i;}

  int compare (const void * a, const void * b)
  {
    return ( mat[*(int*)a][*(int*)a] - mat[*(int*)b][*(int*)b] );
  }
  qsort(index, rank_count, sizeof(int), compare);
  permute_rows(mat, index, rank_count);
  permute_rows(row_t, index, rank_count);
  permute_columns(row_t_inverse, index, nrows, rank_count);

  permute_columns(mat, index, nrows, rank_count);
  permute_columns(col_t, index, ncols, rank_count);

  return rank_count;
}


int** read_matrix(int* nrows, int* ncols, int* power, FILE * fp){
  int** a;
  fscanf(fp, "%d %d %d", nrows, ncols, power);
  //printf("nrows %d ncols %d\n", *nrows, *ncols);

  a = malloc(*nrows * sizeof(int*));
  for(int i = 0; i < *nrows; ++i){
    a[i] = malloc(*ncols * sizeof(int));
  }

  for(int i = 0; i < *nrows; i++){
    for(int j = 0; j < *ncols; j++){
      fscanf(fp, "%d", &a[i][j]);
      //printf("%d %d %d\n", i, j, a[i][j]);
    }
  }
  return a;
}

void write_matrix(int** x, int nrows, int ncols, FILE * fp){
  fprintf(fp, "[ ");
  for(size_t i = 0; i < nrows; ++i){
    fprintf(fp, "[ ");
    for(size_t j = 0; j < ncols; ++j){
            fprintf(fp, "%d, ", x[i][j]);
    }
    fprintf(fp, "], \n");
  }
  fprintf(fp, " ]");
}


int** get_identity_matrix(int nrows){
  int ** a;
  a = malloc(nrows * sizeof(int*));
  for(size_t i = 0; i < nrows; ++i){
    a[i] = malloc(nrows *sizeof(int));
    for(size_t j = 0; j < nrows; ++j){
      if (j==i){
        a[i][j] = 1;
      } else{
        a[i][j] = 0;
      }
    }
  }
  return a;
}


int main(int argc, char** argv){

  int** x;
  int** row_t;
  int** row_t_inverse;
  int** col_t;
  FILE * fp;
  char* output_name;
  int nrows;
  int ncols;
  int power;


  /*
  fp = fopen("/home/alexander/ch2bt/c_code/debug.txt", "r");
  int** a = read_matrix(&nrows, &ncols, &power, fp);
  fclose(fp);

  row_t = get_identity_matrix(nrows);
  row_t_inverse = get_identity_matrix(nrows);
  col_t = get_identity_matrix(ncols);

  print_matrix(a, nrows, ncols);
  */


  output_name = malloc(strlen(argv[1]) + 4);
  strcpy(output_name, argv[1]);
  strcat(output_name, "_out");

  fp = fopen(argv[1], "r");
  int** a = read_matrix(&nrows, &ncols, &power, fp);
  fclose(fp);

  row_t = get_identity_matrix(nrows);
  row_t_inverse = get_identity_matrix(nrows);
  col_t = get_identity_matrix(ncols);

  int rank = snf_transform(a, row_t, row_t_inverse, col_t, nrows, ncols, power);

  /*<----------------
  printf("nrows=%d ncols=%d power=%d\n", nrows, ncols, power);
  print_matrix(a, nrows, ncols);


  printf("rank=%d\n", rank);
  //add_row_mult(x, 0, 1, 3, 2, 10);
  printf("snf:\n");
  print_matrix(a, nrows, ncols);

  printf("row_t:\n");
  print_matrix(row_t, nrows, nrows);

  printf("row_t_inverse:\n");
  print_matrix(row_t_inverse, nrows, nrows);
  <-------------------*/

  fp = fopen(output_name, "w");
  fprintf(fp, "local tr;\n");
  fprintf(fp, "tr:=rec(\n");
  fprintf(fp, "SNF:=");
  write_matrix(a, nrows, ncols, fp);
  fprintf(fp, ",\n");

  fprintf(fp, "rank:= %d,\n", rank);
  fprintf(fp, "power:= %d,\n", power);

  fprintf(fp, "row_t:=");
  write_matrix(row_t, nrows, nrows, fp);
  fprintf(fp, ",\n");

  fprintf(fp, "row_t_inverse:=");
  write_matrix(row_t_inverse, nrows, nrows, fp);
  fprintf(fp, ",\n");

  fprintf(fp, "col_t:=");
  write_matrix(col_t, ncols, ncols, fp);
  fprintf(fp, ");\n return tr;");
  fclose(fp);
}
