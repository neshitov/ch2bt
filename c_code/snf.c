#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <stdbool.h>
#include <omp.h>
//#include <iostream>
//using namespace std;

void switch_rows(int** mat, int i, int j){
  int* temp;
  if (i!=j){
    temp = mat[i];
    mat[i] = mat[j];
    mat[j] = temp;
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

int reduce(int a, int power){
  // a mod order when order is power of 2
  return a & ((1 << power) - 1);
}

int divide(int a, int b, int power){
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
  for (int k = 0; k < ncols; ++k){
    mat[i][k] = reduce(mat[i][k] + c * mat[j][k], power);
  }
}

void add_col_mult(int** mat, int i, int j, int c, int nrows, int power){
  for (int k = 0; k < nrows; ++k){
    mat[k][i] = reduce(mat[k][i] + c * mat[k][j], power);
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
      for (int j = i; j < ncols; j++){
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
      switch_columns(mat, ncols, j_t, t);
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
    for (int i = t + 1; i < nrows; ++i){
      if (mat[i][t] != 0) {
        coef = - divide(mat[i][t], mat[t][t], power);
        add_row_mult(mat, i, t, coef, ncols, power);
        add_row_mult(row_t, i, t, coef, nrows, power);
        add_col_mult(row_t_inverse, i, t, - coef, nrows, power);
      }
    }

    // get rid of all entries in row t
    for (int j = t + 1; j < nrows; ++j){
      if (mat[t][j] != 0) {
        coef = - divide(mat[t][j], mat[t][t], power);
        add_col_mult(mat, j, t, coef, nrows, power);
        add_col_mult(col_t, j, t, coef, ncols, power);
      }
    }

    // Make pivot power of p
    pivot_p_val = two_val(mat[t][t], power);
    int unit = divide((1 << pivot_p_val), mat[t][t], power);
    mat[t][t] = (1 << pivot_p_val);
    multiply_row_by_unit(row_t, t, unit, nrows, power);
    multiply_col_by_unit(row_t_inverse, t, divide(1, unit, power), nrows, power);

  rank_count++;
  }
  return rank_count - 1;
}


int main(){
  int** x;
  int** row_t;
  int** row_t_inverse;
  int** col_t;

  x = malloc(2*sizeof(int));
  row_t = malloc(2*sizeof(int));
  row_t_inverse = malloc(2*sizeof(int));
  col_t = malloc(2*sizeof(int));

  for(int i = 0; i < 2; i++){
    x[i] = malloc(2*sizeof(int));
    row_t[i] = malloc(2*sizeof(int));
    row_t_inverse[i] = malloc(2*sizeof(int));
    col_t[i] = malloc(2*sizeof(int));
  }
  x[0][0] = 2;
  x[0][1] = 2;
  x[1][0] = 3;
  x[1][1] = 4;

  row_t[0][0] = 1;
  row_t[0][1] = 0;
  row_t[1][0] = 0;
  row_t[1][1] = 1;

  row_t_inverse[0][0] = 1;
  row_t_inverse[0][1] = 0;
  row_t_inverse[1][0] = 0;
  row_t_inverse[1][1] = 1;

  col_t[0][0] = 1;
  col_t[0][1] = 0;
  col_t[1][0] = 0;
  col_t[1][1] = 1;

  printf("rank=%d\n", snf_transform(x, row_t, row_t_inverse, col_t, 2, 2, 10));
  //add_row_mult(x, 0, 1, 3, 2, 10);
  //switch_columns(x,2,0,1);

  printf("%d %d\n%d %d\n", x[0][0], x[0][1], x[1][0], x[1][1]);
}
