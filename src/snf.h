//#include <omp.h>

int snf_transform(
  int** mat, int** row_t,
  int** row_t_inverse, int** col_t,
  int nrows, int ncols,
  int power,
  int _num_threads
);

int** read_matrix(
  int* nrows,
  int* ncols,
  int* power
);

void write_matrix(
  int** x,
  int nrows,
  int ncols
);

int** get_identity_matrix(
  int nrows
);
