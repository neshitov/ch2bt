
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char** argv){
  int nrows, ncols, power, num_threads;
  int** a;

  if (strcmp(argv[1], "-num_threads")==0){
    num_threads = atoi(argv[2]);
    //input_name = argv[3];
  }

  scanf("%d %d %d", &nrows, &ncols, &power);
  //printf("nrows %d ncols %d\n", nrows, ncols);

  a = malloc(nrows * sizeof(int*));
  for(int i = 0; i < nrows; ++i){
    a[i] = malloc(ncols * sizeof(int));
  }

  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      scanf("%d", &a[i][j]);
      printf("test progress %d / %d rows\n", i, j);
      sleep(1);
      //printf("\rtest progress %d / %d rows", i, j);
      //fflush(stdout);

      //printf("a[%d][%d]=%d\n", i, j, a[i][j]);
    }
  }

  printf("result\n");
  //printf("c_output:=[");
  //for (int i = 0; i < ncols; ++i){
  //  printf("[%d],", a[i][0]);
  //}
  //printf("];\n");
  printf("local huy; huy:=%d; return huy;", num_threads);
}
