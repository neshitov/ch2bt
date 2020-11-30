// gcc -fopenmp -o go go.c
// ./go

#include <stdio.h>
#include <stdlib.h>

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

//DEBUG:
int** multiply(int** left, int left_nrows, int left_ncols,
               int** right, int right_nrows, int right_ncols, int power){
  int tmp;
  int** prod = allocate_memory(left_nrows, right_ncols);
  for (int i = 0; i < left_nrows; ++i){
    for (int j = 0; j < right_ncols; ++j){
      tmp = 0;
      for (int k = 0; k < left_ncols; ++k){
        tmp = reduce(tmp + left[i][k] * right[k][j], power);
      }
      prod[i][j] = reduce(tmp, power);
      //printf("multiply %d\n", prod[i][j]);
    }
  }
  return prod;
}


int divide(int a, int b, int power){
  // return  a / b
  int val_a = two_val(a, power);
  int val_b = two_val(b, power);
  int b_odd = b / (1 << val_b);
  int a_odd = a / (1 << val_a);
  return reduce(a_odd * inverse(b_odd, power) * (1 << (val_a - val_b)), power);
}


int main(int argc, char **argv){
  int power = 9;
  int x = 511;
  int y = -1;
  printf("%d %d\n", reduce(x, power), reduce(y, power));
  printf("%d %d\n", inverse(x, power), inverse(y, power));
  printf("%d", divide(x,y, power));
}
