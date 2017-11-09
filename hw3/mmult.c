// Author: Eric Kalosa-Kenyon
//
// Sources:
//  1. http://www.programmingsimplified.com/c-program-multiply-matrices
//  2. https://en.wikipedia.org/wiki/C_dynamic_memory_allocation

#include <stdio.h>
#include <stdlib.h>

int main()
{
  int n, i, j, k;
  double sum = 0;

  // get size of matrices from user
  printf("Enter the size (n) of A, B, and hence C (all are n by n):\n");
  scanf("%d", &n);

  // dynamic memory allocation for the matrices
  double *A = malloc((n*n)* sizeof(double));
  double *B = malloc((n*n)* sizeof(double));
  double *C = malloc((n*n)* sizeof(double));
  if ((A == NULL) || (B == NULL) || (C == NULL)) {
    fprintf(stderr, "malloc failed\n");
    return(-1);
  }

  // get A and B from user
  printf("Enter the elements of first matrix\n");
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      scanf("%lf", &A[i + n*j]); // using column-first indexing

  printf("Enter the elements of second matrix\n");
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      scanf("%lf", &B[i + n*j]); // using column-first indexing

  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      for (k = 0; k < n; k++) {
          sum = sum + A[i + n*k]*B[k + n*j];
        }

        C[i + n*j] = sum;
        sum = 0;
      }
  }

  printf("Product C of AB:\n");

  for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++)
        printf("%lf\t", C[i + n*j]);
      printf("\n");
    }

  // ## CLEAN UP
  free(A); // clean up unused arrays that were dynamically allocated
  free(B);
  free(C);
  printf("Done!");

  return 0;
}
