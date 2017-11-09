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

  /* // get A and B from user */
  /* printf("Enter the elements of first matrix\n"); */
  /* for (j = 0; j < n; j++) */
  /*   for (i = 0; i < n; i++) */
  /*     scanf("%lf", &A[i + n*j]); // using column-first indexing */
  /* printf("Enter the elements of second matrix\n"); */
  /* for (j = 0; j < n; j++) */
  /*   for (i = 0; i < n; i++) // TODO: optimize allocation to loop nesting */
  /*     scanf("%lf", &B[i + n*j]); // using column-first indexing */

  // make random matrices A and B
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) { // TODO: optimize allocation to loop nesting
      /* A[i + n*j] = random(); */
      /* B[i + n*j] = random(); */
      A[i + n*j] = i + n*j;
      B[i + n*j] = i + n*j;
    }
  }

  // perform multiplication C = AB
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      for (k = 0; k < n; k++) {
          sum = sum + A[i + n*k]*B[k + n*j];
        }

        C[i + n*j] = sum;
        sum = 0;
      }
  }

  // display A and B
  printf("Matrix A:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
        printf("%lf\t", A[i + n*j]);
      printf("\n");
  }
  printf("Matrix B:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
        printf("%lf\t", B[i + n*j]);
      printf("\n");
  }

  printf("Product C of AB:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
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
