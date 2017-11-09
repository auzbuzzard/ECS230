// Author: Eric Kalosa-Kenyon
//
// Sources:
//  1. http://www.programmingsimplified.com/c-program-multiply-matrices
//  2. https://en.wikipedia.org/wiki/C_dynamic_memory_allocation

#include <stdio.h>
#include <stdlib.h>

int main()
{
  int m, n, p, q, c, d, k, sum = 0; // old author's
  int N, i, j; // mine

  // get size of matrices from user
  printf("Enter the size (n) of A, B, and hence C (all are n by n):\n");
  scanf("%d", &N);

  // dynamic memory allocation for the matrices
  double *A = malloc((N*N)* sizeof(double));
  double *B = malloc((N*N)* sizeof(double));
  double *C = malloc((N*N)* sizeof(double));
  if ((A == NULL) || (B == NULL) || (C == NULL)) {
    fprintf(stderr, "malloc failed\n");
    return(-1);
  }

  // get A from user
  printf("Enter the elements of first matrix\n");
  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
      scanf("%lf", &A[i + n*j]); // using column-first indexing

  /* printf("Enter the number of rows and columns of second matrix\n"); */
  /* scanf("%d%d", &p, &q); */

  /* if (n != p) */
  /*   printf("Matrices with entered orders can't be multiplied with each other.\n"); */
  /* else */
  /* { */
  /*     printf("Enter the elements of second matrix\n"); */

  /*     for (c = 0; c < p; c++) */
  /*       for (d = 0; d < q; d++) */
  /*         scanf("%d", &second[c][d]); */

  /*     for (c = 0; c < m; c++) { */
  /*           for (d = 0; d < q; d++) { */
  /*                   for (k = 0; k < p; k++) { */
  /*                             sum = sum + first[c][k]*second[k][d]; */
  /*                           } */

  /*                   multiply[c][d] = sum; */
  /*                   sum = 0; */
  /*                 } */
  /*         } */

  /*     printf("Product of entered matrices:-\n"); */

  /*     for (c = 0; c < m; c++) { */
  /*           for (d = 0; d < q; d++) */
  /*             printf("%d\t", multiply[c][d]); */
  /*           printf("\n"); */
  /*         } */
  /*   } */

  // ## CLEAN UP
  free(A); // clean up unused arrays that were dynamically allocated
  free(B);
  free(C);
  printf("Done!");

  return 0;
}
