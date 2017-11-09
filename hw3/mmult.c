// Author: Eric Kalosa-Kenyon
//
// Sources:
//  1. http://www.programmingsimplified.com/c-program-multiply-matrices
//  2. https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
//  3. http://web.cs.ucdavis.edu/~fgygi/ecs230/homework/hw3/dotblas.c
//  4. http://web.cs.ucdavis.edu/~fgygi/ecs230/homework/hw2/timing1.c

#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "sys/time.h"
#include "math.h"

volatile double gtod(void) {
  /* Directly from timing1.c on ~fgygi */
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv,&tz);
  return tv.tv_sec + 1.e-6*tv.tv_usec;
}

int main(int argc, char** argv)
{

  // initialize variables
  int n = atoi(argv[1]); // size of mx from commandline
  int i, j, k;
  double sum = 0.0;

  // initialize timing elements
  clock_t clk;
  double t, t_cpu, t_real;
  clk = clock();
  t = gtod();

  // dynamic memory allocation for the matrices
  double *A = (double*)malloc((n*n)* sizeof(double));
  double *B = (double*)malloc((n*n)* sizeof(double));
  double *C = (double*)malloc((n*n)* sizeof(double));
  if ((A == NULL) || (B == NULL) || (C == NULL)) {
    fprintf(stderr, "malloc failed\n");
    return(-1);
  }

  /* // get size of matrices from user */
  /* printf("Enter the size (n) of A, B, and hence C (all are n by n):\n"); */
  /* scanf("%d", &n); */
  /* // get A and B from user */
  /* printf("Enter the elements of first matrix\n"); */
  /* for (j = 0; j < n; j++) */
  /*   for (i = 0; i < n; i++) */
  /*     scanf("%lf", &A[i + n*j]); // using column-first indexing */
  /* printf("Enter the elements of second matrix\n"); */
  /* for (j = 0; j < n; j++) */
  /*   for (i = 0; i < n; i++) // TODO: optimize allocation to loop nesting */
  /*     scanf("%lf", &B[i + n*j]); // using column-first indexing */

  // make arbitrary matrices A and B
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) { // TODO: optimize allocation to loop nesting
      /* A[i + n*j] = random(); */
      /* B[i + n*j] = random(); */
      A[i + n*j] = 1.0/(i + n*j + 1);
      B[i + n*j] = 2.0/(i + n*j + 1);
    }
  }

  // perform multiplication C = AB
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      for (k = 0; k < n; k++) {
          sum = sum + A[i + n*k]*B[k + n*j];
        }

        C[i + n*j] = sum;
        sum = 0.0;
      }
  }

  printf("matrix multiplication finished\n");

  // display timing results
  long long delta_clock = clock() - clk;
  t_cpu = ( (double) delta_clock ) / CLOCKS_PER_SEC;
  t_real = gtod() - t;

  printf(" clock:  %15lld  \n", delta_clock );
  printf(" t cpu:  %15.6f s\n", t_cpu );
  printf(" t real: %15.6f s\n", t_real );

  return 0;

  /* // display A and B */
  /* printf("Matrix A:\n"); */
  /* for (i = 0; i < n; i++) { */
  /*   for (j = 0; j < n; j++) */
  /*       printf("%lf\t", A[i + n*j]); */
  /*     printf("\n"); */
  /* } */
  /* printf("Matrix B:\n"); */
  /* for (i = 0; i < n; i++) { */
  /*   for (j = 0; j < n; j++) */
  /*       printf("%lf\t", B[i + n*j]); */
  /*     printf("\n"); */
  /* } */
  /* // display product C = AB */
  /* printf("Product C of AB:\n"); */
  /* for (i = 0; i < n; i++) { */
  /*   for (j = 0; j < n; j++) */
  /*     printf("%lf\t", C[i + n*j]); */
  /*   printf("\n"); */
  /* } */

  return 0;
}
