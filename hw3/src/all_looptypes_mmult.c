// Author: Eric Kalosa-Kenyon
//
// This program times matrix multiplication using each loop ordering and reports
// timing results.
//
// Takes:
//      n - command line input, size of matrices to multiply
// Returns:
//      stdio - timing results
//
// Sources:
//  1. http://www.programmingsimplified.com/c-program-multiply-matrices
//  2. https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
//  3. http://web.cs.ucdavis.edu/~fgygi/ecs230/homework/hw3/dotblas.c
//  4. http://web.cs.ucdavis.edu/~fgygi/ecs230/homework/hw2/timing1.c
//  5. http://www.cplusplus.com/forum/windows/192252/

#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "sys/time.h"
#include "math.h"

// BEGIN SUBROUTINES
volatile double gtod(void) {
  /* Get time of day */
  /* Directly from timing1.c on ~fgygi */
  /* http://web.cs.ucdavis.edu/~fgygi/ecs230/homework/hw2/timing1.c */
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv,&tz);
  return tv.tv_sec + 1.e-6*tv.tv_usec;
}

long long readTSC(void)
{
  /* Read the time stamp counter on Intel x86 chips */
  /* Directly from timing2.c on ~fgygi */
  /* http://web.cs.ucdavis.edu/~fgygi/ecs230/homework/hw2/timing2.c */
  union { long long complete; unsigned int part[2]; } ticks;
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
    : "=mr" (ticks.part[0]), "=mr" (ticks.part[1])
    : /* no inputs */
    : "eax", "edx");
  return ticks.complete;
}
// END SUBROUTINES

// BEGIN MAIN
int main(int argc, char** argv)
{

  // initialize matrix size and indexing variables
  int n = atoi(argv[1]); // size of mx from commandline
  int r = atoi(argv[2]); // number of times to perform the multiplication
  int i, j, k, s; // indicies for matrix arrays
  double sum = 0.0;

  // dynamic memory allocation for the matrices
  double *A = (double*)malloc((n*n)* sizeof(double));
  double *B = (double*)malloc((n*n)* sizeof(double));
  double *C = (double*)malloc((n*n)* sizeof(double));
  if ((A == NULL) || (B == NULL) || (C == NULL)) {
    fprintf(stderr, "malloc failed\n");
    return(-1);
  }

  // make arbitrary matrices A and B
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      A[i + n*j] = 1.0/(i + n*j + 1);
      B[i + n*j] = 2.0/(i + n*j + 1);
    }
  }

  // initialize timing elements
  clock_t clk;
  double tod, t_cpu, t_real;
  long start, stop, time_elapsed;
  long long delta_clock;

  // perform AB r times
  printf("clocks, time, loop_order\n");
  for (s = 0; s < r; s++) {

    // start timing
    clk = clock();
    start = readTSC();
    tod = gtod();
    // perform multiplication C = AB - using ordering 1 jik
    for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            sum = sum + A[i + n*k]*B[k + n*j];
          }
          C[i + n*j] = sum;
          sum = 0.0;
        }
    }
    // determine and display timing results
    delta_clock = clock() - clk;
    t_real = gtod() - tod;
    stop = readTSC();
    time_elapsed = stop - start;
    printf("%ld, %f, %d\n", time_elapsed, t_real, 1 );

    // start timing
    clk = clock();
    start = readTSC();
    tod = gtod();
    // perform multiplication C = AB - using ordering 2 jki
    for (j = 0; j < n; j++) {
      for (k = 0; k < n; k++) {
        for (i = 0; i < n; i++) {
            sum = sum + A[i + n*k]*B[k + n*j];
          }
          C[i + n*j] = sum;
          sum = 0.0;
        }
    }
    // determine and display timing results
    delta_clock = clock() - clk;
    t_real = gtod() - tod;
    stop = readTSC();
    time_elapsed = stop - start;
    printf("%ld, %f, %d\n", time_elapsed, t_real, 2 );

    // start timing
    clk = clock();
    start = readTSC();
    tod = gtod();
    // perform multiplication C = AB - using ordering 3 ijk
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++) {
            sum = sum + A[i + n*k]*B[k + n*j];
          }
          C[i + n*j] = sum;
          sum = 0.0;
        }
    }
    // determine and display timing results
    delta_clock = clock() - clk;
    t_real = gtod() - tod;
    stop = readTSC();
    time_elapsed = stop - start;
    printf("%ld, %f, %d\n", time_elapsed, t_real, 3 );

    // start timing
    clk = clock();
    start = readTSC();
    tod = gtod();
    // perform multiplication C = AB - using ordering 4 ikj
    for (i = 0; i < n; i++) {
			for (k = 0; k < n; k++) {
				for (j = 0; j < n; j++) {
            sum = sum + A[i + n*k]*B[k + n*j];
          }
          C[i + n*j] = sum;
          sum = 0.0;
        }
    }
    // determine and display timing results
    delta_clock = clock() - clk;
    t_real = gtod() - tod;
    stop = readTSC();
    time_elapsed = stop - start;
    printf("%ld, %f, %d\n", time_elapsed, t_real, 4 );

    // start timing
    clk = clock();
    start = readTSC();
    tod = gtod();
    // perform multiplication C = AB - using ordering 5 kji
    for (k = 0; k < n; k++) {
			for (j = 0; j < n; j++) {
				for (i = 0; i < n; i++) {
            sum = sum + A[i + n*k]*B[k + n*j];
          }
          C[i + n*j] = sum;
          sum = 0.0;
        }
    }
    // determine and display timing results
    delta_clock = clock() - clk;
    t_real = gtod() - tod;
    stop = readTSC();
    time_elapsed = stop - start;
    printf("%ld, %f, %d\n", time_elapsed, t_real, 5 );

    // start timing
    clk = clock();
    start = readTSC();
    tod = gtod();
    // perform multiplication C = AB - using ordering 6 kij
    for (k = 0; k < n; k++) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
            sum = sum + A[i + n*k]*B[k + n*j];
          }
          C[i + n*j] = sum;
          sum = 0.0;
        }
    }
    // determine and display timing results
    delta_clock = clock() - clk;
    t_real = gtod() - tod;
    stop = readTSC();
    time_elapsed = stop - start;
    printf("%ld, %f, %d\n", time_elapsed, t_real, 6 );

  } // END for loop over R (multiply AB multiple times)

  // Remove the matrices from memory
  free(A);
  free(B);
  free(C);

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
// END MAIN
