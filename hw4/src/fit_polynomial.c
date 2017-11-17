// Author: Eric Kalosa-Kenyon
//
// This program performs least squares polynomial fitting using BLAS and LAPACK
//
// Takes:
//      d - command line input, degree of polynomial to fit
//      data.dat - x and y values to fit polynomial to
// Returns:
//      StdIO - fit and intermediate results
//
// Sources:
// 1. https://www.cs.bu.edu/teaching/c/file-io/intro/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "Accelerate/Accelerate.h"
// #include "cblas.h"

// BEGIN MAIN
int main(int argc, char** argv)
{

	// initialize constants and variables
  const char data_fn[] = "../data/data.dat";
  int d = atoi(argv[1]); // degree of fit polynomial from commandline
	float x, y; // for reading x and y from data_fn
	int n = 0; // size of the input data (n by 2 i.e. n xs and n ys)

	// read data into memory from disk
	FILE *ifp; // setup variables
	char *mode = "r";
	ifp = fopen(data_fn, mode);
	if (ifp == NULL) {
	  fprintf(stderr, "Cannot open file\n"); // error if it can't open
	  exit(1);
	}

	// print input contents back out to stdIO and count lines
	while (fscanf(ifp, "%f %f", &x, &y) != EOF) { // for each line in the file
		n++; // count the number of lines
	}
	fclose(ifp);
	ifp = fopen(data_fn, mode);

	// allocate memory for the input data
  double *X = (double*) malloc(n*d* sizeof(double));
	// NOTE: X is col-major indexed i.e. X[i + n*j] = X_(i,j)
  double *Y = (double*) malloc(n* sizeof(double));
  if ((X == NULL) || (Y == NULL)) {
    fprintf(stderr, "malloc failed\n");
    return(1);
  }

	// read input data into the X and Y arrays
	int i = 0, j = 0; // loop counter
	while (fscanf(ifp, "%f %f", &x, &y) != EOF) { // for each line in the file
			for(j = 0; j < d+1; j++){
					X[i + j*n] = pow(x, j); // each col of X is ((x_i)^j)_{i=1..n}, j<=d
			}
			Y[i] = y;
			/* printf("x=%f,\ty=%f\n", x, y); */
			i++;
	}
	fclose(ifp);

	// print the X and Y matrices
	printf("Y\t"); // BEGIN formatting printing header
	for(j = 0; j < d; j++){
			printf("X^%d\t", j);
	}
	printf("\n"); // END formatting printing header

	for(i = 0; i < n; i++){
			printf("%.1f\t", Y[i]);
			for(j = 0; j < d; j++){
				printf("%.2f\t", X[i + j*n]);
			}
			printf("\n");
	}


	// TODO
	// compute A=X^T X

  return 0;
}
// END MAIN
