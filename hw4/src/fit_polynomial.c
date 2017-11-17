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
	ifp = fopen(data_fn, mode); // open file
	if (ifp == NULL) {
	  fprintf(stderr, "Cannot open file\n"); // error if it can't open
	  exit(1);
	}

	// print input contents back out to stdIO and count lines
	printf("X,\tY\n"); // output formatting
	while (fscanf(ifp, "%f %f", &x, &y) != EOF) { // for each line in the file
	  printf("%f, %f\n", x, y); // print the line
		n++; // count the number of lines
	}
	printf("\n");

	// allocate memory for the input data
  double *X = (double*) malloc(n* sizeof(double));
  double *Y = (double*) malloc(n* sizeof(double));
  if ((X == NULL) || (Y == NULL)) {
    fprintf(stderr, "malloc failed\n");
    return(1);
  }

	// read input data into the X and Y arrays
	int i = 0; // loop counter
	while (fscanf(ifp, "%f %f", &x, &y) != EOF) { // for each line in the file
			X[i] = x;
			Y[i] = y;
			i++;
	}

	// close the file
	fclose(ifp);

  return 0;
}
// END MAIN
