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
	float x, y; // for reading x and y from data_fn

	printf("%s\n", data_fn);

	// read data into memory from disk
	FILE *ifp;
	char *mode = "r";

	ifp = fopen(data_fn, mode);

	if (ifp == NULL) {
	  fprintf(stderr, "Cannot open file\n");
	  exit(1);
	}

	printf("X,\tY\n");
	while (fscanf(ifp, "%f %f", &x, &y) != EOF) {
	  printf("%f, %f\n", x, y);
	}

  return 0;
}
// END MAIN
