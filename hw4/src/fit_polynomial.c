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

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "Accelerate/Accelerate.h"
// #include "cblas.h"

// BEGIN MAIN
int main(int argc, char** argv)
{

	// initialize constants and variables
  const char data_path[] = "../data/data.dat";

	printf("%s\n", data_path);

  return 0;
}
// END MAIN
