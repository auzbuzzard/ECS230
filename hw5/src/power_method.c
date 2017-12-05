// Author: Eric Kalosa-Kenyon
//
// This program performs the power method to determine the leading eigenvector
// of a column-stochastic link matrix.
//
// Takes:
//      data.dat - input data (format not explicitly specified here)
// Returns:
//      StdIO - results
//      File - results
//
// Sources:

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "Accelerate/Accelerate.h"
// #include "cblas.h"

// BEGIN PROTOTYPES
void dgemm_(char * transa, char * transb, int * m, int * n, int * k,
    double * alpha, double * A, int * lda,
    double * B, int * ldb, double * beta,
    double *, int * ldc);

void dgemv_(char * trans, int * m, int * n,
    double * alpha, double * A, int * lda,
    double * x, int * incx, double * beta,
    double * y, int * incy);

int dpotrf_(char * uplo, int * n,
    double * A, int * lda, int * info);

void dtrsm_(char * side, char * uplo, char * transa, char * diag,
    int * m, int * n, double * alpha, double * A,
    int * lda, double * B, int * ldb);
// END PROTOTYPES

// BEGIN VARIABLES
// for general configuration
/* const char data_fn[] = "../data/data.dat"; // location of the input data */
const char data_fn[] = "../data/test_A.dat"; // location of the input data
int i, j; // loop counters
int n; // size of graph's vertex set
int firstline; // determines whether reading first line of data.dat or not

// file handling objects
FILE *ifp; // pointer to the file object containing the input data
char *mode = "r"; // opening mode for the file above
FILE * out_raw; // where to write raw output
FILE * gnuplotPipe;

// CBLAS/LAPACK inputs
char TRANSA;
char TRANSB;
char TRANS;
double ALPHA;
double BETA;
int M;
int N;
int K;
int LDA;
int LDB;
int LDC;
int INCB;
int INCY;
int INCP;
char SIDE;
char DIAG;

// END VARIABLES

// BEGIN MAIN
int main(int argc, char** argv)
{
    // read data into memory from disk
    ifp = fopen(data_fn, mode);
    if (ifp == NULL) {
        fprintf(stderr, "Cannot open file\n"); // error if it can't open
        exit(1);
    }

    n = 1;
    firstline = 1;
    double *A = (double*) malloc(n*n* sizeof(double)); // placeholder malloc

    // read file
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * pch;
    i = 0;
    while ((read = getline(&line, &len, ifp)) != -1){
        if(firstline == 1){

            printf("Read first line (len %zu) : ", read);
            printf("%s\n", line);
            n = atoi(line); // n_obs is first line of data.dat
            firstline = 0;
            printf("Using %d observations\n", n);

            // allocate memory for the input data
            // NOTE: A is col-major indexed i.e. A[i + n*j] = A_(i,j)
            A = realloc(A, n*n*sizeof(double));
            if (A == NULL) {
                fprintf(stderr, "malloc failed\n");
                return(1);
            }else{
                printf("Memory allocated for input and design matrices\n\n");
            }

        }else{

            /* printf("Read line (len %zu) : ", read); */
            /* printf("%s", line); */
            pch = strtok(line," "); // parse line into matrix A's j'th col
            double e = 0.0;
            j = 0;
            while (pch != NULL) {
                e = atof(pch);
                printf("%d,%d: %f\t", i, j, e);
                pch = strtok(NULL, " ");
                j++;
            }

            i++;
            printf("\n");
        }
    }
    fclose(ifp);
    printf("\n");

    return 0;
}
// END MAIN
