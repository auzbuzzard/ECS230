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
int i, j, k; // loop counters
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
    double *b = (double*) malloc(n* sizeof(double)); // placeholder malloc

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
            b = realloc(b, n*sizeof(double));

            if (A == NULL) {
                fprintf(stderr, "malloc failed\n");
                return(1);
            }else{
                printf("Memory allocated for input and design matrices\n\n");
            }

            // initialize A as 0s, b as arbitrary
            int i1, i2;
            for(i1=0;i1<n;i1++){
                for(i2=0;i2<n;i2++){
                    A[i1 + n*i2] = 0.0;
                    /* printf(A[i + n*j]); */
                }
                b[i1] = 1.0/(float) (i1+1);
            }

        }else{

            // Calculate and install input into matrix A
            /* printf("Read line (len %zu) : ", read); */
            /* printf("%s", line); */
            pch = strtok(line," "); // parse line into matrix A's j'th col
            k = 0;
            // i is column, k is row
            firstline = 1;
            float nj;
            while (pch != NULL) {
                if(firstline == 1){
                    nj = atof(pch); // number of outlinks for node j
                    firstline = 0;
                    printf("Col %i connected to %d rows\n", i, (int) nj);
                }else{
                    k = atof(pch); // node pointed to by j is k
                    printf("k,i:%d,%d:\n", k,i);
                    A[i + n*(k-1)] = 1.0/nj;
                    printf("A_%d,%d: %f\n", i, k-1, A[i + n*(k-1)]);
                }

                pch = strtok(NULL, " ");
            }

            i++;
            printf("\n");
        }
    }
    fclose(ifp);
    printf("\n");

    // print A and b
    i = 0;
    j = 0;
    printf("A:\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("%f ", A[i + n*j]);
        }
        printf("\n");
    }
    i = 0;
    printf("b:\n");
    for(i=0;i<n;i++){
        printf("%f\t", b[i]);
    }

    // Power method
    // $b_{k+1} = A b_k \over \norm{A b_k}$ until convergence, determined as
    // follows: $\infty_norm{ b_{k+1} - b_k } < \epsilon = 1.0e-6$
    //
    // Using dgemv_(A, b);

    return 0;
}
// END MAIN
