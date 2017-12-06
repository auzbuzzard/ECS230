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

/* double dnrm2_() */
// END PROTOTYPES

// BEGIN VARIABLES
// for general configuration
/* const char data_fn[] = "../data/data.dat"; // location of the input data */
/* const char data_fn[] = "../data/test_A.dat"; // location of the input data */
const char data_fn[] = "../data/bryan_leise_2006_fig1.dat";
int i, j, k, itt; // loop counters
int n; // size of graph's vertex set
int firstline; // determines whether reading first line of data.dat or not

// file handling objects
FILE *ifp; // pointer to the file object containing the input data
char *mode = "r"; // opening mode for the file above
FILE * out_raw; // where to write raw output
FILE * gnuplotPipe;

// CBLAS/LAPACK inputs
char TRANS;
double ALPHA;
double BETA;
int M;
int N;
int LDA;
int LDB;
int LDC;
int INCA;
int INCB;

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
    double *b1 = (double*) malloc(n* sizeof(double)); // placeholder malloc
    double *b2 = (double*) malloc(n* sizeof(double)); // placeholder malloc

    // read file
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * pch;
    i = 0;
    while ((read = getline(&line, &len, ifp)) != -1){
        if(firstline == 1){

            /* printf("Read first line (len %zu) : ", read); */
            /* printf("%s\n", line); */
            n = atoi(line); // n_obs is first line of data.dat
            firstline = 0;
            /* printf("Using %d observations\n", n); */
            printf("n=%d\n", n);

            // allocate memory for the input data
            // NOTE: A is col-major indexed i.e. A[i + n*j] = A_(i,j)
            A = realloc(A, n*n*sizeof(double));
            b1 = realloc(b1, n*sizeof(double));
            b2 = realloc(b2, n*sizeof(double));

            if (A == NULL) {
                fprintf(stderr, "malloc failed\n");
                return(1);
            }else{
                /* printf("Memory allocated for input and design matrices\n\n"); */
            }

            // initialize A as 0s, b as arbitrary
            int i1, i2;
            for(i1=0;i1<n;i1++){
                for(i2=0;i2<n;i2++){
                    A[i1 + n*i2] = 0.0;
                    /* printf(A[i + n*j]); */
                }
                b1[i1] = 1.0;///n;
            }

        }else{

            // Calculate and install input into matrix A
            /* printf("Read line (len %zu) : ", read); */
            /* printf("%s", line); */
            pch = strtok(line," "); // parse line into matrix A's j'th col
            k = 0;
            firstline = 1;
            double nj;
            while (pch != NULL) {
                if(firstline == 1){
                    nj = atof(pch); // number of outlinks for node j
                    firstline = 0;
                    /* printf("Col %i connected to %d rows\n", i, (int) nj); */
                }else{
                    k = atof(pch); // node pointed to by j is k
                    /* printf("k,i:%d,%d:\n", k,i); */
                    A[k-1 + n*i] = 1.0/nj;
                    /* printf("A_%d,%d: %f\n", i, k-1, A[i + n*(k-1)]); */
                }

                pch = strtok(NULL, " ");
            }

            i++;
            /* printf("\n"); */
        }
    }
    fclose(ifp);
    printf("\n");

    // print A
    i = 0;
    j = 0;
    printf("A:\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("%f ", A[i + n*j]);
        }
        printf("\n");
    }

    // normalize b_0
    double norm_b;
    int INCX = 1;
    norm_b = cblas_dnrm2(n, b1, INCX);
    /* printf("\nnorm b:\n%f\tnormalizing...\n", norm_b); */
    for(i=0; i<n; i++){
        b1[i] = b1[i] / norm_b;
    }
    norm_b = cblas_dnrm2(n, b1, INCX);

    // print normalized b_0
    i = 0;
    printf("\nb0:\n");
    for(i=0;i<n;i++){
        printf("%f\n", b1[i]);
    }
    printf("\n");
    /* printf("\nnorm b:\n%f\n", norm_b); */

    // Power method
    // $b_{k+1} = A b_k \over \norm{A b_k}$ until convergence, determined as
    // follows: $\infty_norm{ b_{k+1} - b_k } < \epsilon = 1.0e-6$
    //
    // Using dgemv_(A, b)

    TRANS = 'N'; // transpose the matrix A
    M = n; // rows of A (number of obs)
    N = n; // columns of A (degree of polynom)
    LDA = n; // leading dimension of A
    INCB = 1;
    ALPHA = 1.0;
    BETA = 0.0;

    double epsilon = 1.0e-6;
    /* double epsilon = 1.0e-2; */
    double delta = 1.0;
    itt = 0;
    while(delta > epsilon){

        // normalize b1
        norm_b = cblas_dnrm2(n, b1, INCX);
        for(i=0; i<n; i++){
            b1[i] = b1[i] / norm_b;
        }
        /* printf("\nnorm b%d:%f\n", itt, norm_b); */

        // multiply b2 <- A b1
        dgemv_(&TRANS, &M, &N, &ALPHA,
                A, &LDA,
                b1, &INCB, &BETA,
                b2, &INCB);

        // normalize b2
        norm_b = cblas_dnrm2(n, b2, INCX);
        /* printf("\nb%d:\n", itt); */
        for(i=0;i<n;i++){
            b2[i] = b2[i] / norm_b;
            /* printf("%f\n", b2[i]); */
        }

        // calculate delta
        delta = 0.0;
        for(i=0; i<n; i++){
            delta = delta + fabs(b1[i] - b2[i]);
        }
        printf("delta %d: %f\n", itt, delta);

        // copy b1 <- b2
        for(i=0;i<n;i++){
            b1[i] = b2[i];
            /* printf("%f\n", b1[i]); */
        }

        itt++;
    }

    // scale eigenvector so L1 norm is 1
    norm_b = cblas_dasum(n, b1, 1);
    for(i=0; i<n; i++){
        b1[i] = b1[i] / norm_b;
    }

    printf("\nConvergence after %d iterations to dominant eigenvector b:\n",
            itt);
    for(i=0; i<n; i++){
        printf("%f\n", b1[i]);
    }

    return 0;
}
// END MAIN
