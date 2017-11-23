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
// 2. http://www.netlib.org/clapack/cblas/dgemm.c
// 3. https://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html
// 4  http://www.math.utah.edu/software/lapack/lapack-d/dpotrf.html
// 5. https://stackoverflow.com/questions/3521209/making-c-code-plot-a-graph-automatically

#include "stdio.h"
#include "stdlib.h"
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
const char data_fn[] = "../data/data.dat"; // location of the input data
int dmax = 4, dmin = 3; // degree of fit polynomial from commandline
int d; // counter for the degree of the polynomial to fit
int i, j; // loop counters
float x, y; // for reading x and y from data_fn
int n; // size of the input data (n by 2 i.e. n xs and n ys)

// file handling objects
FILE *ifp; // pointer to the file object containing the input data
FILE * out_raw;
FILE * out_fit;
FILE * out_coef;
FILE * out_designMx;
FILE * gnuplotPipe;
char *mode = "r"; // opening mode for the file above

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

// stats calculations
double mu;
double sse;
double var;
double r;
// END VARIABLES

// BEGIN MAIN
int main(int argc, char** argv)
{

    for(d=dmin; d<dmax+1; d++){

    // read data into memory from disk
    ifp = fopen(data_fn, mode);
    if (ifp == NULL) {
        fprintf(stderr, "Cannot open file\n"); // error if it can't open
        exit(1);
    }

    // print input contents back out to stdIO and count lines
    n = 0;
    while (fscanf(ifp, "%f %f", &x, &y) != EOF) { // for each line in the file
        n++; // count the number of lines
    }
    fclose(ifp);

    // allocate memory for the input data
    double *X = (double*) malloc(n*dmax* sizeof(double));
    // NOTE: X is col-major indexed i.e. X[i + n*j] = X_(i,j)
    double *Y = (double*) malloc(n* sizeof(double));
    if ((X == NULL) || (Y == NULL)) {
        fprintf(stderr, "malloc failed\n");
        return(1);
    }

    // read input data into the X and Y arrays
    ifp = fopen(data_fn, mode);
    i = 0;
    while (fscanf(ifp, "%f %f", &x, &y) != EOF) { // for each line in the file
        for(j = 0; j < dmax+1; j++){
            X[i + j*n] = pow(x, j); // each col of X is ((x_i)^j)_{i=1..n}, j<=d
        }
        Y[i] = y;
        /* printf("x=%f,\ty=%f\n", x, y); */
        i++;
    }
    fclose(ifp);

    // print the X and Y matrices
    printf("Y\t\t"); // BEGIN formatting printing header
    for(j = 0; j < d; j++){
        printf("X^%d\t\t", j);
    }
    printf("\n"); // END formatting printing header

    for(i = 0; i < n; i++){ // BEGIN print X and Y
        printf("%f\t", Y[i]);
        for(j = 0; j < d; j++){
            printf("%f\t", X[i + j*n]);
        }
        printf("\n");
    } // END print X and Y


    // compute A=X^T X using BLAS::dgemm()
    // initialize
    TRANSA = 'T'; // transpose the matrix X for (X^T X)
    TRANSB = 'N';
    ALPHA = 1.0; // dgemm(A,B,C) : C <- alpha*AB + beta*C
    BETA = 0.0;
    M = d; // rows of X^T
    N = d; // columns of X
    K = n; // columns of X^T and rows of X
    LDA = n; // leading dimension of X
    LDB = n; // leading dimension of X
    LDC = d; // leading dimension of A

    // allocate memory for A
    double *A = (double*) malloc(d*d* sizeof(double));
    if (A == NULL) {
        fprintf(stderr, "malloc failed\n");
        return(1);
    }

    // compute A = X^T X
    dgemm_(&TRANSA, &TRANSB,
        &M, &N, &K,
        &ALPHA, X, &LDA, X, &LDB, &BETA,
        A, &LDC);

    // print A
    printf("\nA = X^T X\n");
    for(i = 0; i < d; i++){
        for(j = 0; j < d; j++){
            printf("%f\t", A[i + j*d]);
        }
        printf("\n");
    }

    // Compute P = Xt y using BLAS::dgemm()
    // initialize
    TRANS = 'T'; // transpose the matrix X for (X^T y)
    M = n; // rows of X (number of obs)
    N = d; // columns of X (degree of polynom)
    LDA = n; // leading dimension of X
    INCY = 1; // increment for the input vector
    INCP = 1; // increment for the input vector

    // allocate memory for P = X^T y
    double *P = (double*) malloc(d* sizeof(double));
    if (P == NULL) {
        fprintf(stderr, "malloc failed\n");
        return(1);
    }

    // compute P = X^T y using dgemv
    dgemv_(&TRANS, &M, &N,
        &ALPHA, X, &LDA, Y, &INCY, &BETA,
        P, &INCP);

    // print P
    printf("\nP = X^T y\n");
    for(i = 0; i < d; i++){
        printf("%f\t", P[i]);
    printf("\n");
    }

    // Compute Choelsky decomposition XTX = A = LTL i.e. get L
    // allocate memory for L
    double *L = (double*) malloc(d*d* sizeof(double));
    if (L == NULL) {
        fprintf(stderr, "malloc failed\n");
        return(1);
    }

    // set up L = A for subsequent diagonalization
    for(i=0; i<d*d; i++){
        L[i] = A[i];
    }

    // set up the diagonalization and run it
    char UPLO = 'L'; // lower triangular
    int ORD = d; // L is order d
    int LDD = d; // leading dimension of matrix
    int INFO = 0; // output for diagonalization
    dpotrf_(&UPLO, &ORD, L, &LDD, &INFO);

    if(INFO != 0){
        fprintf(stderr, "Choelsky decomosition failed\n");
        return(1);
    }

    // set the above-diagonal entries to 0 - this is just an aesthetic thing
    for(i=0; i<d; i++){
        for(j=0; j<d; j++){
            if(j>i){
                L[i + d*j] = 0.0;
            }
        }
    }

    // print resultant diagonalized matrix L
    printf("\nL = Chol(A)\n");
    for(i = 0; i < d; i++){
        for(j = 0; j < d; j++){
            printf("%f\t", L[i + j*d]);
        }
        printf("\n");
    }

    // Solve for b (Xb=y) using BLAS::dtrsm()
    //  first, some algebra:
    //  Xb = y -> XtX b = Xt y -> Ab = p -> LLt b = p
    //  let Ltb = q, then solve L q = p using dtrsm_() for q
    //  upon solving for q, solve Ltb = q for b using dtrsm_()

    // set up the first triangular solve for L Q = P and run it
    SIDE = 'L'; // left, i.e. L q = p rather than q L = p
    UPLO = 'L'; // L is lower triangular
    TRANSA = 'N'; // solving L q = b not Lt
    DIAG = 'N'; // L is not unit diagonal (using Choelsky not LDL decomp)
    M = d; // rows of resultant vector P
    N = 1; // columns of resultant vector P
    ALPHA = 1.0; // coef for the rhs vector P
    LDA = d; // leading dimension on the lhs multiplier L
    LDB = d; // leading dimension of the rhs vector P

    // allocate a copy of P for the solution Q
    double *Q = (double*) malloc(d* sizeof(double));
    if (Q == NULL) {
        fprintf(stderr, "malloc failed\n");
        return(1);
    }
    for(i=0; i<d; i++){ // copy P -> Q
        Q[i] = P[i];
    }

    dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, // calculate Q = L^{-1} P
            &M, &N, &ALPHA, L, &LDA, Q, &LDB);

    // print the solution Q = L^{-1} P
    printf("\nQ = L^{-1} P\n");
    for(i = 0; i < d; i++){
        printf("%f\n", Q[i]);
    }

    // solve the second system
    SIDE = 'L'; // left, i.e. Lt
    UPLO = 'L'; // L is lower triangular
    TRANSA = 'T'; // solving Lt b = q
    DIAG = 'N'; // L is not unit diagonal (using Choelsky not LDL decomp)
    M = d; // rows of rhs vector Q
    N = 1; // columns of rhs vector Q
    ALPHA = 1.0; // coef for the rhs vector Q
    LDA = d; // leading dimension on the lhs multiplier Lt
    LDB = d; // leading dimension of the rhs vector Q

    // allocate a copy of Q for the solution B
    double *B = (double*) malloc(d* sizeof(double));
    if (B == NULL) {
        fprintf(stderr, "malloc failed\n");
        return(1);
    }
    for(i=0; i<d; i++){ // copy B <- Q; Lt B = Q
        B[i] = Q[i]; // B is overwritten with result Lt^{-1}Q
    }

    dtrsm_(&SIDE, &UPLO, &TRANSA, &DIAG, // calculate B = Lt^{-1}Q
            &M, &N, &ALPHA, L, &LDA, B, &LDB);

    // print the solution B = Lt^{-1} Q
    printf("\nB = L^T^{-1} Q\n");
    for(i = 0; i < d; i++){
        printf("%f\n", B[i]);
    }

    // Generate predictions using the fitted polynomial represented by B
    // allocate space for the predictions
    double *Yhat = (double*) malloc(n* sizeof(double));
    if (Yhat == NULL) {
        fprintf(stderr, "malloc failed\n");
        return(1);
    }

    TRANS = 'N'; // do not 'N' transpose the matrix X for (XB = yhat)
    M = n; // X is n rows by d columns
    N = d;
    LDA = n; // leading dimension of X
    ALPHA = 1.0; // coefficient for Y in Y := alpha AX + beta y
    BETA = 0.0;
    INCB = 1; // increment for the vector B
    INCY = 1; // increment for the resultant vector Yhat

    dgemv_(&TRANS, &M, &N,
        &ALPHA, X, &LDA, B, &INCB, &BETA,
        Yhat, &INCY);

    // print the solution Yhat = XB
    printf("\nYhat = XB\n");
    for(i = 0; i < n; i++){
        printf("%f\n", Yhat[i]);
    }

    // Write output to disk for subsequent analysis
    out_raw = fopen("../report/poly_raw.dat", "w");
    out_fit = fopen("../report/poly_fit.dat", "w");
    out_coef = fopen("../report/poly_coef.dat", "w");
    out_designMx = fopen("../report/poly_designMx.dat", "w");
    for (i=0; i < n; i++) {
        fprintf(out_raw, "%lf %lf\n", X[i + n], Y[i]);
        fprintf(out_fit, "%lf %lf\n", X[i + n], Yhat[i]);
        for(j=0; j<d; j++){
            fprintf(out_designMx, "%lf ", X[i + j*n]);
        }
        fprintf(out_designMx, "\n");
    }
    for (i=0; i < d; i++){
        fprintf(out_coef, "%lf\n", B[i]);
    }
    fclose(out_raw);
    fclose(out_fit);
    fclose(out_coef);
    fclose(out_designMx);

    // Calculate R^2
    mu = 0.0;
    sse = 0.0;
    var = 0.0;

    for(i=0; i<n; i++){
        mu += Y[i]/n;
        sse += pow((Yhat[i] - Y[i]), 2);
    }
    for(i=0; i<n; i++){
        var += pow((Y[i] - mu), 2);
    }
    r = 1.0 - sse/var;

    // Gnuplot the resulting fit and the raw data
    gnuplotPipe = popen("gnuplot", "w");
    fprintf(gnuplotPipe, "set terminal jpeg\n");
    fprintf(gnuplotPipe, "set output '../report/plot.jpeg'\n");

    fprintf(gnuplotPipe, "set grid\n" );
    fprintf(gnuplotPipe,
        "set title 'Observed data and polynomial fit (d=%d, R2=%f)'\n",
        d, r);
    fprintf(gnuplotPipe, "set key left box\n" );
    fprintf(gnuplotPipe, "set xlabel 'X'\n" );
    fprintf(gnuplotPipe, "set ylabel 'Y'\n" );
    fprintf(gnuplotPipe, "set style data points\n" );
    fprintf(gnuplotPipe, "set pointsize 2\n" );
    fprintf(gnuplotPipe, "plot '../report/poly_raw.dat' title 'Input', ");
    fprintf(gnuplotPipe, "'../report/poly_fit.dat' title 'Fit'\n");

    pclose(gnuplotPipe);

    // Cleanup the matrices
    free(A);
    free(B);
    free(X);
    free(Y);
    free(P);
    free(Q);
    free(L);
    free(Yhat);

    } // END polynomial degree loop (d)

    return 0;
}
// END MAIN
