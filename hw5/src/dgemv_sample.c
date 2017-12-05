// test dgemv
double *X = (double*) malloc(n*n* sizeof(double));
double *y = (double*) malloc(n* sizeof(double));
double *y2 = (double*) malloc(n* sizeof(double));
for(i=0; i<n; i++){
    for(j=0; j<n; j++){
        X[i + j*n] = 1.0;
    }
    y[i] = 1.0;
}
printf("X:\n");
for(i=0; i<n; i++){
    for(j=0; j<n; j++){
        printf("%f\t", X[i + j*n]);
    }
    printf("\n");
}
printf("\ny:\n");
for(i=0; i<n; i++){
    printf("%f\n", y[i]);
}

TRANS = 'N'; // transpose the matrix A
M = n; // rows of A (number of obs)
N = n; // columns of A (degree of polynom)
LDA = n; // leading dimension of A
INCB = 1;
ALPHA = 1.0;
BETA = 0.0;

dgemv_(&TRANS, &M, &N,
    &ALPHA, X, &LDA, y, &INCB, &BETA,
    y2, &INCB);

printf("\ny2:\n");
for(i=0; i<n; i++){
    printf("%f\n", y2[i]);
}
printf("\n");

