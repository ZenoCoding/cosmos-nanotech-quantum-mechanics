#include "utilities.h"


double randunif(double a, double b)
{
    //
    // Returns a random number between a and b.
    // Call like: double a = 0.0, b = 1.0;
    //            double r = randunif(a, b);
    //
    double r;
    r = (double) rand() / RAND_MAX;
    return a + (b - a) * r;
}

void set_seed(int seed)
{
    //
    // Sets the seet of the random number generator.
    // Call like: int seed = 123;
    //            set_seed(seed);
    //
    srand(seed);
}

void starttimer(timertype *timer)
{
    //
    // Starts a timer.
    // Call like: starttimer(&timer)
    //
    clock_gettime(CLOCK_MONOTONIC, &(*timer).start);
}

void endtimer(timertype *timer)
{
    //
    // Ends a timer and records the time.
    // Call like: endtimer(&timer)
    // The elapsed time (in seconds) can be accessed like: timer.time
    //
    clock_gettime(CLOCK_MONOTONIC, &(*timer).end);
    (*timer).time = (*timer).end.tv_sec - (*timer).start.tv_sec +
    ((*timer).end.tv_nsec - (*timer).start.tv_nsec) / 1e9;
}

void dvprint(double *v, int m)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        printf("%12.6lf ", v[i]);
    }
    printf("\n");
}

void dvfprint(double *v, int m, FILE *fout)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        fprintf(fout, "%12.6lf ", v[i]);
    }
    fprintf(fout, "\n");
}

void zvprint(double complex *z, int m)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        printf("%12.6lf %+12.6lfi ", creal(z[i]), cimag(z[i]));
    }
    printf("\n");
}

void zvfprint(double complex *z, int m, FILE *fout)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        fprintf(fout, "%12.6lf %+12.6lfi ", creal(z[i]), cimag(z[i]));
    }
    fprintf(fout, "\n");
}

void dmprint(double *A, int m, int n)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        for (int j = 0 ; j < n ; j = j + 1)
        {
            printf("%12.6lf ", A[i + m*j]);
        }
        printf("\n");
    }
}

void dmfprint(double *A, int m, int n, FILE *fout)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        for (int j = 0 ; j < n ; j = j + 1)
        {
            fprintf(fout, "%12.6lf ", A[i + m*j]);
        }
        fprintf(fout, "\n");
    }
}

void zmprint(double complex *A, int m, int n)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        for (int j = 0 ; j < n ; j = j + 1)
        {
            printf("%12.6lf %+12.6lfi ", creal(A[i + m*j]), cimag(A[i + m*j]));
        }
        printf("\n");
    }
}

void zmfprint(double complex *A, int m, int n, FILE *fout)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        for (int j = 0 ; j < n ; j = j + 1)
        {
            fprintf(fout, "%12.6lf %+12.6lfi ", creal(A[i + m*j]), cimag(A[i + m*j]));
        }
        fprintf(fout, "\n");
    }
}

void dsydiagonalize(double *A, double *V, double *w, int m)
{
    //
    // "Convenient" double symmetric matrix diagonalizer.
    // Diagonalizes the symmetric matrix A of doubles (without altering it)
    // and stores the eigenvectors and eigenvalues in V and w:
    //
    //          eigenvalue w[i] <----> eigenvector column i of V
    //
    // On input:
    // A (double *): m x m matrix of doubles to diagonalize. At least the lower
    //               triangular part of A should be filled in.
    // V (double *): m x m matrix of doubles to store the eigenvectors of A in.
    //               Should already be allocated prior to calling.
    // w (double *): m long vector of doubles to store the eigenvalues of w in.
    //               Should already be allocated prior to calling.
    // m (int)     : dimension of A.
    //
    // On output:
    // A: unchanged.
    // V: column i contains the ith eigenvector of A.
    // w: entry i contains the ith eigenvalue of A.
    // m: unchanged.
    //

    char jobz = 'V', uplo = 'L';
    int info = 0, lwork = (128 + 2) * m;
    double *work = malloc(lwork * sizeof(double));

    // Copy V = lowertri(A)
    dlacpy_(&uplo, &m, &m, A, &m, V, &m);

    // Diagonalize
    dsyev_(&jobz, &uplo, &m, V, &m, w, work, &lwork, &info);

    // Deallocate workspace
    free(work);
}

void zhediagonalize(double complex *A, double complex *V, double *w, int m)
{
    //
    // "Convenient" zomplex Hermitian matrix diagonalizer.
    // Diagonalizes the Hermitian matrix A of zomplexes (without altering it)
    // and stores the eigenvectors and eigenvalues in V and w:
    //
    //          eigenvalue w[i] <----> eigenvector column i of V
    //
    // On input:
    // A (double complex *): m x m matrix of zomplexes to diagonalize. At least the lower
    //                       triangular part of A should be filled in.
    // V (double complex *): m x m matrix of zomplexes to store the eigenvectors of A in.
    //                       Should already be allocated prior to calling.
    // w (double *)        : m long vector of doubles to store the eigenvalues of w in.
    //                       Should already be allocated prior to calling.
    // m (int)             : dimension of A.
    //
    // On output:
    // A: Unchanged.
    // V: Column i contains the ith eigenvector of A.
    // w: Entry i contains the ith eigenvalue of A.
    // m: Unchanged.
    //

    char jobz = 'V', uplo = 'L';
    int info = 0, lwork = (128 + 2) * m;
    double complex *work = malloc(lwork * sizeof(double complex));
    double *rwork = malloc((3*m - 2) * sizeof(double));

    // Copy V = lowertri(A)
    zlacpy_(&uplo, &m, &m, A, &m, V, &m);

    // Diagonalize
    zheev_(&jobz, &uplo, &m, V, &m, w, work, &lwork, rwork, &info);

    // Deallocate workspaces
    free(work); free(rwork);
}

void dsydiagexpmv(double *V, double *w, double a, double *v, int m, double *work)
{
    //
    // Updates:
    //
    //     v = exp(a*A) * v
    //
    // where A (a symmetric matrix of doubles) has an eigendecomposition given by V and w:
    //
    //     A = V * w * V**T
    //
    // (see dsydiagonalize, V and w calculated from there work here).
    //
    // On input:
    // V (double *)   : m x m matrix of doubles containing the eigenvectors of A.
    //                  Column i is eigenvector i of A.
    // w (double *)   : m long vector of doubles containing the eigenvalues of A.
    //                  Entry i is the ith eigenvalue of A.
    // a (double)     : Scalar a in exp(a*A).
    // v (double *)   : m long vector of doubles to update v = exp(a*A) * v.
    // m (int)        : dimension of A and v.
    // work (double *): Workspace vector of doubles of length m.
    //                  Should already be allocated prior to calling.
    //
    // On output:
    // V   : Unchanged.
    // w   : Unchanged.
    // a   : Unchanged.
    // v   : Updated to v = exp(a*A) * v
    // m   : Unchanged.
    // work: Changed.
    //

    int one = 1;
    char trans;
    double done = 1.0, dzero = 0.0;

    // 1. work = V**T * v
    trans = 'T';
    dgemv_(&trans, &m, &m, &done, V, &m, v, &one, &dzero, work, &one);

    // 2. work = exp(a * w) * work
    for (int i = 0 ; i < m ; i = i + 1)
    {
        work[i] = exp(a * w[i]) * work[i];
    }

    // 3. v = V * work
    trans = 'N';
    dgemv_(&trans, &m, &m, &done, V, &m, work, &one, &dzero, v, &one);
}

void zhediagexpmv(double complex *V, double *w, double complex a, double complex *v, int m, double complex *work)
{
    //
    // Updates:
    //
    //     v = exp(a*A) * v
    //
    // where A (a Hermitian matrix of zomplexes) has an eigendecomposition given by V and w:
    //
    //     A = V * w * V**H
    //
    // (see zhediagonalize, V and w calculated from there work here).
    //
    // On input:
    // V (double complex *)   : m x m matrix of doubles containing the eigenvectors of A.
    //                          Column i is eigenvector i of A.
    // w (double *)           : m long vector of doubles containing the eigenvalues of A.
    //                          Entry i is the ith eigenvalue of A.
    // a (double complex)     : Scalar a in exp(a*A).
    // v (double complex *)   : m long vector of doubles to update v = exp(a*A) * v.
    // m (int)                : dimension of A and v.
    // work (double complex *): Workspace vector of doubles of length m.
    //                          Should already be allocated prior to calling.
    //
    // On output:
    // V   : Unchanged.
    // w   : Unchanged.
    // a   : Unchanged.
    // v   : Updated to v = exp(a*A) * v
    // m   : Unchanged.
    // work: Changed.
    //

    int one = 1;
    char trans;
    double complex zone = 1.0, zzero = 0.0;

    // 1. work = V**H * v
    trans = 'C';
    zgemv_(&trans, &m, &m, &zone, V, &m, v, &one, &zzero, work, &one);

    // 2. work = exp(a * w) * work
    for (int i = 0 ; i < m ; i = i + 1)
    {
        work[i] = exp(a * w[i]) * work[i];
    }

    // 3. v = V * work
    trans = 'N';
    zgemv_(&trans, &m, &m, &zone, V, &m, work, &one, &zzero, v, &one);
}

void dnormalize(double *v, int m)
{
    //
    // Normalizes the vector v of doubles.
    //
    // On input:
    // v (double *): Vector of doubles to normalize.
    // m (int)     : Length of the vector v.
    //
    // On output:
    // v: Changed. Normalized.
    // m: Unchanged.
    //

    int one = 1;
    double factor = 1.0 / dnrm2_(&m, v, &one);
    dscal_(&m, &factor, v, &one);
}

void znormalize(double complex *v, int m)
{
    //
    // Normalizes the vector v of zomplexes.
    //
    // On input:
    // v (double complex *): Vector of zomplexes to normalize.
    // m (int)             : Length of the vector v.
    //
    // On output:
    // v: Changed. Normalized.
    // m: Unchanged.
    //
    int one = 1;
    double factor = 1.0 / dznrm2_(&m, v, &one);
    zdscal_(&m, &factor, v, &one);
}