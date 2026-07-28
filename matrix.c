#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h>
#include "matrix.h"

// function prototype for BLAS's zscal
// notice the trailing _
void zscal_(int *n, double complex *za, double complex *zx, int *incx);
double dznrm2_(int *n, double complex *zx, int *incx);
void zgemv_ (
    char *trans,
    int *m,
    int *n,
    double complex *alpha,
    double complex *a,
    int *lda,
    double complex *x,
    int *incx,
    double complex *beta,
    double complex *Y,
    int *incy);

void zgemm_(
    char *transa, char *transb,
    int *m, int *n, int *k,
    double complex *alpha,
    double complex *a, int *lda,
    double complex *b, int *ldb, 
    double complex *beta,
    double complex *c, int *ldc);

void zheev_(char *jobz, char *uplo,
    int *n,
    double complex *a, int *lda,
    double *w, 
    double complex *work, int *lwork, double* rwork, int *info);

double normalize(double complex *vector, int length){
    int one = 1;
    return dznrm2_(&length, vector, &one);
}

double complex *mv_mult(double complex *matrix, double complex *vector, int m, int n){
    char trans = 'n';
    int inc = 1;
    double complex alpha = 1.0;
    double complex beta = 0.0;
    double complex *result = calloc(m, sizeof(double complex));
    zgemv_(&trans, &m, &n, &alpha, matrix, &m, vector, &inc, &beta, result, &inc);
    return result;
}

double complex *mm_mult(double complex *matrix1, double complex *matrix2, int m, int n, int k){
    char trans = 'n';
    double complex alpha = 1.0;
    double complex beta = 0.0;
    double complex *result = calloc(m*n, sizeof(double complex));
    zgemm_(&trans, &trans, &m, &n, &k, &alpha, matrix1, &m, matrix2, &k, &beta, result, &m);
    return result;
}

double complex *eigen(double complex *matrix, int n){
    char jobz = 'V';
    char uplo = 'U';
    int info;
    int lwork = -1;
    double complex work_query;
    double complex *work = NULL;
    double *w;
    double *rwork;
    double complex *result;

    if (matrix == NULL || n < 1) {
        return NULL;
    }

    w = calloc(n, sizeof(double));
    rwork = calloc(n > 1 ? 3*n - 2 : 1, sizeof(double));
    result = calloc(n, sizeof(double complex));
    if (w == NULL || rwork == NULL || result == NULL) {
        free(w);
        free(rwork);
        free(result);
        return NULL;
    }

    /*
     * Ask LAPACK for the optimal complex workspace size.  On success,
     * zheev overwrites matrix with orthonormal eigenvectors and writes the
     * (real) eigenvalues to w in ascending order.
     */
    zheev_(&jobz, &uplo, &n, matrix, &n, w, &work_query, &lwork,
           rwork, &info);
    if (info != 0) {
        fprintf(stderr, "zheev workspace query failed (info = %d)\n", info);
        free(w);
        free(rwork);
        free(result);
        return NULL;
    }

    lwork = (int) creal(work_query);
    work = calloc(lwork, sizeof(double complex));
    if (work == NULL) {
        free(w);
        free(rwork);
        free(result);
        return NULL;
    }

    zheev_(&jobz, &uplo, &n, matrix, &n, w, work, &lwork, rwork, &info);
    free(work);
    free(rwork);

    if (info != 0) {
        if (info < 0) {
            fprintf(stderr, "zheev: argument %d is invalid\n", -info);
        } else {
            fprintf(stderr, "zheev failed to converge (info = %d)\n", info);
        }
        free(w);
        free(result);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        result[i] = w[i];
    }
    free(w);
    return result;

}
