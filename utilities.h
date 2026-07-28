#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h>
#include <time.h>
#include <unistd.h>

typedef struct
{
    struct timespec start;
    struct timespec end;
    double time;
} timertype;

double randunif(double a, double b);
void set_seed(int seed);
void starttimer(timertype *timer);
void endtimer(timertype *timer);
void dvprint(double *v, int m);
void dmprint(double *A, int m, int n);
void zvprint(double complex *z, int m);
void zmprint(double complex *A, int m, int n);
void dvfprint(double *v, int m, FILE *fout);
void zvfprint(double complex *z, int m, FILE *fout);
void dmfprint(double *A, int m, int n, FILE *fout);
void zmfprint(double complex *A, int m, int n, FILE *fout);
void dsydiagonalize(double *A, double *V, double *w, int m);
void zhediagonalize(double complex *A, double complex *V, double *w, int m);
void dsydiagexpmv(double *V, double *w, double a, double *v, int m, double *work);
void zhediagexpmv(double complex *V, double *w, double complex a, double complex *v, int m, double complex *work);
void dnormalize(double *v, int m);
void znormalize(double complex *v, int m);

// BLAS function prototypes -----------------------------------
void dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
void zgemv_(char *trans, int *m, int *n, double complex *alpha, double complex *a, int *lda, double complex *x, int *incx, double complex *beta, double complex *y, int *incy);
void dscal_(int *n, double *da, double *dx, int *incx);
double dnrm2_(int *n, double *x, int *incx);
double dznrm2_(int *n, double complex *x, int *incx);
void zdscal_(int *n, double *da, double complex *zx, int *incx);
// ------------------------------------------------------------

// LAPACK function prototypes ---------------------------------
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
void zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda, double *w, double complex *work, int *lwork, double *rwork, int *info);
void dlacpy_(char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);
void zlacpy_(char *uplo, int *m, int *n, double complex *a, int *lda, double complex *b, int *ldb);
// ------------------------------------------------------------

#endif