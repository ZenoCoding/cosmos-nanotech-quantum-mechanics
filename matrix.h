#ifndef MATRIX_H
#define MATRIX_H

#include <complex.h>

double normalize(double complex *vector, int length);

/*
 * Returns matrix*vector for an m-by-n matrix.
 * The caller must free the returned m-element vector.
 */
double complex *mv_mult(
    double complex *matrix,
    double complex *vector,
    int m,
    int n
);

/*
 * Multiplies an m-by-k matrix by a k-by-n matrix.
 * The caller must free the returned m-by-n matrix.
 */
double complex *mm_mult(
    double complex *matrix1,
    double complex *matrix2,
    int m,
    int n,
    int k
);

/*
 * Diagonalizes an n-by-n Hermitian matrix.
 * The input matrix is overwritten with eigenvectors stored by column.
 * The caller must free the returned eigenvalue array.
 */
double complex *eigen(double complex *matrix, int n);

#endif
