#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "matrix.h"
#include "utilities.h"

int main(void)
{
    /*
     * Write your tight-binding / atomic-overlap analysis here.
     *
     * Matrix helpers available from matrix.h:
     *   normalize(vector, length)
     *   mv_mult(matrix, vector, rows, columns)
     *   mm_mult(matrix1, matrix2, rows, columns, shared_dimension)
     *   eigen(hermitian_matrix, size)
     *
     * Printing helpers available from utilities.h:
     *   dvprint, zvprint, dmprint, zmprint
     *
     * Matrices use column-major storage:
     *   matrix[row + number_of_rows * column]
     */

    int N = 1000;
    double EA = -5.0;
    double EB = -4.0;
    double t = .6;

    double complex *matrix = calloc(N*N, sizeof *matrix);
    if (matrix == NULL) {
        fprintf(stderr, "Could not allocate the Hamiltonian.\n");
        return EXIT_FAILURE;
    }

    for(int i = 0; i < N; i++){
        if(i > 0) matrix[i+i*N-1] = -t;
        matrix[i+i*N] = i % 2 == 0 ? EA : EB;
        if(i < N - 1) matrix[i+i*N+1] = -t;
    }

    // zmprint(matrix, N, N);

    printf("\n");

    double complex *eigenvalues = eigen(matrix, N);
    if (eigenvalues == NULL) {
        free(matrix);
        return EXIT_FAILURE;
    }

    // zvprint(eigenvalues, N);

    const char *eigenvalue_file = "labs/12/data/eigenvalues2.dat";
    FILE *output = fopen(eigenvalue_file, "w");
    if (output == NULL) {
        perror(eigenvalue_file);
        free(matrix);
        free(eigenvalues);
        return EXIT_FAILURE;
    }

    for (int i = 0; i < N; i++) {
        fprintf(output, "%.17g\n", creal(eigenvalues[i]));
    }

    if (fclose(output) != 0) {
        perror(eigenvalue_file);
        free(matrix);
        free(eigenvalues);
        return EXIT_FAILURE;
    }

    printf("\nWrote %d eigenvalues to %s\n", N, eigenvalue_file);

    free(matrix);
    free(eigenvalues);

    return EXIT_SUCCESS;
}
