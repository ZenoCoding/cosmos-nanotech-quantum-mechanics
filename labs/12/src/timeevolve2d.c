//
// From the repository root, compile with:
//
//     make lab12
//
// Run like:
//
//     ./build/12/timeevolve2d 8
//
// The optional command-line argument is the side length of the square lattice.
// Each line of timeevolution2d.txt contains one time slice, flattened in
// row-major lattice order:
//
//     (0,0), (1,0), ..., (L-1,0), (0,1), ..., (L-1,L-1).
//
#include "utilities.h"
#include <errno.h>
#include <limits.h>
#include <stdint.h>

double modulussq(double complex z);
void probability(const double complex *psi, double *prob, int n);

int main(int argc, char **argv)
{
    //
    // L     : side length of the L-by-L square lattice.
    // J     : nearest-neighbor hopping strength.
    // steps : number of time steps.
    // tmax  : final simulation time.
    //
    int          L     = 8;
    const double J     = 1.0;
    const int    steps = 100;
    const double tmax  = 10.0 / J;
    const char  *fname = "labs/12/data/timeevolution2d.txt";

    if (argc > 2)
    {
        fprintf(stderr, "Usage: %s [lattice-side-length]\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (argc == 2)
    {
        char *end = NULL;
        errno = 0;
        long requested_L = strtol(argv[1], &end, 10);

        if (errno != 0 || end == argv[1] || *end != '\0' ||
            requested_L < 2 || requested_L > INT_MAX)
        {
            fprintf(stderr, "The lattice side length must be an integer >= 2.\n");
            return EXIT_FAILURE;
        }
        L = (int) requested_L;
    }

    if ((size_t) L > SIZE_MAX / (size_t) L)
    {
        fprintf(stderr, "The requested lattice is too large.\n");
        return EXIT_FAILURE;
    }

    size_t site_count = (size_t) L * (size_t) L;
    if (site_count > INT_MAX ||
        site_count > SIZE_MAX / site_count ||
        site_count * site_count > SIZE_MAX / sizeof(double complex))
    {
        fprintf(stderr, "The requested lattice is too large.\n");
        return EXIT_FAILURE;
    }

    int N = (int) site_count;
    const double dtime = tmax / steps;

    double complex *H    = calloc(site_count * site_count, sizeof *H);
    double complex *V    = calloc(site_count * site_count, sizeof *V);
    double         *w    = calloc(site_count, sizeof *w);
    double complex *psi  = calloc(site_count, sizeof *psi);
    double         *prob = calloc(site_count, sizeof *prob);
    double complex *work = calloc(site_count, sizeof *work);

    if (H == NULL || V == NULL || w == NULL || psi == NULL ||
        prob == NULL || work == NULL)
    {
        fprintf(stderr, "Could not allocate arrays for a %d-by-%d lattice.\n",
                L, L);
        free(H); free(V); free(w); free(psi); free(prob); free(work);
        return EXIT_FAILURE;
    }

    //
    // Build the square-lattice tight-binding Hamiltonian with open boundaries:
    //
    //              H = -J sum_<i,j> (|i><j| + |j><i|).
    //
    // Only right and upper bonds are visited so that each bond is added once.
    // Matrices use column-major storage: H[row + N*column].
    //
    for (int y = 0; y < L; y++)
    {
        for (int x = 0; x < L; x++)
        {
            int site = x + L*y;

            if (x + 1 < L)
            {
                int right = site + 1;
                H[site  + N*right] = -J;
                H[right + N*site]  = -J;
            }

            if (y + 1 < L)
            {
                int up = site + L;
                H[site + N*up]   = -J;
                H[up   + N*site] = -J;
            }
        }
    }

    zhediagonalize(H, V, w, N);

    FILE *fout = fopen(fname, "w");
    if (fout == NULL)
    {
        perror(fname);
        free(H); free(V); free(w); free(psi); free(prob); free(work);
        return EXIT_FAILURE;
    }

    //
    // Begin with one excitation localized at the center of the lattice.
    // For even L, this selects the upper-right site of the central four.
    //
    int center = (L / 2) + L*(L / 2);
    psi[center] = 1.0;

    znormalize(psi, N);
    probability(psi, prob, N);
    dvfprint(prob, N, fout);

    for (int step = 1; step <= steps; step++)
    {
        zhediagexpmv(V, w, -I*dtime, psi, N, work);
        znormalize(psi, N);
        probability(psi, prob, N);
        dvfprint(prob, N, fout);
    }

    printf("Square lattice: %d x %d (%d sites), J = %.6f\n", L, L, N, J);
    printf("Initial site: (%d, %d)\n", L / 2, L / 2);
    printf("Evolved from t = 0 to t = %.6f in %d steps\n", tmax, steps);
    printf("Wrote %d time slices to %s\n", steps + 1, fname);

    fclose(fout);
    free(H); free(V); free(w); free(psi); free(prob); free(work);

    return EXIT_SUCCESS;
}

double modulussq(double complex z)
{
    return creal(z * conj(z));
}

void probability(const double complex *psi, double *prob, int n)
{
    for (int i = 0; i < n; i++)
    {
        prob[i] = modulussq(psi[i]);
    }
}
