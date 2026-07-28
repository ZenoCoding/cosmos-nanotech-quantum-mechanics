/*
 * Simulated annealing for state transfer on a finite honeycomb lattice.
 *
 * The lattice uses the same brick-wall indexing as timehex.c.  All original
 * nearest-neighbor bonds remain present and their magnitudes are constrained
 * to [J_MIN, J_MAX], so the optimizer cannot reduce the lattice to a chain.
 *
 * For an even side length L, reflection through the vertical midline is a
 * graph automorphism.  Couplings on reflected bonds are kept equal, and the
 * transfer is optimized from (0, 0) to (L - 1, 0).
 *
 * Compile:
 *
 *     gcc -O3 -std=c11 -D_POSIX_C_SOURCE=200809L \
 *         annealhex.c utilities.c -llapack -lblas -lm -o annealhex
 *
 * Run (all arguments are optional):
 *
 *     ./annealhex [L] [iterations] [seed] [output_file]
 *
 * Example:
 *
 *     ./annealhex 6 20000 123 annealed_couplings.dat
 *
 * The output contains one bond per line:
 *
 *     x1 y1 x2 y2 coupling
 */

#include "utilities.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>

#define DEFAULT_L          6
#define DEFAULT_ITERATIONS 10000
#define DEFAULT_SEED       123
#define DEFAULT_OUTPUT     "labs/12/data/annealed_couplings.dat"

#define J_MIN          0.25
#define J_MAX          2.00
#define INITIAL_J      1.00
#define INITIAL_TIME   6.00
#define TIME_MIN       0.05
#define TIME_MAX       30.00
#define INITIAL_TEMP   0.05
#define FINAL_TEMP     1.0e-5
#define J_STEP         0.35
#define TIME_STEP      1.00
#define PROGRESS_LINES 20

typedef struct {
    int a;
    int b;
    int orbit;
} Bond;

typedef struct {
    int L;
    int nsites;
    int nbonds;
    int norbits;
    int source;
    int target;
    Bond *bonds;
    double *hamiltonian;
    double *eigenvectors;
    double *eigenvalues;
} TransferProblem;

static void usage(const char *program);
static int parse_int(const char *text, int minimum, int *value);
static int site_index(int x, int y, int L);
static int reflected_site(int site, int L);
static void ordered_pair(int *a, int *b);
static int find_bond(const Bond *bonds, int nbonds, int a, int b);
static int initialize_problem(TransferProblem *problem, int L);
static void free_problem(TransferProblem *problem);
static void build_hamiltonian(TransferProblem *problem,
                              const double *couplings);
static double transfer_fidelity(TransferProblem *problem,
                                const double *couplings,
                                double time);
static double clamp(double value, double lower, double upper);
static int write_result(const char *filename,
                        const TransferProblem *problem,
                        const double *couplings,
                        double time,
                        double fidelity,
                        int seed,
                        int iterations);

int main(int argc, char **argv)
{
    int L = DEFAULT_L;
    int iterations = DEFAULT_ITERATIONS;
    int seed = DEFAULT_SEED;
    const char *output_file = DEFAULT_OUTPUT;

    if (argc > 5 ||
        (argc > 1 && !parse_int(argv[1], 2, &L)) ||
        (argc > 2 && !parse_int(argv[2], 1, &iterations)) ||
        (argc > 3 && !parse_int(argv[3], 0, &seed))) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (argc > 4) {
        output_file = argv[4];
    }
    if (L % 2 != 0) {
        fprintf(stderr,
                "L must be even so left-right reflection maps the lattice "
                "onto itself.\n");
        return EXIT_FAILURE;
    }

    TransferProblem problem = {0};
    if (!initialize_problem(&problem, L)) {
        free_problem(&problem);
        return EXIT_FAILURE;
    }

    double *current = malloc((size_t) problem.norbits * sizeof *current);
    double *best = malloc((size_t) problem.norbits * sizeof *best);
    if (current == NULL || best == NULL) {
        fprintf(stderr, "Could not allocate coupling arrays.\n");
        free(current);
        free(best);
        free_problem(&problem);
        return EXIT_FAILURE;
    }

    set_seed(seed);
    for (int i = 0; i < problem.norbits; ++i) {
        current[i] = INITIAL_J;
    }

    double current_time = clamp(INITIAL_TIME, TIME_MIN, TIME_MAX);
    double current_fidelity =
        transfer_fidelity(&problem, current, current_time);
    double best_time = current_time;
    double best_fidelity = current_fidelity;
    memcpy(best, current, (size_t) problem.norbits * sizeof *best);

    printf("Honeycomb: L=%d, sites=%d, bonds=%d, coupling variables=%d\n",
           L, problem.nsites, problem.nbonds, problem.norbits);
    printf("Transfer: (%d,%d) -> (%d,%d), seed=%d, iterations=%d\n",
           problem.source % L, problem.source / L,
           problem.target % L, problem.target / L, seed, iterations);
    printf("All bond magnitudes remain in [%.2f, %.2f].\n", J_MIN, J_MAX);
    printf("initial: fidelity=%.12f at time=%.8f\n",
           current_fidelity, current_time);

    int progress_interval = iterations / PROGRESS_LINES;
    if (progress_interval < 1) {
        progress_interval = 1;
    }

    for (int step = 0; step < iterations; ++step) {
        double fraction = iterations == 1
                              ? 1.0
                              : (double) step / (double) (iterations - 1);
        double temperature =
            INITIAL_TEMP * pow(FINAL_TEMP / INITIAL_TEMP, fraction);
        double scale = sqrt(temperature / INITIAL_TEMP);

        int change_time = randunif(0.0, 1.0) < 0.15;
        int changed_orbit = -1;
        double previous_value = 0.0;
        double proposed_time = current_time;

        if (change_time) {
            proposed_time = clamp(current_time +
                                      randunif(-TIME_STEP, TIME_STEP) *
                                          (0.15 + 0.85 * scale),
                                  TIME_MIN, TIME_MAX);
        } else {
            changed_orbit =
                (int) randunif(0.0, (double) problem.norbits);
            if (changed_orbit == problem.norbits) {
                --changed_orbit;
            }
            previous_value = current[changed_orbit];
            current[changed_orbit] =
                clamp(previous_value +
                          randunif(-J_STEP, J_STEP) *
                              (0.15 + 0.85 * scale),
                      J_MIN, J_MAX);
        }

        double proposed_fidelity =
            transfer_fidelity(&problem, current, proposed_time);
        double improvement = proposed_fidelity - current_fidelity;
        int accept = improvement >= 0.0;
        if (!accept) {
            double probability = exp(improvement / temperature);
            accept = randunif(0.0, 1.0) < probability;
        }

        if (accept) {
            current_time = proposed_time;
            current_fidelity = proposed_fidelity;
            if (current_fidelity > best_fidelity) {
                best_fidelity = current_fidelity;
                best_time = current_time;
                memcpy(best, current,
                       (size_t) problem.norbits * sizeof *best);
            }
        } else if (changed_orbit >= 0) {
            current[changed_orbit] = previous_value;
        }

        if ((step + 1) % progress_interval == 0 ||
            step + 1 == iterations) {
            printf("%7d/%d  T=%9.3g  current=%.10f  best=%.10f"
                   "  best_time=%.6f\n",
                   step + 1, iterations, temperature, current_fidelity,
                   best_fidelity, best_time);
        }
    }

    if (!write_result(output_file, &problem, best, best_time,
                      best_fidelity, seed, iterations)) {
        free(current);
        free(best);
        free_problem(&problem);
        return EXIT_FAILURE;
    }

    printf("Wrote best solution to %s\n", output_file);
    printf("best fidelity = %.15f at time %.12f\n",
           best_fidelity, best_time);

    free(current);
    free(best);
    free_problem(&problem);
    return EXIT_SUCCESS;
}

static void usage(const char *program)
{
    fprintf(stderr,
            "Usage: %s [even_L] [iterations] [seed] [output_file]\n",
            program);
}

static int parse_int(const char *text, int minimum, int *value)
{
    char *end = NULL;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        parsed < minimum || parsed > INT_MAX) {
        return 0;
    }
    *value = (int) parsed;
    return 1;
}

static int site_index(int x, int y, int L)
{
    return x + L * y;
}

static int reflected_site(int site, int L)
{
    int x = site % L;
    int y = site / L;
    return site_index(L - 1 - x, y, L);
}

static void ordered_pair(int *a, int *b)
{
    if (*a > *b) {
        int temporary = *a;
        *a = *b;
        *b = temporary;
    }
}

static int find_bond(const Bond *bonds, int nbonds, int a, int b)
{
    ordered_pair(&a, &b);
    for (int i = 0; i < nbonds; ++i) {
        if (bonds[i].a == a && bonds[i].b == b) {
            return i;
        }
    }
    return -1;
}

static int initialize_problem(TransferProblem *problem, int L)
{
    int nsites = L * L;
    int maximum_bonds = 2 * nsites;
    Bond *bonds = calloc((size_t) maximum_bonds, sizeof *bonds);
    if (bonds == NULL) {
        fprintf(stderr, "Could not allocate the bond list.\n");
        return 0;
    }

    int nbonds = 0;
    for (int y = 0; y < L; ++y) {
        for (int x = 0; x < L; ++x) {
            int site = site_index(x, y, L);
            if (y + 1 < L) {
                bonds[nbonds++] =
                    (Bond) {site, site_index(x, y + 1, L), -1};
            }
            if (x + 1 < L && (x + y) % 2 == 0) {
                bonds[nbonds++] =
                    (Bond) {site, site_index(x + 1, y, L), -1};
            }
        }
    }

    int norbits = 0;
    for (int i = 0; i < nbonds; ++i) {
        if (bonds[i].orbit >= 0) {
            continue;
        }
        int reflected_a = reflected_site(bonds[i].a, L);
        int reflected_b = reflected_site(bonds[i].b, L);
        int reflected =
            find_bond(bonds, nbonds, reflected_a, reflected_b);
        if (reflected < 0) {
            fprintf(stderr,
                    "Internal error: reflection is not a lattice "
                    "automorphism for bond %d.\n", i);
            free(bonds);
            return 0;
        }
        bonds[i].orbit = norbits;
        bonds[reflected].orbit = norbits;
        ++norbits;
    }

    problem->L = L;
    problem->nsites = nsites;
    problem->nbonds = nbonds;
    problem->norbits = norbits;
    problem->source = site_index(0, 0, L);
    problem->target = site_index(L - 1, 0, L);
    problem->bonds = bonds;
    problem->hamiltonian =
        calloc((size_t) nsites * nsites, sizeof *problem->hamiltonian);
    problem->eigenvectors =
        calloc((size_t) nsites * nsites, sizeof *problem->eigenvectors);
    problem->eigenvalues =
        calloc((size_t) nsites, sizeof *problem->eigenvalues);

    if (problem->hamiltonian == NULL ||
        problem->eigenvectors == NULL ||
        problem->eigenvalues == NULL) {
        fprintf(stderr, "Could not allocate eigensystem arrays.\n");
        return 0;
    }
    return 1;
}

static void free_problem(TransferProblem *problem)
{
    free(problem->bonds);
    free(problem->hamiltonian);
    free(problem->eigenvectors);
    free(problem->eigenvalues);
    *problem = (TransferProblem) {0};
}

static void build_hamiltonian(TransferProblem *problem,
                              const double *couplings)
{
    size_t elements = (size_t) problem->nsites * problem->nsites;
    memset(problem->hamiltonian, 0,
           elements * sizeof *problem->hamiltonian);

    for (int i = 0; i < problem->nbonds; ++i) {
        int a = problem->bonds[i].a;
        int b = problem->bonds[i].b;
        double hopping = -couplings[problem->bonds[i].orbit];
        problem->hamiltonian[a + problem->nsites * b] = hopping;
        problem->hamiltonian[b + problem->nsites * a] = hopping;
    }
}

static double transfer_fidelity(TransferProblem *problem,
                                const double *couplings,
                                double time)
{
    build_hamiltonian(problem, couplings);
    dsydiagonalize(problem->hamiltonian, problem->eigenvectors,
                   problem->eigenvalues, problem->nsites);

    double real_part = 0.0;
    double imaginary_part = 0.0;
    for (int k = 0; k < problem->nsites; ++k) {
        double overlap =
            problem->eigenvectors[problem->source +
                                  problem->nsites * k] *
            problem->eigenvectors[problem->target +
                                  problem->nsites * k];
        double phase = problem->eigenvalues[k] * time;
        real_part += overlap * cos(phase);
        imaginary_part -= overlap * sin(phase);
    }
    return real_part * real_part + imaginary_part * imaginary_part;
}

static double clamp(double value, double lower, double upper)
{
    if (value < lower) {
        return lower;
    }
    if (value > upper) {
        return upper;
    }
    return value;
}

static int write_result(const char *filename,
                        const TransferProblem *problem,
                        const double *couplings,
                        double time,
                        double fidelity,
                        int seed,
                        int iterations)
{
    FILE *output = fopen(filename, "w");
    if (output == NULL) {
        perror(filename);
        return 0;
    }

    fprintf(output, "# simulated-annealing honeycomb couplings\n");
    fprintf(output, "# L %d\n", problem->L);
    fprintf(output, "# source 0 0\n");
    fprintf(output, "# target %d 0\n", problem->L - 1);
    fprintf(output, "# time %.17g\n", time);
    fprintf(output, "# fidelity %.17g\n", fidelity);
    fprintf(output, "# seed %d\n", seed);
    fprintf(output, "# iterations %d\n", iterations);
    fprintf(output, "# columns: x1 y1 x2 y2 coupling\n");

    for (int i = 0; i < problem->nbonds; ++i) {
        int a = problem->bonds[i].a;
        int b = problem->bonds[i].b;
        fprintf(output, "%d %d %d %d %.17g\n",
                a % problem->L, a / problem->L,
                b % problem->L, b / problem->L,
                couplings[problem->bonds[i].orbit]);
    }

    if (fclose(output) != 0) {
        perror(filename);
        return 0;
    }
    return 1;
}
