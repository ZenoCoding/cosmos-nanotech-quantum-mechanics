/*
 * Simulated annealing for end-to-end state transfer through the 1D barrier
 * used in timeevolve.c.
 *
 * The chain has 20 sites, a fixed on-site energy E=-4, and a fixed barrier
 * on site m/2.  All 19 nearest-neighbor hopping magnitudes and the arrival
 * time are optimized.  The hopping bound prevents the optimizer from making
 * the barrier negligible merely by taking arbitrarily large couplings.
 *
 * Compile:
 *
 *     gcc -O3 -std=c11 -D_POSIX_C_SOURCE=200809L \
 *         annealbarrier.c utilities.c -llapack -lblas -lm \
 *         -o annealbarrier
 *
 * Run:
 *
 *     ./annealbarrier [barrier_height] [iterations] [seed] [output_file]
 *
 * Example matching timeevolve.c:
 *
 *     ./annealbarrier 2.0 50000 123 barrier_couplings.dat
 */

#include "utilities.h"

#include <errno.h>
#include <limits.h>
#include <string.h>

#define SITE_COUNT         20
#define ONSITE_ENERGY      -4.0
#define DEFAULT_BARRIER    2.0
#define DEFAULT_ITERATIONS 50000
#define DEFAULT_SEED       123
#define DEFAULT_OUTPUT     "barrier_couplings.dat"

#define COUPLING_MIN  0.10
#define COUPLING_MAX  2.00
#define TIME_MIN      0.05
#define TIME_MAX      25.132741228718345 /* 8*pi, as in timeevolve.c */
#define INITIAL_TEMP  0.03
#define FINAL_TEMP    1.0e-7
#define COUPLING_STEP 0.25
#define TIME_STEP     0.50
#define PI            3.14159265358979323846

typedef struct {
    double *hamiltonian;
    double *eigenvectors;
    double *eigenvalues;
} Workspace;

static void usage(const char *program);
static int parse_int(const char *text, int minimum, int *value);
static int parse_double(const char *text, double minimum, double *value);
static double clamp(double value, double lower, double upper);
static void build_hamiltonian(Workspace *workspace,
                              const double *couplings,
                              double barrier);
static double fidelity(Workspace *workspace, const double *couplings,
                       double barrier, double time);
static double spectral_fidelity(const Workspace *workspace, double time);
static double best_uniform_fidelity(Workspace *workspace, double barrier,
                                    double *best_time);
static int write_result(const char *filename, const double *couplings,
                        double barrier, double time, double result_fidelity,
                        int seed, int iterations);

int main(int argc, char **argv)
{
    double barrier = DEFAULT_BARRIER;
    int iterations = DEFAULT_ITERATIONS;
    int seed = DEFAULT_SEED;
    const char *output_file = DEFAULT_OUTPUT;

    if (argc > 5 ||
        (argc > 1 && !parse_double(argv[1], 0.0, &barrier)) ||
        (argc > 2 && !parse_int(argv[2], 1, &iterations)) ||
        (argc > 3 && !parse_int(argv[3], 0, &seed))) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (argc > 4) {
        output_file = argv[4];
    }

    Workspace workspace = {
        .hamiltonian =
            calloc(SITE_COUNT * SITE_COUNT, sizeof(double)),
        .eigenvectors =
            calloc(SITE_COUNT * SITE_COUNT, sizeof(double)),
        .eigenvalues = calloc(SITE_COUNT, sizeof(double))
    };
    if (workspace.hamiltonian == NULL ||
        workspace.eigenvectors == NULL ||
        workspace.eigenvalues == NULL) {
        fprintf(stderr, "Could not allocate eigensystem arrays.\n");
        free(workspace.hamiltonian);
        free(workspace.eigenvectors);
        free(workspace.eigenvalues);
        return EXIT_FAILURE;
    }

    double uniform_time = 0.0;
    double uniform_fidelity =
        best_uniform_fidelity(&workspace, barrier, &uniform_time);

    double current[SITE_COUNT - 1];
    double best[SITE_COUNT - 1];

    /*
     * Start at the exactly solvable mirror-inverting chain, scaled so its
     * largest hopping equals COUPLING_MAX:
     *
     *     J_n = scale * sqrt(n * (N - n)).
     */
    double largest = sqrt((SITE_COUNT / 2.0) *
                          (SITE_COUNT - SITE_COUNT / 2.0));
    double scale = COUPLING_MAX / largest;
    for (int n = 1; n < SITE_COUNT; ++n) {
        current[n - 1] = scale * sqrt(n * (SITE_COUNT - n));
    }
    double current_time =
        clamp(PI / (2.0 * scale), TIME_MIN, TIME_MAX);
    double current_fidelity =
        fidelity(&workspace, current, barrier, current_time);
    double best_time = current_time;
    double best_fidelity = current_fidelity;
    memcpy(best, current, sizeof best);

    set_seed(seed);
    printf("20-site chain, barrier site=%d, barrier height=%.8g\n",
           SITE_COUNT / 2, barrier);
    printf("Coupling bounds [%.2f, %.2f], time bounds [%.2f, %.2f]\n",
           COUPLING_MIN, COUPLING_MAX, TIME_MIN, TIME_MAX);
    printf("uniform chain:   fidelity=%.12f at time=%.6f\n",
           uniform_fidelity, uniform_time);
    printf("initial profile: fidelity=%.12f at time=%.6f\n",
           current_fidelity, current_time);

    int progress_interval = iterations / 10;
    if (progress_interval < 1) {
        progress_interval = 1;
    }

    for (int step = 0; step < iterations; ++step) {
        double fraction = iterations == 1
                              ? 1.0
                              : (double) step / (iterations - 1);
        double temperature =
            INITIAL_TEMP * pow(FINAL_TEMP / INITIAL_TEMP, fraction);
        double step_scale = 0.05 + 0.95 * sqrt(temperature / INITIAL_TEMP);

        int change_time = randunif(0.0, 1.0) < 0.12;
        int changed_coupling = -1;
        double old_value = 0.0;
        double proposed_time = current_time;

        if (change_time) {
            proposed_time =
                clamp(current_time +
                          randunif(-TIME_STEP, TIME_STEP) * step_scale,
                      TIME_MIN, TIME_MAX);
        } else {
            changed_coupling =
                (int) randunif(0.0, SITE_COUNT - 1.0);
            if (changed_coupling == SITE_COUNT - 1) {
                --changed_coupling;
            }
            old_value = current[changed_coupling];
            current[changed_coupling] =
                clamp(old_value +
                          randunif(-COUPLING_STEP, COUPLING_STEP) *
                              step_scale,
                      COUPLING_MIN, COUPLING_MAX);
        }

        double proposed_fidelity =
            fidelity(&workspace, current, barrier, proposed_time);
        double improvement = proposed_fidelity - current_fidelity;
        int accept = improvement >= 0.0;
        if (!accept) {
            accept = randunif(0.0, 1.0) <
                     exp(improvement / temperature);
        }

        if (accept) {
            current_time = proposed_time;
            current_fidelity = proposed_fidelity;
            if (current_fidelity > best_fidelity) {
                best_fidelity = current_fidelity;
                best_time = current_time;
                memcpy(best, current, sizeof best);
            }
        } else if (changed_coupling >= 0) {
            current[changed_coupling] = old_value;
        }

        if ((step + 1) % progress_interval == 0 ||
            step + 1 == iterations) {
            printf("%7d/%d  T=%9.3g  current=%.10f  best=%.10f"
                   "  time=%.6f\n",
                   step + 1, iterations, temperature, current_fidelity,
                   best_fidelity, best_time);
        }
    }

    if (!write_result(output_file, best, barrier, best_time,
                      best_fidelity, seed, iterations)) {
        free(workspace.hamiltonian);
        free(workspace.eigenvectors);
        free(workspace.eigenvalues);
        return EXIT_FAILURE;
    }

    printf("optimized chain: fidelity=%.15f at time=%.9f\n",
           best_fidelity, best_time);
    printf("Wrote %s\n", output_file);

    free(workspace.hamiltonian);
    free(workspace.eigenvectors);
    free(workspace.eigenvalues);
    return EXIT_SUCCESS;
}

static void usage(const char *program)
{
    fprintf(stderr,
            "Usage: %s [barrier>=0] [iterations] [seed] [output_file]\n",
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

static int parse_double(const char *text, double minimum, double *value)
{
    char *end = NULL;
    errno = 0;
    double parsed = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' ||
        !isfinite(parsed) || parsed < minimum) {
        return 0;
    }
    *value = parsed;
    return 1;
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

static void build_hamiltonian(Workspace *workspace,
                              const double *couplings,
                              double barrier)
{
    memset(workspace->hamiltonian, 0,
           SITE_COUNT * SITE_COUNT * sizeof(double));
    for (int site = 0; site < SITE_COUNT; ++site) {
        workspace->hamiltonian[site + SITE_COUNT * site] =
            ONSITE_ENERGY;
    }
    workspace->hamiltonian[SITE_COUNT / 2 +
                           SITE_COUNT * (SITE_COUNT / 2)] += barrier;

    for (int site = 0; site < SITE_COUNT - 1; ++site) {
        double hopping = -couplings[site];
        workspace->hamiltonian[site + SITE_COUNT * (site + 1)] =
            hopping;
        workspace->hamiltonian[site + 1 + SITE_COUNT * site] =
            hopping;
    }
}

static double fidelity(Workspace *workspace, const double *couplings,
                       double barrier, double time)
{
    build_hamiltonian(workspace, couplings, barrier);
    dsydiagonalize(workspace->hamiltonian, workspace->eigenvectors,
                   workspace->eigenvalues, SITE_COUNT);
    return spectral_fidelity(workspace, time);
}

static double spectral_fidelity(const Workspace *workspace, double time)
{
    double real_part = 0.0;
    double imaginary_part = 0.0;
    for (int k = 0; k < SITE_COUNT; ++k) {
        double overlap =
            workspace->eigenvectors[SITE_COUNT * k] *
            workspace->eigenvectors[SITE_COUNT - 1 + SITE_COUNT * k];
        double phase = workspace->eigenvalues[k] * time;
        real_part += overlap * cos(phase);
        imaginary_part -= overlap * sin(phase);
    }
    return real_part * real_part + imaginary_part * imaginary_part;
}

static double best_uniform_fidelity(Workspace *workspace, double barrier,
                                    double *best_time)
{
    double couplings[SITE_COUNT - 1];
    for (int i = 0; i < SITE_COUNT - 1; ++i) {
        couplings[i] = 1.0;
    }
    build_hamiltonian(workspace, couplings, barrier);
    dsydiagonalize(workspace->hamiltonian, workspace->eigenvectors,
                   workspace->eigenvalues, SITE_COUNT);

    double best = 0.0;
    *best_time = TIME_MIN;
    const int samples = 30000;
    for (int sample = 0; sample <= samples; ++sample) {
        double time = TIME_MIN +
                      (TIME_MAX - TIME_MIN) * sample / samples;
        double candidate = spectral_fidelity(workspace, time);
        if (candidate > best) {
            best = candidate;
            *best_time = time;
        }
    }
    return best;
}

static int write_result(const char *filename, const double *couplings,
                        double barrier, double time, double result_fidelity,
                        int seed, int iterations)
{
    FILE *output = fopen(filename, "w");
    if (output == NULL) {
        perror(filename);
        return 0;
    }

    fprintf(output, "# annealed 1D barrier-chain couplings\n");
    fprintf(output, "# sites %d\n", SITE_COUNT);
    fprintf(output, "# onsite_energy %.17g\n", ONSITE_ENERGY);
    fprintf(output, "# barrier_site %d\n", SITE_COUNT / 2);
    fprintf(output, "# barrier_height %.17g\n", barrier);
    fprintf(output, "# time %.17g\n", time);
    fprintf(output, "# fidelity %.17g\n", result_fidelity);
    fprintf(output, "# seed %d\n", seed);
    fprintf(output, "# iterations %d\n", iterations);
    fprintf(output, "# columns: left_site right_site hopping_magnitude\n");
    for (int site = 0; site < SITE_COUNT - 1; ++site) {
        fprintf(output, "%d %d %.17g\n",
                site, site + 1, couplings[site]);
    }

    if (fclose(output) != 0) {
        perror(filename);
        return 0;
    }
    return 1;
}
