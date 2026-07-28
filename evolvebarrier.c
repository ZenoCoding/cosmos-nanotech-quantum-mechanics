/*
 * Evolve a 1D state using couplings produced by annealbarrier.c.
 *
 * Compile:
 *
 *     gcc -O3 -std=c11 -D_POSIX_C_SOURCE=200809L \
 *         evolvebarrier.c utilities.c -llapack -lblas -lm \
 *         -o evolvebarrier
 *
 * Run:
 *
 *     ./evolvebarrier [coupling_file] [evolution_file] [frames]
 *
 * Example:
 *
 *     ./evolvebarrier barrier_couplings.dat \
 *         annealed_barrier_evolution.txt 240
 */

#include "utilities.h"

#include <errno.h>
#include <limits.h>
#include <string.h>

#define DEFAULT_INPUT  "barrier_couplings.dat"
#define DEFAULT_OUTPUT "annealed_barrier_evolution.txt"
#define DEFAULT_FRAMES 240
#define LINE_LENGTH     512

typedef struct {
    int sites;
    int barrier_site;
    double onsite_energy;
    double barrier_height;
    double transfer_time;
    double predicted_fidelity;
} Metadata;

static void usage(const char *program);
static int parse_positive_int(const char *text, int *value);
static int read_metadata(const char *filename, Metadata *metadata);
static int read_couplings(const char *filename, int sites,
                          double *couplings);
static double modulus_squared(double complex value);

int main(int argc, char **argv)
{
    const char *input_file = DEFAULT_INPUT;
    const char *output_file = DEFAULT_OUTPUT;
    int frames = DEFAULT_FRAMES;

    if (argc > 4 ||
        (argc > 3 && !parse_positive_int(argv[3], &frames))) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (argc > 1) {
        input_file = argv[1];
    }
    if (argc > 2) {
        output_file = argv[2];
    }
    if (frames < 2) {
        fprintf(stderr, "frames must be at least 2.\n");
        return EXIT_FAILURE;
    }

    Metadata metadata = {
        .sites = 0,
        .barrier_site = -1,
        .onsite_energy = 0.0,
        .barrier_height = -1.0,
        .transfer_time = 0.0,
        .predicted_fidelity = -1.0
    };
    if (!read_metadata(input_file, &metadata)) {
        return EXIT_FAILURE;
    }

    int n = metadata.sites;
    double *couplings = calloc((size_t) (n - 1), sizeof *couplings);
    double complex *hamiltonian =
        calloc((size_t) n * n, sizeof *hamiltonian);
    double complex *eigenvectors =
        calloc((size_t) n * n, sizeof *eigenvectors);
    double *eigenvalues = calloc((size_t) n, sizeof *eigenvalues);
    double complex *psi = calloc((size_t) n, sizeof *psi);
    double complex *work = calloc((size_t) n, sizeof *work);
    double *probabilities = calloc((size_t) n, sizeof *probabilities);

    if (couplings == NULL || hamiltonian == NULL ||
        eigenvectors == NULL || eigenvalues == NULL || psi == NULL ||
        work == NULL || probabilities == NULL) {
        fprintf(stderr, "Could not allocate evolution arrays.\n");
        free(couplings);
        free(hamiltonian);
        free(eigenvectors);
        free(eigenvalues);
        free(psi);
        free(work);
        free(probabilities);
        return EXIT_FAILURE;
    }
    if (!read_couplings(input_file, n, couplings)) {
        free(couplings);
        free(hamiltonian);
        free(eigenvectors);
        free(eigenvalues);
        free(psi);
        free(work);
        free(probabilities);
        return EXIT_FAILURE;
    }

    for (int site = 0; site < n; ++site) {
        hamiltonian[site + n * site] = metadata.onsite_energy;
    }
    hamiltonian[metadata.barrier_site +
                n * metadata.barrier_site] += metadata.barrier_height;
    for (int site = 0; site < n - 1; ++site) {
        double hopping = -couplings[site];
        hamiltonian[site + n * (site + 1)] = hopping;
        hamiltonian[site + 1 + n * site] = hopping;
    }

    zhediagonalize(hamiltonian, eigenvectors, eigenvalues, n);
    psi[0] = 1.0;

    FILE *output = fopen(output_file, "w");
    if (output == NULL) {
        perror(output_file);
        free(couplings);
        free(hamiltonian);
        free(eigenvectors);
        free(eigenvalues);
        free(psi);
        free(work);
        free(probabilities);
        return EXIT_FAILURE;
    }

    double timestep = metadata.transfer_time / (frames - 1);
    for (int frame = 0; frame < frames; ++frame) {
        if (frame > 0) {
            zhediagexpmv(eigenvectors, eigenvalues, -I * timestep,
                         psi, n, work);
        }
        for (int site = 0; site < n; ++site) {
            probabilities[site] = modulus_squared(psi[site]);
        }
        dvfprint(probabilities, n, output);
    }

    double evolved_fidelity = probabilities[n - 1];
    double norm = 0.0;
    for (int site = 0; site < n; ++site) {
        norm += probabilities[site];
    }

    if (fclose(output) != 0) {
        perror(output_file);
        free(couplings);
        free(hamiltonian);
        free(eigenvectors);
        free(eigenvalues);
        free(psi);
        free(work);
        free(probabilities);
        return EXIT_FAILURE;
    }

    printf("Loaded %d sites and %d hoppings from %s\n",
           n, n - 1, input_file);
    printf("Barrier: site=%d, height=%.8g\n",
           metadata.barrier_site, metadata.barrier_height);
    printf("Generated %d frames from t=0 to t=%.12f in %s\n",
           frames, metadata.transfer_time, output_file);
    if (metadata.predicted_fidelity >= 0.0) {
        printf("Annealer fidelity: %.15f\n",
               metadata.predicted_fidelity);
    }
    printf("Evolved fidelity:  %.15f\n", evolved_fidelity);
    printf("Final norm:        %.15f\n", norm);
    printf("Animate with:\n");
    printf("  python3 animatedwave.py %s --tmax %.12f "
           "--save annealed_barrier.gif\n",
           output_file, metadata.transfer_time);

    free(couplings);
    free(hamiltonian);
    free(eigenvectors);
    free(eigenvalues);
    free(psi);
    free(work);
    free(probabilities);
    return EXIT_SUCCESS;
}

static void usage(const char *program)
{
    fprintf(stderr,
            "Usage: %s [coupling_file] [evolution_file] [frames]\n",
            program);
}

static int parse_positive_int(const char *text, int *value)
{
    char *end = NULL;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        parsed < 1 || parsed > INT_MAX) {
        return 0;
    }
    *value = (int) parsed;
    return 1;
}

static int read_metadata(const char *filename, Metadata *metadata)
{
    FILE *input = fopen(filename, "r");
    if (input == NULL) {
        perror(filename);
        return 0;
    }

    char line[LINE_LENGTH];
    while (fgets(line, sizeof line, input) != NULL) {
        int integer_value;
        double real_value;
        if (sscanf(line, "# sites %d", &integer_value) == 1) {
            metadata->sites = integer_value;
        } else if (sscanf(line, "# onsite_energy %lf",
                          &real_value) == 1) {
            metadata->onsite_energy = real_value;
        } else if (sscanf(line, "# barrier_site %d",
                          &integer_value) == 1) {
            metadata->barrier_site = integer_value;
        } else if (sscanf(line, "# barrier_height %lf",
                          &real_value) == 1) {
            metadata->barrier_height = real_value;
        } else if (sscanf(line, "# time %lf", &real_value) == 1) {
            metadata->transfer_time = real_value;
        } else if (sscanf(line, "# fidelity %lf", &real_value) == 1) {
            metadata->predicted_fidelity = real_value;
        }
    }

    if (ferror(input)) {
        perror(filename);
        fclose(input);
        return 0;
    }
    fclose(input);

    if (metadata->sites < 2 ||
        metadata->barrier_site < 0 ||
        metadata->barrier_site >= metadata->sites ||
        metadata->barrier_height < 0.0 ||
        metadata->transfer_time <= 0.0) {
        fprintf(stderr, "%s has missing or invalid metadata.\n",
                filename);
        return 0;
    }
    return 1;
}

static int read_couplings(const char *filename, int sites,
                          double *couplings)
{
    FILE *input = fopen(filename, "r");
    if (input == NULL) {
        perror(filename);
        return 0;
    }

    unsigned char *seen = calloc((size_t) (sites - 1), sizeof *seen);
    if (seen == NULL) {
        fprintf(stderr, "Could not allocate coupling validation data.\n");
        fclose(input);
        return 0;
    }

    char line[LINE_LENGTH];
    int count = 0;
    while (fgets(line, sizeof line, input) != NULL) {
        if (line[0] == '#') {
            continue;
        }

        int left;
        int right;
        double coupling;
        char extra;
        int fields = sscanf(line, "%d %d %lf %c",
                            &left, &right, &coupling, &extra);
        if (fields != 3 || left < 0 || left >= sites - 1 ||
            right != left + 1 || coupling <= 0.0 || seen[left]) {
            fprintf(stderr, "Invalid coupling line in %s: %s",
                    filename, line);
            free(seen);
            fclose(input);
            return 0;
        }
        couplings[left] = coupling;
        seen[left] = 1;
        ++count;
    }

    if (ferror(input)) {
        perror(filename);
        free(seen);
        fclose(input);
        return 0;
    }
    fclose(input);

    if (count != sites - 1) {
        fprintf(stderr, "%s contains %d hoppings; expected %d.\n",
                filename, count, sites - 1);
        free(seen);
        return 0;
    }

    free(seen);
    return 1;
}

static double modulus_squared(double complex value)
{
    return creal(value * conj(value));
}
