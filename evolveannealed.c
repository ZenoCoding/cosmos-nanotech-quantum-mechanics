/*
 * Evolve a state using couplings produced by annealhex.c.
 *
 * Compile:
 *
 *     gcc -O3 -std=c11 -D_POSIX_C_SOURCE=200809L \
 *         evolveannealed.c utilities.c -llapack -lblas -lm \
 *         -o evolveannealed
 *
 * Run:
 *
 *     ./evolveannealed [coupling_file] [evolution_file] [frames]
 *
 * Example:
 *
 *     ./evolveannealed annealed_couplings.dat annealed_evolution.txt 240
 */

#include "utilities.h"

#include <errno.h>
#include <limits.h>
#include <string.h>

#define DEFAULT_INPUT  "annealed_couplings.dat"
#define DEFAULT_OUTPUT "annealed_evolution.txt"
#define DEFAULT_FRAMES 240
#define LINE_LENGTH    512

typedef struct {
    int x1;
    int y1;
    int x2;
    int y2;
    double coupling;
} Coupling;

static void usage(const char *program);
static int parse_positive_int(const char *text, int *value);
static int read_metadata(const char *filename, int *L, double *transfer_time,
                         double *predicted_fidelity);
static int read_couplings(const char *filename, int L,
                          Coupling **couplings, int *count);
static int expected_bond_count(int L);
static int site_index(int x, int y, int L);
static int are_honeycomb_neighbors(const Coupling *bond);
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

    int L = 0;
    double transfer_time = 0.0;
    double predicted_fidelity = -1.0;
    if (!read_metadata(input_file, &L, &transfer_time,
                       &predicted_fidelity)) {
        return EXIT_FAILURE;
    }

    Coupling *couplings = NULL;
    int coupling_count = 0;
    if (!read_couplings(input_file, L, &couplings, &coupling_count)) {
        return EXIT_FAILURE;
    }

    int nsites = L * L;
    double complex *hamiltonian =
        calloc((size_t) nsites * nsites, sizeof *hamiltonian);
    double complex *eigenvectors =
        calloc((size_t) nsites * nsites, sizeof *eigenvectors);
    double *eigenvalues = calloc((size_t) nsites, sizeof *eigenvalues);
    double complex *psi = calloc((size_t) nsites, sizeof *psi);
    double complex *work = calloc((size_t) nsites, sizeof *work);
    double *probabilities = calloc((size_t) nsites, sizeof *probabilities);

    if (hamiltonian == NULL || eigenvectors == NULL ||
        eigenvalues == NULL || psi == NULL || work == NULL ||
        probabilities == NULL) {
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

    for (int i = 0; i < coupling_count; ++i) {
        int a = site_index(couplings[i].x1, couplings[i].y1, L);
        int b = site_index(couplings[i].x2, couplings[i].y2, L);
        double hopping = -couplings[i].coupling;
        hamiltonian[a + nsites * b] = hopping;
        hamiltonian[b + nsites * a] = hopping;
    }

    zhediagonalize(hamiltonian, eigenvectors, eigenvalues, nsites);
    psi[site_index(0, 0, L)] = 1.0;

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

    double timestep = transfer_time / (frames - 1);
    for (int frame = 0; frame < frames; ++frame) {
        if (frame > 0) {
            zhediagexpmv(eigenvectors, eigenvalues, -I * timestep,
                         psi, nsites, work);
        }
        for (int site = 0; site < nsites; ++site) {
            probabilities[site] = modulus_squared(psi[site]);
        }
        dvfprint(probabilities, nsites, output);
    }

    int target = site_index(L - 1, 0, L);
    double final_fidelity = probabilities[target];
    double norm = 0.0;
    for (int site = 0; site < nsites; ++site) {
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

    printf("Loaded %d sites and %d honeycomb bonds from %s\n",
           nsites, coupling_count, input_file);
    printf("Generated %d frames from t=0 to t=%.12f in %s\n",
           frames, transfer_time, output_file);
    if (predicted_fidelity >= 0.0) {
        printf("Annealer fidelity: %.15f\n", predicted_fidelity);
    }
    printf("Evolved fidelity:  %.15f\n", final_fidelity);
    printf("Final norm:        %.15f\n", norm);
    printf("Animate with:\n");
    printf("  python3 animatedwavehex.py %s --tmax %.12f "
           "--save annealed_transfer.gif\n",
           output_file, transfer_time);

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

static int read_metadata(const char *filename, int *L, double *transfer_time,
                         double *predicted_fidelity)
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
        if (sscanf(line, "# L %d", &integer_value) == 1) {
            *L = integer_value;
        } else if (sscanf(line, "# time %lf", &real_value) == 1) {
            *transfer_time = real_value;
        } else if (sscanf(line, "# fidelity %lf", &real_value) == 1) {
            *predicted_fidelity = real_value;
        }
    }

    if (ferror(input)) {
        perror(filename);
        fclose(input);
        return 0;
    }
    fclose(input);

    if (*L < 2 || *transfer_time <= 0.0) {
        fprintf(stderr,
                "%s is missing a valid '# L' or '# time' entry.\n",
                filename);
        return 0;
    }
    if (*L > 46340) {
        fprintf(stderr, "L in %s is too large.\n", filename);
        return 0;
    }
    return 1;
}

static int read_couplings(const char *filename, int L,
                          Coupling **couplings, int *count)
{
    FILE *input = fopen(filename, "r");
    if (input == NULL) {
        perror(filename);
        return 0;
    }

    int expected = expected_bond_count(L);
    Coupling *result = calloc((size_t) expected, sizeof *result);
    if (result == NULL) {
        fprintf(stderr, "Could not allocate the coupling list.\n");
        fclose(input);
        return 0;
    }

    char line[LINE_LENGTH];
    int used = 0;
    while (fgets(line, sizeof line, input) != NULL) {
        if (line[0] == '#') {
            continue;
        }

        Coupling bond;
        char extra;
        int fields = sscanf(line, "%d %d %d %d %lf %c",
                            &bond.x1, &bond.y1, &bond.x2, &bond.y2,
                            &bond.coupling, &extra);
        if (fields != 5) {
            fprintf(stderr, "Malformed coupling line in %s: %s",
                    filename, line);
            free(result);
            fclose(input);
            return 0;
        }
        if (bond.x1 < 0 || bond.x1 >= L ||
            bond.x2 < 0 || bond.x2 >= L ||
            bond.y1 < 0 || bond.y1 >= L ||
            bond.y2 < 0 || bond.y2 >= L ||
            bond.coupling <= 0.0 ||
            !are_honeycomb_neighbors(&bond)) {
            fprintf(stderr, "Invalid honeycomb bond in %s: %s",
                    filename, line);
            free(result);
            fclose(input);
            return 0;
        }
        if (used >= expected) {
            fprintf(stderr, "%s contains too many bonds.\n", filename);
            free(result);
            fclose(input);
            return 0;
        }
        result[used++] = bond;
    }

    if (ferror(input)) {
        perror(filename);
        free(result);
        fclose(input);
        return 0;
    }
    fclose(input);

    if (used != expected) {
        fprintf(stderr, "%s contains %d bonds; expected %d for L=%d.\n",
                filename, used, expected, L);
        free(result);
        return 0;
    }

    *couplings = result;
    *count = used;
    return 1;
}

static int expected_bond_count(int L)
{
    int vertical = L * (L - 1);
    int horizontal = L * (L - 1) / 2;
    return vertical + horizontal;
}

static int site_index(int x, int y, int L)
{
    return x + L * y;
}

static int are_honeycomb_neighbors(const Coupling *bond)
{
    int dx = bond->x2 - bond->x1;
    int dy = bond->y2 - bond->y1;

    if (dx == 0 && dy == 1) {
        return 1;
    }
    return dx == 1 && dy == 0 && (bond->x1 + bond->y1) % 2 == 0;
}

static double modulus_squared(double complex value)
{
    return creal(value * conj(value));
}
