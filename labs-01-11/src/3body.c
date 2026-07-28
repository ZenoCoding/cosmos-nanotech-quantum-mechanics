#include <math.h>
#include <stdio.h>

/* Units: AU, solar masses, and years. */
#define PI 3.14159265358979323846
#define N_BODIES 3
#define G (4.0 * PI * PI)
#define DT 0.0005
#define END_TIME 30.0
#ifndef PERTURBATION
#define PERTURBATION 0.005
#endif

enum body { BODY_A, BODY_B, BODY_C };

/* Fill ax and ay with the acceleration of every body from every other body. */
static void calculate_accelerations(const double mass[], const double x[],
                                    const double y[], double ax[], double ay[])
{
    for (int i = 0; i < N_BODIES; i++) {
        ax[i] = 0.0;
        ay[i] = 0.0;
    }

    for (int i = 0; i < N_BODIES; i++) {
        for (int j = i + 1; j < N_BODIES; j++) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            double distance_squared = dx * dx + dy * dy;
            double inverse_distance_cubed =
                1.0 / (distance_squared * sqrt(distance_squared));

            ax[i] += G * mass[j] * dx * inverse_distance_cubed;
            ay[i] += G * mass[j] * dy * inverse_distance_cubed;
            ax[j] -= G * mass[i] * dx * inverse_distance_cubed;
            ay[j] -= G * mass[i] * dy * inverse_distance_cubed;
        }
    }
}

int main(void)
{
    /*
     * A slightly perturbed equal-mass figure-eight. The perturbation breaks the
     * exact repetition while preserving the otherwise bounded starting orbit.
     */
    const double mass[N_BODIES] = {1.0, 1.0, 1.0};
    double x[N_BODIES] = {0.97000436, -0.97000436, 0.0};
    double y[N_BODIES] = {-0.24308753, 0.24308753, 0.0};
    double vx[N_BODIES] = {0.466203685 * 2.0 * PI, 0.466203685 * 2.0 * PI, 0.0};
    double vy[N_BODIES] = {
        0.432365730 * 2.0 * PI,
        (0.432365730 + PERTURBATION) * 2.0 * PI,
        0.0
    };
    double ax[N_BODIES];
    double ay[N_BODIES];
    double next_ax[N_BODIES];
    double next_ay[N_BODIES];

    /* Give body C the opposite momentum, keeping the system barycenter fixed. */
    vx[BODY_C] = -(mass[BODY_A] * vx[BODY_A] + mass[BODY_B] * vx[BODY_B]) /
                 mass[BODY_C];
    vy[BODY_C] = -(mass[BODY_A] * vy[BODY_A] + mass[BODY_B] * vy[BODY_B]) /
                 mass[BODY_C];

    FILE *fileout = fopen("data/3body", "w");
    if (fileout == NULL) {
        perror("data/3body");
        return 1;
    }

    fprintf(fileout, "# t body_a_x body_a_y body_b_x body_b_y body_c_x body_c_y\n");
    calculate_accelerations(mass, x, y, ax, ay);

    for (double t = 0.0; t <= END_TIME; t += DT) {
        fprintf(fileout, "%.6f %.9f %.9f %.9f %.9f %.9f %.9f\n", t,
                x[BODY_A], y[BODY_A], x[BODY_B], y[BODY_B], x[BODY_C], y[BODY_C]);

        /* Velocity-Verlet integration is much more stable than Euler for orbits. */
        for (int i = 0; i < N_BODIES; i++) {
            x[i] += vx[i] * DT + 0.5 * ax[i] * DT * DT;
            y[i] += vy[i] * DT + 0.5 * ay[i] * DT * DT;
        }

        calculate_accelerations(mass, x, y, next_ax, next_ay);
        for (int i = 0; i < N_BODIES; i++) {
            vx[i] += 0.5 * (ax[i] + next_ax[i]) * DT;
            vy[i] += 0.5 * (ay[i] + next_ay[i]) * DT;
            ax[i] = next_ax[i];
            ay[i] = next_ay[i];
        }
    }

    fclose(fileout);
    return 0;
}
