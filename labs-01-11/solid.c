#include <stdio.h>
#include <math.h>
#include <time.h>    // for seeding rand()
#include <stdlib.h>  // for rand(), RAND_MAX
#include <stdbool.h>

typedef struct {
    double m;     // mass
    double F[2];  // force vector [Fx, Fy]
    double x[2];  // position vector [x, y]
    double v[2];  // velocity vector [vx, vy]
    double a[2];  // acceleration vector [ax, ay]
} Particle;

int num_atoms = 200;
double dt = 1e-8;
int steps = 1e4;

double box_size = 6.0;
double mass = 1.0;
double v_init = 5e2;
double epsilon = 1e4;
double sigma = 0.5;
double max_force = 1e12;
double Temperature_change = 0.99;

float sign(double a) {
    if (a > 0) return 1.0;
    return -1.0;
}

// Force on a from b
void Leonard_Jones_Force(Particle a, Particle b, double *fx, double *fy) {
    double dx = a.x[0] - b.x[0];
    double dy = a.x[1] - b.x[1];
    double r2 = dx * dx + dy * dy;
    double r = sqrt(r2);

    // Lennard-Jones force magnitude
    double sr2 = (sigma * sigma) / r2;
    double sr6 = sr2 * sr2 * sr2;
    double sr12 = sr6 * sr6;
    double F_mag = 24 * epsilon * (2 * sr12 - sr6) / r2;

    double Fx = F_mag * (dx / r);
    double Fy = F_mag * (dy / r);
    if (isnan(F_mag) || F_mag > max_force) {
        Fx = max_force * sign(dx);
        Fy = max_force * sign(dy);
    }
    // Unit vector components
    *fx = Fx;
    *fy = Fy;
}

// Compute Lennard-Jones forces in 2D
void compute_forces(Particle atoms[], int num_atoms) {
    // Reset all forces
    for (int i = 0; i < num_atoms; ++i) {
        atoms[i].F[0] = 0.0;
        atoms[i].F[1] = 0.0;
    }

    // Compute pairwise forces
    double fx = 0.0;
    double fy = 0.0;
    for (int i = 0; i < num_atoms - 1; ++i) {
        for (int j = i + 1; j < num_atoms; ++j) {
            if (i == j)
                continue;
            Leonard_Jones_Force(atoms[i], atoms[j], &fx, &fy);
            atoms[i].F[0] += fx;
            atoms[i].F[1] += fy;
            atoms[j].F[0] -= fx;
            atoms[j].F[1] -= fy;
        }
    }
}

// Bounce atoms off walls of a square box
void apply_wall_bounce(Particle *p, double box_size) {
    for (int d = 0; d < 2; ++d) {
        double half_box = box_size / 2.0;
        if (p->x[d] >= half_box) {
            p->x[d] = half_box;
            p->v[d] *= -1;
        } else if (p->x[d] <= -half_box) {
            p->x[d] = -half_box;
            p->v[d] *= -1;
        }
    }
}

// Run simulation and save 2D positions to CSV
void simulate_particle_dynamics(Particle atoms[], int num_atoms, double dt, int steps, const char *filename) {
    FILE *out = fopen(filename, "w");
    if (!out) {
        perror("File open error");
        return;
    }

    // Write CSV header
    fprintf(out, "step");
    for (int i = 0; i < num_atoms; ++i) {
        fprintf(out, ",p%d_x,p%d_y", i, i);
    }
    fprintf(out, "\n");

    compute_forces(atoms, num_atoms);
    for (int i = 0; i < num_atoms; ++i) {
        atoms[i].a[0] = atoms[i].F[0] / atoms[i].m;
        atoms[i].a[1] = atoms[i].F[1] / atoms[i].m;
    }

    for (int step = 0; step < steps; ++step) {
        compute_forces(atoms, num_atoms);

        fprintf(out, "%d", step);
        for (int i = 0; i < num_atoms; ++i) {
            atoms[i].a[0] = atoms[i].F[0] / atoms[i].m;
            atoms[i].a[1] = atoms[i].F[1] / atoms[i].m;

            atoms[i].v[0] += atoms[i].a[0] * dt;
            atoms[i].v[1] += atoms[i].a[1] * dt;
            int N = 0.0009 / dt;
            if (i % 1 == 0) {
                atoms[i].v[0] *= Temperature_change;
                atoms[i].v[1] *= Temperature_change;
            }

            atoms[i].x[0] += atoms[i].v[0] * dt;
            atoms[i].x[1] += atoms[i].v[1] * dt;

            apply_wall_bounce(&atoms[i], box_size);  // bounce after position update

            fprintf(out, ",%lf,%lf", atoms[i].x[0], atoms[i].x[1]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
}

// Initialize atoms with random 2D positions
void init_atoms(Particle atoms[], int num_atoms) {
    srand(time(NULL));

    for (int i = 0; i < num_atoms; ++i) {
        // Generate random position in box
        double x = box_size * ((double)rand() / RAND_MAX - 0.5);
        double y = box_size * ((double)rand() / RAND_MAX - 0.5);
        double vx = v_init * ((double)rand() / RAND_MAX - 0.5);
        double vy = v_init * ((double)rand() / RAND_MAX - 0.5);
        atoms[i].m = mass;
        atoms[i].x[0] = x;
        atoms[i].x[1] = y;
        atoms[i].v[0] = vx, atoms[i].v[1] = vy;
        atoms[i].F[0] = atoms[i].F[1] = 0.0;
        atoms[i].a[0] = atoms[i].a[1] = 0.0;
    }
}

int main() {
    Particle atoms[num_atoms];

    init_atoms(atoms, num_atoms);


    simulate_particle_dynamics(atoms, num_atoms, dt, steps, "positions.csv");
    return 0;
}
