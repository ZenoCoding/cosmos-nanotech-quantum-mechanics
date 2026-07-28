#include <stdio.h>
#include <math.h>

int main(){
    double x, v, t, k, b, dt, F, m, E, g=-9.8;
    int drag = 0;
    int N;
    FILE * fileout;
    fileout = fopen("labs/01-11/data/mass_spring", "w");
    printf("\nEnter mass m:   ");
    scanf("%lf", &m);
    printf("\nEnter starting position x and velocity v:    ");
    scanf("%lf %lf", &x, &v);
    printf("\nEnter spring constant k and drag constant b:    ");
    scanf("%lf %lf", &k, &b);
    printf("\nEnter time step dt:");
    scanf("%lf", &dt);
    t = 0.;
    E = .5 *  m * v * v + .5*k*x;
    fprintf(fileout, "\n   %12.6lf %12.6lf %12.6lf", t, x, E);
    for(;t <= 5.0;t+=dt) {
        x += v * dt;
        F = -k*x - b*v;
        v += F/m * dt;
        E = .5 * m * v * v + .5 * k * x * x;
        //printf("KE: %g SPE: %g E: %g\n", .5 * m * v * v, .5 * k * x * x, E);
        fprintf(fileout, "\n   %12.6lf %12.6lf %12.6lf", t, x, E);
    }
    fclose(fileout);
}
