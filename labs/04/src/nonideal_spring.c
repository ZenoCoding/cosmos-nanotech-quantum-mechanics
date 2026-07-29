#include <stdio.h>
#include <math.h>

int main(){
    double x, v = 0.0, t = 0.0, k, k_a = 20, dt = 0.001, F, m;
    int drag = 0;
    int N;
    FILE * fileout;
    fileout = fopen("labs/04/data/mass_spring", "w");
    printf("\nEnter mass m:   ");
    scanf("%lf", &m);
    printf("\nEnter starting position x, and spring constant k: ");
    scanf("%lf %lf", &x, &k);
    fprintf(fileout, "\n   %12.6lf %12.6lf", t, x);
    for(t=0.;t <= 1.0;t+=dt) {
        x += v * dt;
        F = -k*x - k_a * x * x * x;
        v += F/m * dt;
        fprintf(fileout, "\n   %12.6lf %12.6lf", t, x);
    }
    fclose(fileout);
}
