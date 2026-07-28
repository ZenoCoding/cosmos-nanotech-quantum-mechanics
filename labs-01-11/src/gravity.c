#include <stdio.h>
#include <math.h>

int main(){
    double y, v, t, dt, F, m, g=-9.8;
    int j, N;
    FILE * fileout;
    fileout = fopen("data/gravity", "w");
    printf("\nEnter mass m:   ");
    scanf("%lf", &m);
    printf("\nEnter starting height y and velocity v:    ");
    scanf("%lf %lf", &y, &v);
    printf("\nEnter time step dt and number of steps N:    ");
    scanf("%lf %i", &dt, &N);
    t = 0.;
    fprintf(fileout, "\n   %12.6lf %12.6lf", t, y);
    for (j = 0; j < N; j++){
        t += dt;
        y += v * dt;
        F = m*g;
        v += F/m * dt;
        fprintf(fileout, "\n   %12.6lf %12.6lf", t, y);
    }
    fclose(fileout);
}
