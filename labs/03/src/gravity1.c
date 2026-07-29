#include <stdio.h>
#include <math.h>

int main(){
    double y, v, t, dt, F, m, E, g=-9.8;
    int N;
    FILE * fileout;
    fileout = fopen("labs/03/data/gravity", "w");
    printf("\nEnter mass m:   ");
    scanf("%lf", &m);
    printf("\nEnter starting height y and velocity v:    ");
    scanf("%lf %lf", &y, &v);
    printf("\nEnter time step dt:");
    scanf("%lf", &dt);
    t = 0.;
    E = .5 *  m * v * v + m * g * y;
    fprintf(fileout, "\n   %12.6lf %12.6lf %12.6lf", t, y, E);
    for(;y>=0;t+=dt) {
        y += v * dt;
        F = m*g;
        v += F/m * dt;
        E = .5 *  m * v * v + m * g * y;
        fprintf(fileout, "\n   %12.6lf %12.6lf %12.6lf", t, y, E);
    }
    fclose(fileout);
}
