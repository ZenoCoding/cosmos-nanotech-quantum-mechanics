#include <stdio.h>
#include <math.h>

int main(){
    double x = 0.0, y = 0.0, vx, vy, t, dt, Fx, Fy, m, g=-9.8;
    FILE * fileout;
    fileout = fopen("../data/projmotion", "w");
    printf("\nEnter mass m:   ");
    scanf("%lf", &m);
    printf("\nEnter starting height y and velocity vx vy:    ");
    scanf("%lf %lf %lf", &y, &vx, &vy);
    printf("\nEnter time step dt:");
    scanf("%lf", &dt);
    t = 0.;
    fprintf(fileout, "\n   %12.6lf %12.6lf %12.6lf", t, x, y);
    for(;y>=0;t+=dt) {
        // Position Updates
        x += vx * dt;
        y += vy * dt;

        // Force Updates
        Fx = 0.0;
        Fy = m*g;

        // Velocity Updates
        vx += Fx/m * dt;
        vy += Fy/m * dt;
        fprintf(fileout, "\n   %12.6lf %12.6lf %12.6lf", t, x, y);
    }
    fclose(fileout);
}
