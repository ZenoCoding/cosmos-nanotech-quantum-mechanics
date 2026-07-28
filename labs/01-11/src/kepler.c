#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

int main(){
    // Use AU, YRS
    double MASS_EARTH = 3.0e-6, MASS_SUN = 1.0, EARTH_RADIUS = 1.0, DT=0.001, G=4.0*PI*PI;
    double x=EARTH_RADIUS, y=0, r=EARTH_RADIUS, vx=0, vy=1*PI, t=0.0, dt=DT, Fx, Fy, m=MASS_EARTH, b=0, g=-9.8;
    FILE * fileout;
    fileout=fopen("labs/01-11/data/kepler1", "w");

    fprintf(fileout, "%12.6lf %12.6lf %12.6lf", t, x, y);
    for(;t<1.0; t+=dt){
        x += vx * dt;
        y += vy * dt;
        r = sqrt(x*x + y*y);

        Fx = -G*MASS_SUN * m * x / r / r / r;
        Fy = -G*MASS_SUN * m * y / r / r / r;

        vx += Fx / m * dt;
        vy += Fy / m * dt;
        fprintf(fileout, "\n   %12.6lf %12.6lf %12.6lf", t, x, y);
    }
    
}
