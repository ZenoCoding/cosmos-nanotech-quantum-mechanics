#include <stdio.h>
#include <math.h>

int main(void) {
    FILE * fileout;
    int x, t, Nt, Nx=1000;
    double rho[Nx], newrho[Nx], D, dt, dx, Ddttodx2, norm, realx, realt, BIAS = .0;
    
    fileout=fopen("labs/08/data/diffusion2.txt", "w");

    printf("\nEnter D, dt, dx: ");
    scanf("%lf %lf %lf", &D, &dt, &dx);
    Ddttodx2 = D*dt /(dx*dx);
    printf("\nEnter number of time steps: ");
    scanf("%d", &Nt);

    for (x=0; x<Nx; x++){
        double xfrac = ((double) x) / Nx;
        rho[x] = 0;
         rho[x] = xfrac/500.0;
        //rho[x] = xfrac * xfrac * 3 /
    }

    //rho[Nx/2] = 1.0;

    for(t=0; t<Nt; t++)
    {   
        newrho[0] = rho[0] + Ddttodx2*(rho[1]-rho[0]);
        for(x=1; x<Nx-1; x++){
            newrho[x]=rho[x]+Ddttodx2*(rho[x+1]*(1.0 - BIAS) + rho[x-1]*(1.0 + BIAS)-2.0*rho[x]);
        }
        newrho[Nx-1] = rho[Nx-1] + Ddttodx2*(rho[Nx-2]-rho[Nx-1]);

        norm=0.0;
        for (x=0; x<Nx; x++){
            rho[x]=newrho[x];
            norm += rho[x];
        }

        if(t%100 == 0){
            for(x=0; x<Nx; x++){
                realx = dx*x;
                fprintf(fileout, "%8.4lf ",rho[x]/dx);
            }
        }
        fprintf(fileout, "\n");
        printf("\n checking no particles are lost! t = %8.3lf norm = %8.3lf  ", t*dt, norm);
    }
    
    fclose(fileout);
    printf("\n");
    return 0;
}
