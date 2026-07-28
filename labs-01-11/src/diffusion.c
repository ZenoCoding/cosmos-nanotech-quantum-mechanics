#include <stdio.h>
#include <math.h>

int main(void) {
    FILE * fileout;
    int x, t, Nt, Nx=1000;
    double rho[Nx], newrho[Nx], D, dt, dx, Ddttodx2, norm, realx, realt;
    
    fileout=fopen("data/diffusion.txt", "w");

    printf("\nEnter D, dt, dx: ");
    scanf("%lf %lf %lf", &D, &dt, &dx);
    Ddttodx2 = D*dt /(dx*dx);
    printf("\nEnter number of time steps: ");
    scanf("%d", &Nt);

    for (x=0; x<Nx; x++){
        rho[x] = ((double) x)/Nx/500.0;
    }

    for(t=0; t<Nt; t++)
    {
        for(x=1; x<Nx-1; x++){
            newrho[x]=rho[x]+Ddttodx2*(rho[x+1]+rho[x-1]-2.0*rho[x]);
        }
        norm=0.0;
        for (x=1; x<Nx-1; x++){
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
        printf("\n checking no particles are lost! norm = %8.3lf  ", norm);
    }
    
    fclose(fileout);
    printf("\n");
    return 0;
}