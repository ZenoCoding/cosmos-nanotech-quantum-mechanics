//
// From the repository root, compile with:
//
//     make lab12
//
#include "utilities.h"

double modulussq(double complex z);
void probability(double complex *psi, double *prob, int m);

int main(void)
{
    //
    // m      : number of sites in 1d chain.
    // timemax: amount of time to time evolve for.
    // dtime  : time spacing.
    // fname  : file name to output results to.
    //
    int     m       = 20;
    double  timemax = 8.0 * M_PI;
    double  dtime   = 0.1;
    char   *fname   = "labs/12/data/timeevolution.txt";

    double complex *H    = calloc(m*m, sizeof(double complex));
    double complex *V    = calloc(m*m, sizeof(double complex));
    double         *w    = calloc(m  , sizeof(double        ));
    double complex *psi  = calloc(m  , sizeof(double complex));
    double         *prob = calloc(m  , sizeof(double        ));
    double complex *work = calloc(m  , sizeof(double complex));

    //
    // Setting up the 1d chain Hamiltonian.
    //
    double t = 1.0;
    double E = -4.0;
    double barrier_height = 2.0 * t;
    int barrier_site = m / 2;

    H[0 + m*0] = E ; H[0 + m*1] = -t;
    for (int i = 1 ; i < m-1 ; i = i + 1)
    {
        H[i + m*(i-1)] = -t ; H[i + m* i] =  E ; H[i + m*(i+1)] = -t;  
    }
    H[(m-1) + m*((m-1)-1)] = -t ; H[(m-1) + m*(m-1)] = E;

    // Central on-site barrier (about 50% transmission near the band center).
    H[barrier_site + m*barrier_site] += barrier_height;

    //
    // Diagonalize the Hamiltonian.
    //
    zhediagonalize(H, V, w, m);

    //
    // Open the file to output results into.
    //
    FILE *fout = fopen(fname, "w");


    //
    // Initialize wavefunction.
    //
    psi[0] = 1.0;


    //
    // Time 0 is automatic (no time evolution).
    // Carrying out the 0th step:
    //      ensure the wavefunction is normalized.
    //      compute the probability distribution from the wavefunction.
    //      print the probability distribution to the file.
    //
    znormalize(psi, m);
    probability(psi, prob, m);
    dvfprint(prob, m, fout);

    //
    // Time loop.
    // At each time step, time evolve the wavefunction a time dtime forward,
    // by updating psi = exp(-i*dtime*H) * psi. The rest of each iteration is
    // the same as the 0th one.
    //
    for (double time = dtime ; time <= timemax ; time = time + dtime)
    {
        zhediagexpmv(V, w, -I*dtime, psi, m, work);
        znormalize(psi, m);
        probability(psi, prob, m);
        dvfprint(prob, m, fout);
    }

    //
    // Clean everything up.
    //
    fclose(fout);
    free(H) ; free(V) ; free(w) ; free(psi) ; free(prob) ; free(work) ;

    return 0;
}

double modulussq(double complex z)
{
    return creal(z * conj(z));
}

void probability(double complex *psi, double *prob, int m)
{
    for (int i = 0 ; i < m ; i = i + 1)
    {
        prob[i] = modulussq(psi[i]);
    }
}
