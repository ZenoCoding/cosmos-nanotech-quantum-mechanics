#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(){
    double R;
    unsigned int seed;
    int iter, moments=5;
    printf("Enter nonnegative seed, and iteration ct: ");
    scanf("%u %i", &seed, &iter);
    srand(seed);
    double totals[5] = {0};
    for(int i = 0; i < iter; i++){
        R = (double) rand()/RAND_MAX;
        printf("R: %lf\n", R);
        for(int j = 0; j < moments; j++){
            totals[j] += pow(R, j+1);
        }
    }

    for(int j = 0; j < moments; j++){
        printf("moment %i: %g \n", j+1, totals[j]/(double) iter);
    }
    return 0;
}