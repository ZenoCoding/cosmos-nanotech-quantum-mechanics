#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void){
    unsigned int seed, steps, iter;
    printf("Enter seed, steps, and iterations: ");
    scanf("%u %u %u", &seed, &steps, &iter);
    for(int s = steps/10; s <= steps; s+= steps/10){
        int total = 0;
        for(int j = 0; j < iter; j++){
            int pos = 0;
            for(int i = 0; i < s; i++){
                pos += (double) rand()/RAND_MAX > .5 ? 1 : -1;
            }
            total += pos*pos;
        }
        printf("Average for %i steps    %g\n", s, (double) total / (double) iter);
    }
    return 0;
}