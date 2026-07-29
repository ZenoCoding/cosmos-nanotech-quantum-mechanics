#include <stdio.h>
#include <math.h>

int main (void){
    int j; 
    for(j = 1; j < 20; j = j + 2){
        printf("\n%i %i", j, j*j);
    }
    return 0;
}