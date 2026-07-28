#include <math.h>
#include <stdio.h>

int main(void){
    double x, sum;
    long int fact;
    long double power;
    unsigned int iter;
    printf("Enter x, and # of iterations: ");
    scanf("%lf %u", &x, &iter);
    for(int i = 0; i < iter; i++){

    }

    sum = 1.;
    fact = 1;
    power = 1;

    for(int j = 1; j < iter; j++){
        fact *= j;
        power *= x;
        sum = sum + power/fact;
        printf("\n %5i %12.8lf", j, sum);
    }
    printf("\n");
    return 0;
}

double sin(int x){
    unsigned int 
}