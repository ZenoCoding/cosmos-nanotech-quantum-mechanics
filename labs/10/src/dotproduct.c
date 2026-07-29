#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h>

int main(void){
    int m = 3;
    int n = 1;
    double *x = calloc(m, sizeof(double));
    double *y = calloc(m, sizeof(double));

    for(int i = 0; i < m; i++){
        printf("Input x%d: ", i);
        scanf("%lf", &x[i]);
    }
    printf("\n");
    
    for(int i = 0; i < m; i++){
        printf("Input y%d: ", i);
        scanf("%lf", &y[i]);
    }
    printf("\n");

    double result = 0.;
    for(int i = 0; i < m; i++){
        result += x[i]*y[i];
    }

    printf("Dot product: %lf", result);
}