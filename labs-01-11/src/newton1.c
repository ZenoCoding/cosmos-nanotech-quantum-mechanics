#include <stdio.h>
#include <math.h>

int main(){
    int j;
    double x=1., A;
    printf("Enter A: ");
    scanf("%lf", &A);
    for(j=0; j<20; j=j+1){
        x=x/2. + A/(2.*x);
        printf("%5i %g \n", j, x);
    }
    return 0;
}