#include <stdio.h>
#include <math.h>

int main(){
    double a, b, c;
    printf("Enter a, b, c: ");
    scanf("%lf %lf %lf", &a, &b, &c);
    double disc = (b*b) - (4.0 * a * c);
    if (disc < 0) {
        printf("no valid solutions");
        return 1;
    }
    double sq = sqrt(disc);
    double x1 = (-b + sq)/(2.0*a);
    double x2 = (-b - sq)/(2.0*a);
    printf("x1 = %g, x2 = %g", x1, x2);
    return 0;
}