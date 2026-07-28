#include <math.h>
#include <stdio.h>

int main(void)
{
    double x;
    unsigned int iterations;

    printf("Enter x, and # of iterations: ");
    if (scanf("%lf %u", &x, &iterations) != 2) {
        fprintf(stderr, "Expected a number and an iteration count.\n");
        return 1;
    }

    double sum = 1.0;
    double term = 1.0;

    for (unsigned int j = 1; j < iterations; ++j) {
        term *= x / j;
        sum += term;
        printf("%5u %12.8f\n", j, sum);
    }

    printf("exp(%g) = %.8f (libm: %.8f)\n", x, sum, exp(x));
    return 0;
}
