#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h>

void dvprint(double *v, int m);
void dmprint(double *matrix, int m, int n);
double* mmult(double *matrix1, double *matrix2, int m1, int p, int n2);

int main(void){
    double *m = calloc(6, sizeof(double));
    m[0] = 1.0; m[1] = 4.0; m[2] = 2.0;
    m[3] = 5.0; m[4] = 3.0; m[5] = 6.0;

    double *v = calloc(3, sizeof(double));
    v[0] = -1.2; v[1] = 0.8; v[2] = 3.6;

    //dmprint(m, 2, 3);

    double *r = mmult(m, v, 2, 3, 1);
    dmprint(r, 2, 1);
}

double* mmult(double *matrix1, double *matrix2, int m1, int p, int n2){
    double *r = calloc(m1 * n2, sizeof(double));

    for(int i = 0; i < m1; i++){
        for(int j = 0; j < n2; j++){
            for(int k = 0; k < p; k++){
                r[i + m1*j] += matrix1[i + m1*k] * matrix2[k + p*j];
            }
        }
    }

    return r;
}

void dvprint(double *v, int m)
{
    for (int i = 0; i < m; i++){
        printf("%12.6lf ", v[i]);
    }
    printf("\n");
}

void dmprint(double *matrix, int m, int n)
{
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            printf("%12.6lf ", matrix[i + m*j]);
        }
        printf("\n");
    }
}