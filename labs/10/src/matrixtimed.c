#include <stdlib.h>
#include <stdio.h>
#include <tgmath.h>
#include <time.h>

typedef struct {
    struct timespec start;
    struct timespec end;
    double time;
} timertype;

void starttimer(timertype *timer);
void endtimer(timertype *timer);

void dvprint(double *v, int m);
void dmprint(double *matrix, int m, int n);
double* mmult(double *matrix1, double *matrix2, int m1, int p, int n2);

int main(void){
    timertype timer;

    int size = 2000;

    double *m = calloc(size*size, sizeof(double));
    double *v = calloc(size*size, sizeof(double));

    srand(100);

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            m[i + size*j] = (double) rand()/RAND_MAX;
            v[i + size*j] = (double) rand()/RAND_MAX;
        }
    }

    starttimer(&timer);

    double *r = mmult(m, v, size, size, size);
    endtimer(&timer);

    printf("That took %g seconds\n", timer.time);

    //dmprint(r, size, size);
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

void starttimer(timertype *timer){
    clock_gettime(CLOCK_MONOTONIC, &(*timer).start);
}

void endtimer(timertype *timer){
    clock_gettime(CLOCK_MONOTONIC, &(*timer).end);
    (*timer).time = (*timer).end.tv_sec - (*timer).start.tv_sec +
    ((*timer).end.tv_nsec - (*timer).start.tv_nsec) / 1e9;
}