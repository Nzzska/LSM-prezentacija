#include <stdio.h>
#include <math.h>

int main(){
    return 0;
}

int jacobi(
    MATRIX_T A          /* in */,
    float    x[]       /* out */,
    float    b[]        /* in */,
    int      n          /* in */,
    float    tol        /* in */,
    int      max_iter   /* in */) {
        int   i,j;
        int   iter_num;
        float x_old[MAX_DIM];

        float distance(float x[], float y[], int n);

        /* Init x */
        iter_num = 0;
        do {
            iter_num++;

            /* saving old solution */
            for (i = 0; i < n; i++) {
                x_old[i] = x[i];
            }

            for (i = 0; i < n; i++) {
                x[i] = b[i];
                for (j = 0; j < i; j++) {
                    x[i] = x[i] - A[i][j] * x_old[j];
                }
                for (j = i + 1; j < n; j++) {
                    x[i] = x[i] - A[i][j] * x_old[j];
                }
                x[i] = x[i]/A[i][j];
            }
        } while ((iter_num < max_iter) &&
                 distance(x, x_old, n) >= tol));

        if (distance(x, x_old, n) < tol) {
            return 1;
        } else {
            return 0;
        }
    }

float distance(float x[], float y[], int n){
    int i;
    float sum = 0.0;

    for (i = 0; i < n; i++) {
        sum = sum + (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(sum);
}