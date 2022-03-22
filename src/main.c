#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ArrayAlgos/my_arrays.h"

int main() {
    srand(time(NULL));
    int size = 90;
    double *vector = (double*) malloc(size * sizeof(double));
    double_array_fill_random_uniform(vector, size, -2, 2);

    dv example1 = init_dv(vector, size);
    print_dv(example1);

    free_dv(example1);

    return 0;
}