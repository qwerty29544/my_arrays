#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "ArrayAlgos/my_arrays.h"
#include "DiscreteShapes/cubed_cube.h"

int cube_test_case() {
    double center[3] = {0., 0., 0.};
    CUBE cube = init(center, 3., 20);
    print_cube(cube);
    cube_destructor(cube);    
    return 0;
}

int rcube_test_case() {
    double center[3] = {2., 4., 6.};
    double h = 0.25;
    double N_xyz[3] = {8, 10, 12};
    rCUBE rcube = init_rCUBE(center, h, N_xyz);
    print_rCUBE(rcube);
    printf("\n Distance between (0, 0, 0) and (7, 9, 11): %lf\n", rcube.distance[0][959]);
    destructor_rCUBE(rcube);    
    return 0;
}

int random_test_case() {
    int size = 4000;
    printf("max threads: %d\n", omp_get_max_threads());
    
    double complex *array = complex128_array_allocate_memory(size);
    complex128_array_fill_random_uniform(array, size, -5., 5., 123);
    complex128_array_print(array, size);
    complex128_array_free_memory(array);

    printf("\n ------------------------------- \n");

    float **matrix = float_matrix_allocate_memory(size, size);
    float_matrix_fill_random_uniform(matrix, size, size, -5, 5, 123);
    float_matrix_print(matrix, size, size);
    float_matrix_free_memory(matrix, size);

    printf("\n ------------------------------- \n");

    double complex **matrix_c = complex128_matrix_allocate_memory(size, size);
    complex128_matrix_fill_random_uniform(matrix_c, size, size, -5, 5, 123);
    complex128_matrix_print(matrix_c, size, size);
    complex128_matrix_free_memory(matrix_c, size);

    printf("\n ------------------------------- \n");

    float *array1 = float_array_allocate_memory(size);
    float *array2 = float_array_allocate_memory(size);
    float *array_res = float_array_allocate_memory(size);
    float_array_fill_random_uniform(array1, size, -2, 2, 123);
    float_array_fill_random_uniform(array2, size, -2, 2, 123);
    float_array_print(array1, size);
    float_array_print(array2, size);
    float_array_prod_array(array1, array2, array_res, size);     
    float_array_print(array_res, size);

    float_array_free_memory(array1);
    float_array_free_memory(array2);
    float_array_free_memory(array_res);


    return 0;
}

int main() {
    //cube_test_case();
    //rcube_test_case();
    random_test_case();
    return 0;
}