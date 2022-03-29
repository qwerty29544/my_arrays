#include "my_arrays.h"
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

// ----------------------------------------------------------------__
float *float_array_allocate_memory(int size) {
    float *array = (float*) malloc(size * sizeof(float));
    return array;
};

double *double_array_allocate_memory(int size) {
    double *array = (double*) malloc(size * sizeof(double));
    return array;
};

float complex *complex64_array_allocate_memory(int size) {
    float complex *array = (float complex*) malloc(size * sizeof(float complex));
    return array;
};


double complex *complex128_array_allocate_memory(int size) {
    double complex *array = (double complex*) malloc(size * sizeof(double complex));
    return array;
};


float **float_matrix_allocate_memory(int nrows, int ncols) {
    float **matrix = (float **) malloc(nrows * sizeof(float *));
    for (int index = 0; index < nrows; index++) {
        matrix[index] = (float *) malloc(ncols * sizeof(float));
    }
    return matrix;
};

double **double_matrix_allocate_memory(int nrows, int ncols) {
    double **matrix = (double **) malloc(nrows * sizeof(double *));
    for (int index = 0; index < nrows; index++) {
        matrix[index] = (double *) malloc(ncols * sizeof(double));
    }
    return matrix;
};

float complex **complex64_matrix_allocate_memory(int nrows, int ncols) {
    float complex **matrix = (float complex **) malloc(nrows * sizeof(float complex *));
    for (int index = 0; index < nrows; index++) {
        matrix[index] = (float complex *) malloc(ncols * sizeof(float complex));
    }
    return matrix;
};

double complex **complex128_matrix_allocate_memory(int nrows, int ncols) {
    double complex **matrix = (double complex **) malloc(nrows * sizeof(double complex *));
    for (int index = 0; index < nrows; index++) {
        matrix[index] = (double complex *) malloc(ncols * sizeof(double complex));
    }
    return matrix;
};

float ***float_3tensor_allocate_memory(int nrows, int ncols, int nmats) {
    float ***tensor = (float ***) malloc(nmats * sizeof(float **));
    for (int matrix = 0; matrix < nmats; matrix++) {
        tensor[matrix] = (float **) malloc(nrows * sizeof(float *));
        for (int row = 0; row < nrows; row++) {
            tensor[matrix][row] = (float *) malloc(ncols * sizeof(float));
        }
    }
    return tensor;
};

double ***double_3tensor_allocate_memory(int nrows, int ncols, int nmats) {
    double ***tensor = (double ***) malloc(nmats * sizeof(double **));
    for (int matrix = 0; matrix < nmats; matrix++) {
        tensor[matrix] = (double **) malloc(nrows * sizeof(double *));
        for (int row = 0; row < nrows; row++) {
            tensor[matrix][row] = (double *) malloc(ncols * sizeof(double));
        }
    }
    return tensor;
};

float complex ***complex64_3tensor_allocate_memory(int nrows, int ncols, int nmats) {
    float complex ***tensor = (float complex ***) malloc(nmats * sizeof(float complex **));
    for (int matrix = 0; matrix < nmats; matrix++) {
        tensor[matrix] = (float complex **) malloc(nrows * sizeof(float complex *));
        for (int row = 0; row < nrows; row++) {
            tensor[matrix][row] = (float complex *) malloc(ncols * sizeof(float complex));
        }
    }
    return tensor;
};

double complex ***complex128_3tensor_allocate_memory(int nrows, int ncols, int nmats) {
    double complex ***tensor = (double complex ***) malloc(nmats * sizeof(double complex **));
    for (int matrix = 0; matrix < nmats; matrix++) {
        tensor[matrix] = (double complex **) malloc(nrows * sizeof(double complex *));
        for (int row = 0; row < nrows; row++) {
            tensor[matrix][row] = (double complex *) malloc(ncols * sizeof(double complex));
        }
    }
    return tensor;
};

// -------------------------------------------------------------
// Free memory -------------------------------------------------
void float_array_free_memory(float *array) {
    free(array);
};

void double_array_free_memory(double *array) {
    free(array);
};

void complex64_array_free_memory(float complex *array) {
    free(array);
};

void complex128_array_free_memory(double complex *array) {
    free(array);
};

void float_matrix_free_memory(float **array, int nrows) {
    for (int index = 0; index < nrows; index++) {
        free(array[index]);
    }
    free(array);
};

void double_matrix_free_memory(double **array, int nrows) {
    for (int index = 0; index < nrows; index++) {
        free(array[index]);
    }
    free(array);
};

void complex64_matrix_free_memory(float complex **array, int nrows) {
    for (int index = 0; index < nrows; index++) {
        free(array[index]);
    }
    free(array);
};

void complex128_matrix_free_memory(double complex **array, int nrows) {
    for (int index = 0; index < nrows; index++) {
        free(array[index]);
    }
    free(array);
};

void float_3tensor_free_memory(float ***tensor, int nmats, int nrows) {
    for (int matrix = 0; matrix < nmats; matrix++) {
        for (int row = 0; row < nrows; row++) {
            free(tensor[matrix][row]);
        }
        free(tensor[matrix]);
    }
    free(tensor);
};

void double_3tensor_free_memory(double ***tensor, int nmats, int nrows) {
    for (int matrix = 0; matrix < nmats; matrix++) {
        for (int row = 0; row < nrows; row++) {
            free(tensor[matrix][row]);
        }
        free(tensor[matrix]);
    }
    free(tensor);
};

void complex64_3tensor_free_memory(float complex ***tensor, int nmats, int nrows) {
    for (int matrix = 0; matrix < nmats; matrix++) {
        for (int row = 0; row < nrows; row++) {
            free(tensor[matrix][row]);
        }
        free(tensor[matrix]);
    }
    free(tensor);
};

void complex128_3tensor_free_memory(double complex ***tensor, int nmats, int nrows) {
    for (int matrix = 0; matrix < nmats; matrix++) {
        for (int row = 0; row < nrows; row++) {
            free(tensor[matrix][row]);
        }
        free(tensor[matrix]);
    }
    free(tensor);
};

// -------------------------------------------------------------
// Apply func to two arrays ------------------------------------
void float_array_apply_array(float *array_1, float *array_2, float *array_res, int size, float (*function)(float, float)) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = function(array_1[index], array_2[index]);
    }
};

void double_array_apply_array(double *array_1, double *array_2, double *array_res, int size, double (*function)(double, double)) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = function(array_1[index], array_2[index]);
    }
};

void complex64_array_apply_array(float complex *array_1, float complex *array_2, float complex *array_res, int size, float complex (*function)(float complex, float complex)) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = function(array_1[index], array_2[index]);
    }
};

void complex128_array_apply_array(double complex *array_1, double complex *array_2, double complex *array_res, int size, double complex (*function)(double complex, double complex)) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = function(array_1[index], array_2[index]);
    }
};


void float_matrix_apply_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols, float (*function)(float, float)) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        float_array_apply_array(array_1[row], array_2[row], array_res[row], ncols, function);
    }
};

void double_matrix_apply_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols, double (*function)(double, double)) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        double_array_apply_array(array_1[row], array_2[row], array_res[row], ncols, function);
    }
};

void complex64_matrix_apply_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols, float complex (*function)(float complex, float complex)) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex64_array_apply_array(array_1[row], array_2[row], array_res[row], ncols, function);
    }
};

void complex128_matrix_apply_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols, double complex (*function)(double complex, double complex)) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex128_array_apply_array(array_1[row], array_2[row], array_res[row], ncols, function);
    }
};

// -------------------------------------------------------------

// -------------------------------------------------------------
// Sum two arrays ----------------------------------------------
void float_array_sum_array(float *array_1, float *array_2, float *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] + array_2[index];
    }
};

void double_array_sum_array(double *array_1, double *array_2, double *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] + array_2[index];
    }
};

void complex64_array_sum_array(float complex *array_1, float complex *array_2, float complex *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] + array_2[index];
    }
};

void complex128_array_sum_array(double complex *array_1, double complex *array_2, double complex *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] + array_2[index];
    }
};

void float_matrix_sum_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        float_array_sum_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void double_matrix_sum_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        double_array_sum_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex64_matrix_sum_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex64_array_sum_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex128_matrix_sum_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex128_array_sum_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

// -------------------------------------------------------------

// -------------------------------------------------------------
// Prod two arrays ---------------------------------------------
void float_array_prod_array(float *array_1, float *array_2, float *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] * array_2[index];
    }
};

void double_array_prod_array(double *array_1, double *array_2, double *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] * array_2[index];
    }
};

void complex64_array_prod_array(float complex *array_1, float complex *array_2, float complex *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] * array_2[index];
    }
};

void complex128_array_prod_array(double complex *array_1, double complex *array_2, double complex *array_res, int size) {
    int index;
    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] * array_2[index];
    }
};

void float_matrix_prod_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        float_array_prod_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void double_matrix_prod_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        double_array_prod_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex64_matrix_prod_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex64_array_prod_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex128_matrix_prod_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex128_array_prod_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

// -------------------------------------------------------------


// -------------------------------------------------------------
// Div two arrays ----------------------------------------------
void float_array_div_array(float *array_1, float *array_2, float *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] / array_2[index];
    }
};

void double_array_div_array(double *array_1, double *array_2, double *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] / array_2[index];
    }
};

void complex64_array_div_array(float complex *array_1, float complex *array_2, float complex *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] / array_2[index];
    }
};

void complex128_array_div_array(double complex *array_1, double complex *array_2, double complex *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] / array_2[index];
    }
};



void float_matrix_div_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        float_array_div_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void double_matrix_div_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        double_array_div_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex64_matrix_div_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex64_array_div_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex128_matrix_div_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex128_array_div_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

// -------------------------------------------------------------

// -------------------------------------------------------------
// Diff two arrays ---------------------------------------------
void float_array_diff_array(float *array_1, float *array_2, float *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] - array_2[index];
    }
};

void double_array_diff_array(double *array_1, double *array_2, double *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] - array_2[index];
    }
};

void complex64_array_diff_array(float complex *array_1, float complex *array_2, float complex *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] - array_2[index];
    }
};

void complex128_array_diff_array(double complex *array_1, double complex *array_2, double complex *array_res, int size) {
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array_res[index] = array_1[index] - array_2[index];
    }
};


void float_matrix_diff_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        float_array_diff_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void double_matrix_diff_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        double_array_diff_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex64_matrix_diff_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex64_array_diff_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

void complex128_matrix_diff_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols) {
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array_1, array_2, array_res) private(row)
    for (row = 0; row < nrows; row++) {
        complex128_array_diff_array(array_1[row], array_2[row], array_res[row], ncols);
    }
};

// -------------------------------------------------------------

// -------------------------------------------------------------
// array fill random value -------------------------------------
void float_array_fill_random_uniform(float *array, int size, float from, float to, int seed) {
    srand(seed);
    float normalizer =  (to - from) / (__INT_MAX__ * 1.);
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array[index] = rand() * normalizer + from;
    }
};

void double_array_fill_random_uniform(double *array, int size, double from, double to, int seed) {
    srand(seed);
    double normalizer =  (to - from) / (__INT_MAX__ * 1.);
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array[index] = rand() * normalizer + from;
    }
};

void complex64_array_fill_random_uniform(float complex *array, int size, float from, float to, int seed) {
    srand(seed);
    float normalizer =  (to - from) / (__INT_MAX__ * 1.);
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array[index] = (rand() * normalizer + from) + I * (rand() * normalizer + from);
    }
};

void complex128_array_fill_random_uniform(double complex *array, int size, double from, double to, int seed) {
    srand(seed);
    double normalizer =  (to - from) / (__INT_MAX__ * 1.);
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        array[index] = (rand() * normalizer + from) + I * (rand() * normalizer + from);
    }
};

void float_matrix_fill_random_uniform(float **array, int nrows, int ncols, float from, float to, int seed) {
    srand(seed);
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array) private(row)
    for (row = 0; row < nrows; row++) {
        float_array_fill_random_uniform(array[row], ncols, from, to, row + 1);
    }
};

void double_matrix_fill_random_uniform(double **array, int nrows, int ncols, double from, double to, int seed) {
    srand(seed);
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array) private(row)
    for (row = 0; row < nrows; row++) {
        double_array_fill_random_uniform(array[row], ncols, from, to, row + 1);
    }
};

void complex64_matrix_fill_random_uniform(float complex **array, int nrows, int ncols, float from, float to, int seed) {
    srand(seed);
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array) private(row)
    for (row = 0; row < nrows; row++) {
        complex64_array_fill_random_uniform(array[row], ncols, from, to, row + 1);
    }
};

void complex128_matrix_fill_random_uniform(double complex **array, int nrows, int ncols, double from, double to, int seed) {
    srand(seed);
    int row;

    omp_set_num_threads(1);
    #pragma omp parallel for shared(array) private(row)
    for (row = 0; row < nrows; row++) {
        complex128_array_fill_random_uniform(array[row], ncols, from, to, row + 1);
    }
};

// ----------------------------------------------------------------__
float float_array_dot_array(float *array_1, float *array_2, int size) {
    float result = 0.;
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

double double_array_dot_array(double *array_1, double *array_2, int size) {
    double result = 0.;
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

float complex complex64_array_dot_array(float complex *array_1, float complex *array_2, int size) {
    float complex result = 0.;
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

double complex complex128_array_dot_array(double complex *array_1, double complex *array_2, int size) {
    double complex result = 0.;
    int index;

    #pragma omp for
    for (index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

void float_matrix_array_multiplication(float **matrix, float *array, float *result, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);

    #pragma omp parallel for shared(matrix, array, result) private(row)
    for (row = 0; row < nrows; row++) {
        result[row] = float_array_dot_array(matrix[row], array, ncols);
    }
};

void double_matrix_array_multiplication(double **matrix, double *array, double *result, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);

    #pragma omp parallel for shared(matrix, array, result) private(row)
    for (row = 0; row < nrows; row++) {
        result[row] = double_array_dot_array(matrix[row], array, ncols);
    }
};

void complex64_matrix_array_multiplication(float complex **matrix, float complex *array, float complex *result, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);

    #pragma omp parallel for shared(matrix, array, result) private(row)
    for (row = 0; row < nrows; row++) {
        result[row] = complex64_array_dot_array(matrix[row], array, ncols);
    }
};

void complex128_matrix_array_multiplication(double complex **matrix, double complex *array, double complex *result, int nrows, int ncols) {
    int row;
    omp_set_num_threads(1);

    #pragma omp parallel for shared(matrix, array, result) private(row)
    for (row = 0; row < nrows; row++) {
        result[row] = complex128_array_dot_array(matrix[row], array, ncols);
    }
};

// ----------------------------------------------------------------__
     
void float_array_print(float *array, int size) {
    if (size > 0) {
        printf("[%.5f", array[0]);
        if (size <= 10) {
            for (int i = 1; i < size; i++) {
                printf(", %.5f", array[i]);
            }
        } else {
            for (int i = 1; i < 3; i++) {
                printf(", %.5f", array[i]);
            }
            printf(", ...");
            for (int i = size - 3; i < size; i++) {
                printf(", %.5f", array[i]);
            }
        }
        printf("]\n");
    } else {
        printf("None array \n");
    }
};

void double_array_print(double *array, int size) {
    if (size > 0) {
        printf("[%.5lf", array[0]);
        if (size <= 10) {
            for (int i = 1; i < size; i++) {
                printf(", %.5lf", array[i]);
            }
        } else {
            for (int i = 1; i < 3; i++) {
                printf(", %.5lf", array[i]);
            }
            printf(", ...");
            for (int i = size - 3; i < size; i++) {
                printf(", %.5lf", array[i]);
            }
        }
        printf("]\n");
    } else {
        printf("None array \n");
    }
};

void complex64_array_print(float complex *array, int size) {
    if (size > 0) {
        printf("[(%.4lf) + i(%.4lf)", crealf(array[0]), cimagf(array[0]));
        if (size <= 8) {
            for (int i = 1; i < size; i++) {
                printf(", (%.4lf) + i(%.4lf)", crealf(array[i]), cimagf(array[i]));
            }
        } else {
            for (int i = 1; i < 2; i++) {
                printf(", (%.4lf) + i(%.4lf)", crealf(array[i]), cimagf(array[i]));
            }
            printf(", ...");
            for (int i = size - 2; i < size; i++) {
                printf(", (%.4lf) + i(%.4lf)", crealf(array[i]), cimagf(array[i]));
            }
        }
        printf("]\n");
    } else {
        printf("None array \n");
    }
};

void complex128_array_print(double complex *array, int size) {
    if (size > 0) {
        printf("[(%.4lf) + i(%.4lf)", creal(array[0]), cimag(array[0]));
        if (size <= 8) {
            for (int i = 1; i < size; i++) {
                printf(", (%.4lf) + i(%.4lf)", creal(array[i]), cimag(array[i]));
            }
        } else {
            for (int i = 1; i < 2; i++) {
                printf(", (%.4lf) + i(%.4lf)", creal(array[i]), cimag(array[i]));
            }
            printf(", ...");
            for (int i = size - 2; i < size; i++) {
                printf(", (%.4lf) + i(%.4lf)", creal(array[i]), cimag(array[i]));
            }
        }
        printf("]\n");
    } else {
        printf("None array \n");
    }
};

void float_matrix_print(float **matrix, int nrows, int ncols) {
    if (nrows > 0 || ncols > 0) {
        printf("\nfloat matrix:\n");
        if (nrows <= 10) {
            for (int row = 0; row < nrows; row++) {
                printf("| ");
                float_array_print(matrix[row], ncols);
            }
        } else {
            for (int row = 0; row < 3; row++) {
                printf("| ");
                float_array_print(matrix[row], ncols);
            }
            printf("| ...\n");   
            for (int row = nrows - 3; row < nrows; row++) {
                printf("| ");
                float_array_print(matrix[row], ncols);
            }
        }
        printf("\n");
    
    } else {
        printf("\nNone matrix\n");
    }
};

void double_matrix_print(double **matrix, int nrows, int ncols) {
    if (nrows > 0 || ncols > 0) {
        printf("\ndouble matrix:\n");
        if (nrows <= 10) {
            for (int row = 0; row < nrows; row++) {
                printf("| ");
                double_array_print(matrix[row], ncols);
            }
        } else {
            for (int row = 0; row < 3; row++) {
                printf("| ");
                double_array_print(matrix[row], ncols);
            }
            printf("| ...\n");   
            for (int row = nrows - 3; row < nrows; row++) {
                printf("| ");
                double_array_print(matrix[row], ncols);
            }
        }
        printf("\n");
    
    } else {
        printf("\nNone matrix\n");
    }
};

void complex64_matrix_print(float complex **matrix, int nrows, int ncols) {
    if (nrows > 0 || ncols > 0) {
        printf("\ncmplx64 matrix:\n");
        if (nrows <= 8) {
            for (int row = 0; row < nrows; row++) {
                printf("| ");
                complex64_array_print(matrix[row], ncols);
            }
        } else {
            for (int row = 0; row < 3; row++) {
                printf("| ");
                complex64_array_print(matrix[row], ncols);
            }
            printf("| ...\n");   
            for (int row = nrows - 3; row < nrows; row++) {
                printf("| ");
                complex64_array_print(matrix[row], ncols);
            }
        }
        printf("\n");
    
    } else {
        printf("\nNone matrix\n");
    }
};

void complex128_matrix_print(double complex **matrix, int nrows, int ncols) {
    if (nrows > 0 || ncols > 0) {
        printf("\ncmplx128 matrix:\n");
        if (nrows <= 8) {
            for (int row = 0; row < nrows; row++) {
                printf("| ");
                complex128_array_print(matrix[row], ncols);
            }
        } else {
            for (int row = 0; row < 3; row++) {
                printf("| ");
                complex128_array_print(matrix[row], ncols);
            }
            printf("| ...\n");   
            for (int row = nrows - 3; row < nrows; row++) {
                printf("| ");
                complex128_array_print(matrix[row], ncols);
            }
        }
        printf("\n");
    
    } else {
        printf("\nNone matrix\n");
    }
};