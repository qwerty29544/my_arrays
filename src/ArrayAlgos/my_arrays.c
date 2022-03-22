#include "my_arrays.h"
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

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

// ----------------------------------------------------------------__
float float_array_dot_array(float *array_1, float *array_2, int size) {
    float result = 0.;
    for (int index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

double double_array_dot_array(double *array_1, double *array_2, int size) {
    double result = 0.;
    for (int index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

float complex complex64_array_dot_array(float complex *array_1, float complex *array_2, int size) {
    float complex result = 0.;
    for (int index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

double complex complex128_array_dot_array(double complex *array_1, double complex *array_2, int size) {
    double complex result = 0.;
    for (int index = 0; index < size; index++) {
        result += array_1[index] * array_2[index];
    }
    return result;
};

void float_matrix_array_multiplication(float **matrix, float *array, float *result, int nrows, int ncols) {
    for (int row = 0; row < nrows; row++) {
        result[row] = float_array_dot_array(matrix[row], array, ncols);
    }
};

void double_matrix_array_multiplication(double **matrix, double *array, double *result, int nrows, int ncols) {
    for (int row = 0; row < nrows; row++) {
        result[row] = double_array_dot_array(matrix[row], array, ncols);
    }
};

void complex64_matrix_array_multiplication(float complex **matrix, float complex *array, float complex *result, int nrows, int ncols) {
    for (int row = 0; row < nrows; row++) {
        result[row] = complex64_array_dot_array(matrix[row], array, ncols);
    }
};

void complex128_matrix_array_multiplication(double complex **matrix, double complex *array, double complex *result, int nrows, int ncols) {
    for (int row = 0; row < nrows; row++) {
        result[row] = complex128_array_dot_array(matrix[row], array, ncols);
    }
};

// ----------------------------------------------------------------__
dv init_dv(double* vector, int len);
void print_dv(dv obj);
void free_dv(dv str);
dv arange_dv(double start, double stop, double step);
dv fill_dv(double fill_value, int len);        