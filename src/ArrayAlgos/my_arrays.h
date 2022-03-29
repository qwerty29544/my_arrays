#include <complex.h>


typedef struct float_vector
{
    int size;
    float *array;
} fv;


typedef struct double_vector
{
    int size;
    double *array;
} dv;


typedef struct complex_vector64
{
    int size;
    float complex *array;
} cv64;


typedef struct complex_vector128
{
    int size;
    double complex *array;
} cv128;


typedef struct float_matrix 
{
    int nrows;
    int ncols;
    float **matrix;
} fm;


typedef struct double_matrix 
{
    int nrows;
    int ncols;
    double **matrix;
} dm;


typedef struct complex_matrix64
{
    int nrows;
    int ncols;
    float complex **matrix;
} cm64;


typedef struct complex_matrix128
{
    int nrows;
    int ncols;
    double complex **matrix;
} cm128;


// -------------------------------------------------------------
// Allocate memory ---------------------------------------------
float *float_array_allocate_memory(int size);
double *double_array_allocate_memory(int size);
float complex *complex64_array_allocate_memory(int size);
double complex *complex128_array_allocate_memory(int size);

float **float_matrix_allocate_memory(int nrows, int ncols);
double **double_matrix_allocate_memory(int nrows, int ncols);
float complex **complex64_matrix_allocate_memory(int nrows, int ncols);
double complex **complex128_matrix_allocate_memory(int nrows, int ncols);

float ***float_3tensor_allocate_memory(int nrows, int ncols, int nmats);
double ***double_3tensor_allocate_memory(int nrows, int ncols, int nmats);
float complex ***complex64_3tensor_allocate_memory(int nrows, int ncols, int nmats);
double complex ***complex128_3tensor_allocate_memory(int nrows, int ncols, int nmats);
//--------------------------------------------------------------

// -------------------------------------------------------------
// Free memory -------------------------------------------------
void float_array_free_memory(float *array);
void double_array_free_memory(double *array);
void complex64_array_free_memory(float complex *array);
void complex128_array_free_memory(double complex *array);

void float_matrix_free_memory(float **array, int nrows);
void double_matrix_free_memory(double **array, int nrows);
void complex64_matrix_free_memory(float complex **array, int nrows);
void complex128_matrix_free_memory(double complex **array, int nrows);

void float_3tensor_free_memory(float ***tensor, int nmats, int nrows);
void double_3tensor_free_memory(double ***tensor, int nmats, int nrows);
void complex64_3tensor_free_memory(float complex ***tensor, int nmats, int nrows);
void complex128_3tensor_free_memory(double complex ***tensor, int nmats, int nrows);
// -------------------------------------------------------------

// -------------------------------------------------------------
// Apply func to two arrays ------------------------------------
void float_array_apply_array(float *array_1, float *array_2, float *array_res, int size, float (*function)(float, float));
void double_array_apply_array(double *array_1, double *array_2, double *array_res, int size, double (*function)(double, double));
void complex64_array_apply_array(float complex *array_1, float complex *array_2, float complex *array_res, int size, float complex (*function)(float complex, float complex));
void complex128_array_apply_array(double complex *array_1, double complex *array_2, double complex *array_res, int size, double complex (*function)(double complex, double complex));

void float_matrix_apply_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols, float (*function)(float, float));
void double_matrix_apply_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols, double (*function)(double, double));
void complex64_matrix_apply_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols, float complex (*function)(float complex, float complex));
void complex128_matrix_apply_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols, double complex (*function)(double complex, double complex));
// -------------------------------------------------------------

// -------------------------------------------------------------
// Sum two arrays ----------------------------------------------
void float_array_sum_array(float *array_1, float *array_2, float *array_res, int size);
void double_array_sum_array(double *array_1, double *array_2, double *array_res, int size);
void complex64_array_sum_array(float complex *array_1, float complex *array_2, float complex *array_res, int size);
void complex128_array_sum_array(double complex *array_1, double complex *array_2, double complex *array_res, int size);

void float_matrix_sum_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols);
void double_matrix_sum_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols);
void complex64_matrix_sum_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols);
void complex128_matrix_sum_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols);
// -------------------------------------------------------------

// -------------------------------------------------------------
// Prod two arrays ---------------------------------------------
void float_array_prod_array(float *array_1, float *array_2, float *array_res, int size);
void double_array_prod_array(double *array_1, double *array_2, double *array_res, int size);
void complex64_array_prod_array(float complex *array_1, float complex *array_2, float complex *array_res, int size);
void complex128_array_prod_array(double complex *array_1, double complex *array_2, double complex *array_res, int size);

void float_matrix_prod_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols);
void double_matrix_prod_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols);
void complex64_matrix_prod_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols);
void complex128_matrix_prod_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols);
// -------------------------------------------------------------


// -------------------------------------------------------------
// Div two arrays ----------------------------------------------
void float_array_div_array(float *array_1, float *array_2, float *array_res, int size);
void double_array_div_array(double *array_1, double *array_2, double *array_res, int size);
void complex64_array_div_array(float complex *array_1, float complex *array_2, float complex *array_res, int size);
void complex128_array_div_array(double complex *array_1, double complex *array_2, double complex *array_res, int size);

void float_matrix_div_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols);
void double_matrix_div_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols);
void complex64_matrix_div_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols);
void complex128_matrix_div_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols);
// -------------------------------------------------------------

// -------------------------------------------------------------
// Diff two arrays ---------------------------------------------
void float_array_diff_array(float *array_1, float *array_2, float *array_res, int size);
void double_array_diff_array(double *array_1, double *array_2, double *array_res, int size);
void complex64_array_diff_array(float complex *array_1, float complex *array_2, float complex *array_res, int size);
void complex128_array_diff_array(double complex *array_1, double complex *array_2, double complex *array_res, int size);

void float_matrix_diff_matrix(float **array_1, float **array_2, float **array_res, int nrows, int ncols);
void double_matrix_diff_matrix(double **array_1, double **array_2, double **array_res, int nrows, int ncols);
void complex64_matrix_diff_matrix(float complex **array_1, float complex **array_2, float complex **array_res, int nrows, int ncols);
void complex128_matrix_diff_matrix(double complex **array_1, double complex **array_2, double complex **array_res, int nrows, int ncols);
// -------------------------------------------------------------

// -------------------------------------------------------------
// array fill random value -------------------------------------
void float_array_fill_random_uniform(float *array, int size, float from, float to, int seed);
void double_array_fill_random_uniform(double *array, int size, double from, double to, int seed);
void complex64_array_fill_random_uniform(float complex *array, int size, float from, float to, int seed);
void complex128_array_fill_random_uniform(double complex *array, int size, double from, double to, int seed);

void float_matrix_fill_random_uniform(float **array, int nrows, int ncols, float from, float to, int seed);
void double_matrix_fill_random_uniform(double **array, int nrows, int ncols, double from, double to, int seed);
void complex64_matrix_fill_random_uniform(float complex **array, int nrows, int ncols, float from, float to, int seed);
void complex128_matrix_fill_random_uniform(double complex **array, int nrows, int ncols, double from, double to, int seed);

// -------------------------------------------------------------
// Array array dot product -------------------------------------
float float_array_dot_array(float *array_1, float *array_2, int size);
double double_array_dot_array(double *array_1, double *array_2, int size);
float complex complex64_array_dot_array(float complex *array_1, float complex *array_2, int size);
double complex complex128_array_dot_array(double complex *array_1, double complex *array_2, int size);

// Matrix vector multiplication --------------------------------
void float_matrix_array_multiplication(float **matrix, float *array, float *result, int nrow, int ncol);
void double_matrix_array_multiplication(double **matrix, double *array, double *result, int nrow, int ncol);
void complex64_matrix_array_multiplication(float complex **matrix, float complex *array, float complex *result, int nrow, int ncol);
void complex128_matrix_array_multiplication(double complex **matrix, double complex *array, double complex *result, int nrow, int ncol);
// -------------------------------------------------------------
     
void float_array_print(float *array, int size);
void double_array_print(double *array, int size);
void complex64_array_print(float complex *array, int size);
void complex128_array_print(double complex *array, int size);

void float_matrix_print(float **matrix, int nrows, int ncols);
void double_matrix_print(double **matrix, int nrows, int ncols);
void complex64_matrix_print(float complex **matrix, int nrows, int ncols);
void complex128_matrix_print(double complex **matrix, int nrows, int ncols);