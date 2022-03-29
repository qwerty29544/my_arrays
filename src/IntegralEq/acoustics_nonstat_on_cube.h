#include <complex.h>

typedef struct time_vector {
    double complex ***elements;
    int tau_max;
    int size;
    int num_coords;
    double *timeline_vector;
    double delta_time;
    double time_zero;
    int size_of_timeline;
    int current_index;
    char *name;
} TIME_VECTOR;

TIME_VECTOR init_TIME_VECTOR(int tau_max, int size, int num_coords, double delta_time, double time_zero, int size_of_timeline, char *name);
void print_TIME_VECTOR(TIME_VECTOR tv);
void destructor_TIME_VECTOR(TIME_VECTOR tv);
double *get_element_TIME_VECTOR(TIME_VECTOR tv, int component_number, double tau);