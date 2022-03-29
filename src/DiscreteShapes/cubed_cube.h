#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>


typedef struct cubed_cube {
    double center[3];
    double lenght;
    double cubes_volume;
    int n_discr;
    double colloc_dist;
    double *****cubes;
    double ****colloc;
} CUBE;

CUBE init(double center[3], double lenght, int n_discr);
void print_cube(CUBE cube);
void cube_destructor(CUBE cube);

// ---------------------------------------------

typedef struct cubed_rectangular {
    int N_xyz[3];
    double center[3];
    double h;
    double length_xyz[3];
    double cubes_volume;
    double *****cubes;
    double ****collocations;
    double **distance;
} rCUBE;


rCUBE init_rCUBE(double center[3], double h, double N_xyz[3]);
void print_rCUBE(rCUBE rcube);
void destructor_rCUBE(rCUBE rcube);
