#include "cubed_cube.h"



CUBE init(double center[3], double lenght, int n_discr) {
    CUBE cube;
    int points_in_cube = 8;
    int coords_in_point = 3;
    double inc_by_coord = lenght/n_discr;
    double collocation_distances = sqrt(3 * inc_by_coord * inc_by_coord);

    for (int index = 0; index < 3; index++) {
        cube.center[index] = center[index];
    }

    cube.lenght = lenght;
    printf("\n%lf\n", cube.lenght);
    cube.n_discr = n_discr;
    printf("\n%d\n", cube.n_discr);
    cube.cubes_volume = inc_by_coord * inc_by_coord * inc_by_coord;
    cube.colloc_dist = collocation_distances;

    // Create cubes tensor n*n*n*8*3
    double *****array_cubes = (double *****) malloc(n_discr * sizeof(double ****));
    
    for (int x = 0; x < n_discr; x++) {
        
        array_cubes[x] = (double ****) malloc(n_discr * sizeof(double ***));
        
        for (int y = 0; y < n_discr; y++) {
            
            array_cubes[x][y] = (double ***) malloc(n_discr * sizeof(double **));
            
            for (int z = 0; z < n_discr; z++) {
                
                array_cubes[x][y][z] = (double **) malloc(points_in_cube * sizeof(double *));
                
                for (int point = 0; point < points_in_cube; point++) {
                    
                    array_cubes[x][y][z][point] = (double *) malloc(coords_in_point * sizeof(double));
                    
                    if (point == 0) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * x;
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * y;
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * z;
                    }

                    if (point == 1) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * y;
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * z;
                    }

                    if (point == 2) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * z;
                    }

                    if (point == 3) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * x;
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * z;
                    }

                    if (point == 4) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * x;
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * y;
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * (z + 1);
                    }

                    if (point == 5) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * y;
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * (z + 1);
                    }

                    if (point == 6) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * (z + 1);
                    }

                    if (point == 7) {
                        array_cubes[x][y][z][point][0] = center[0] - lenght / 2. + inc_by_coord * x;
                        array_cubes[x][y][z][point][1] = center[1] - lenght / 2. + inc_by_coord * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - lenght / 2. + inc_by_coord * (z + 1);
                    }
                
                }
            
            }
        
        }
    
    }    

    // assign tensor and struct field
    cube.cubes = array_cubes;

    double ****array_collocations = (double ****) malloc(n_discr * sizeof(double ***));
    
    for (int x = 0; x < n_discr; x++) {

        array_collocations[x] = (double ***) malloc(n_discr * sizeof(double **));

        for (int y = 0; y < n_discr; y++) {

            array_collocations[x][y] = (double **) malloc(n_discr * sizeof(double *));

            for (int z = 0; z < n_discr; z++) {

                array_collocations[x][y][z] = (double *) malloc(coords_in_point * sizeof(double));
                
                for (int coord = 0; coord < coords_in_point; coord++) {

                    array_collocations[x][y][z][coord] = array_cubes[x][y][z][0][coord] + inc_by_coord / 2.;

                }

            }

        }

    }

    cube.colloc = array_collocations;
    return cube;
}


void print_cube(CUBE cube) {
    int total_cubes = cube.n_discr * cube.n_discr * cube.n_discr;
    printf("\n ----------------------------------------------------------------------- \n");
    printf("Cube object has %d sub cubes with 8 points with 3 coords\n", total_cubes);
    printf("Total discrete shape size in memory: %.2lf MB\n", total_cubes * 8 * 3 * sizeof(double) * 1. / 1024. / 1024.);
    printf("Cubes volume is %.5lf\n\n", cube.cubes_volume);

    printf("Cube [%d, %d, %d]:\n", 0, 0, 0);
    for (int point = 0; point < 8; point++) {
        printf("\t(%.2lf, %.2lf, %.2lf)\n", 
                cube.cubes[0][0][0][point][0],
                cube.cubes[0][0][0][point][1],
                cube.cubes[0][0][0][point][2]);
    }

    printf("\n........\n");

    int N = cube.n_discr - 1;
    printf("Cube [%d, %d, %d]:\n", N, N, N);
    for (int point = 0; point < 8; point++) {
        printf("\t(%.2lf, %.2lf, %.2lf)\n", 
                cube.cubes[N][N][N][point][0],
                cube.cubes[N][N][N][point][1],
                cube.cubes[N][N][N][point][2]);
    }
    printf("\n ----------------------------------------------------------------------- \n");
};


void cube_destructor(CUBE cube) {
    for (int x = 0; x < cube.n_discr; x++) {
        for (int y = 0; y < cube.n_discr; y++) {
            for (int z = 0; z < cube.n_discr; z++) {
                for (int point = 0; point < 8; point++) {
                    free(cube.cubes[x][y][z][point]);
                }
                free(cube.cubes[x][y][z]);
                free(cube.colloc[x][y][z]);
            }
            free(cube.cubes[x][y]);
            free(cube.colloc[x][y]);
        }
        free(cube.cubes[x]);
        free(cube.colloc[x]);
    }
    free(cube.cubes);
    free(cube.colloc);
    printf("\nDestruction completed\n");
};

// -----------------------------------------------------------------------

rCUBE init_rCUBE(double center[3], double h, double N_xyz[3]) {
    rCUBE rcube;
    
    for (int dim = 0; dim < 3; dim++) {
        rcube.center[dim] = center[dim];
        rcube.N_xyz[dim] = N_xyz[dim];
        rcube.length_xyz[dim] = N_xyz[dim] * h;
    }

    rcube.h = h;
    rcube.cubes_volume = h * h * h;

    double *****array_cubes = (double *****) malloc(N_xyz[0] * sizeof(double ****)) ;
    double ****array_collocations = (double ****) malloc(N_xyz[0] * sizeof(double ***));
    
    for (int x = 0; x < N_xyz[0]; x++) {
    
        array_cubes[x] = (double ****) malloc(N_xyz[1] * sizeof(double ***));
        array_collocations[x] = (double ***) malloc(N_xyz[1] * sizeof(double **));
        
        for (int y = 0; y < N_xyz[1]; y++) {
    
            array_cubes[x][y] = (double ***) malloc(N_xyz[2] * sizeof(double **));
            array_collocations[x][y] = (double **) malloc(N_xyz[2] * sizeof(double *));

            for (int z = 0; z < N_xyz[2]; z++) {
    
                array_cubes[x][y][z] = (double **) malloc(8 * sizeof(double *));
                array_collocations[x][y][z] = (double *) malloc(3 * sizeof(double));

                for (int point = 0; point < 8; point++) {
    
                    array_cubes[x][y][z][point] = (double *) malloc(3 * sizeof(double));

                    if (point == 0) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * x;
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * y;
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * z;

                        array_collocations[x][y][z][0] = array_cubes[x][y][z][point][0] + h / 2;
                        array_collocations[x][y][z][1] = array_cubes[x][y][z][point][1] + h / 2;
                        array_collocations[x][y][z][2] = array_cubes[x][y][z][point][2] + h / 2;
    
                    }

                    if (point == 1) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * y;
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * z;
    
                    }

                    if (point == 2) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * z;
    
                    }

                    if (point == 3) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * x;
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * z;
    
                    }

                    if (point == 4) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * x;
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * y;
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * (z + 1);
    
                    }

                    if (point == 5) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * y;
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * (z + 1);
    
                    }

                    if (point == 6) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * (x + 1);
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * (z + 1);
    
                    }

                    if (point == 7) {
    
                        array_cubes[x][y][z][point][0] = center[0] - rcube.length_xyz[0] / 2. + h * x;
                        array_cubes[x][y][z][point][1] = center[1] - rcube.length_xyz[1] / 2. + h * (y + 1);
                        array_cubes[x][y][z][point][2] = center[2] - rcube.length_xyz[2] / 2. + h * (z + 1);
    
                    }

                }

            }

        }

    }

    double **array_distances = (double **) malloc(N_xyz[0] * N_xyz[1] * N_xyz[2] * sizeof(double *));

    for (int x = 0; x < N_xyz[0]; x++) {
        for (int y = 0; y < N_xyz[1]; y++) {
            for (int z = 0; z < N_xyz[2]; z++) {
                
                int position = z + N_xyz[2] * y + N_xyz[2] * N_xyz[1] * x;
                array_distances[position] = (double *) malloc(N_xyz[0] * N_xyz[1] * N_xyz[2] * sizeof(double));
                
                for (int x1 = 0; x1 < N_xyz[0]; x1++) {
                    for (int y1 = 0; y1 < N_xyz[1]; y1++) {
                        for (int z1 = 0; z1 < N_xyz[2]; z1++) {
                            
                            int position1 = z1 + N_xyz[2] * y1 + N_xyz[2] * N_xyz[1] * x1;
                            double distance_xyz[3] = {array_collocations[x][y][z][0] - array_collocations[x1][y1][z1][0],
                                                      array_collocations[x][y][z][1] - array_collocations[x1][y1][z1][1],
                                                      array_collocations[x][y][z][2] - array_collocations[x1][y1][z1][2]};
                            array_distances[position][position1] = sqrt(distance_xyz[0] * distance_xyz[0] + 
                                                                        distance_xyz[1] * distance_xyz[1] + 
                                                                        distance_xyz[2] * distance_xyz[2]);

                        }
                    }
                }

            }
        }
    }

    rcube.cubes = array_cubes;
    rcube.collocations = array_collocations;
    rcube.distance = array_distances;
    
    return rcube;
}


void print_rCUBE(rCUBE rcube) {
        int total_cubes = rcube.N_xyz[0] * rcube.N_xyz[1] * rcube.N_xyz[2];
    printf("\n ----------------------------------------------------------------------- \n");
    printf("Cube object has %d sub cubes with 8 points with 3 coords\n", total_cubes);
    printf("Total discrete shape size in memory: %.2lf MB\n", total_cubes * 8 * 3 * sizeof(double) * 1. / 1024. / 1024.);
    printf("Cubes volume is %.5lf\n\n", rcube.cubes_volume);

    printf("Cube [%d, %d, %d]:\n", 0, 0, 0);
    for (int point = 0; point < 8; point++) {
        printf("\t(%.2lf, %.2lf, %.2lf)\n", 
                rcube.cubes[0][0][0][point][0],
                rcube.cubes[0][0][0][point][1],
                rcube.cubes[0][0][0][point][2]);
    }

    printf("\n........\n");

    int Nx = rcube.N_xyz[0] - 1, 
        Ny = rcube.N_xyz[1] - 1, 
        Nz = rcube.N_xyz[2] - 1;
    printf("Cube [%d, %d, %d]:\n", Nx, Ny, Nz);
    for (int point = 0; point < 8; point++) {
        printf("\t(%.2lf, %.2lf, %.2lf)\n", 
                rcube.cubes[Nx][Ny][Nz][point][0],
                rcube.cubes[Nx][Ny][Nz][point][1],
                rcube.cubes[Nx][Ny][Nz][point][2]);
    }
    printf("\n ----------------------------------------------------------------------- \n");
}


void destructor_rCUBE(rCUBE rcube) {
    
    for (int x = 0; x < rcube.N_xyz[0]; x++) {
        
        for (int y = 0; y < rcube.N_xyz[1]; y++) {
            
            for (int z = 0; z < rcube.N_xyz[2]; z++) {
                
                for (int point = 0; point < 8; point++) {
                
                    free(rcube.cubes[x][y][z][point]);
                
                }

                free(rcube.cubes[x][y][z]);
                free(rcube.collocations[x][y][z]);
                free(rcube.distance[x * rcube.N_xyz[2] * rcube.N_xyz[1] + y * rcube.N_xyz[2] + z]);
            }

            free(rcube.cubes[x][y]);
            free(rcube.collocations[x][y]);
        }

        free(rcube.cubes[x]);
        free(rcube.collocations[x]);
    }

    free(rcube.cubes);
    free(rcube.collocations);
    free(rcube.distance);
    
    printf("\nDestruction completed\n");
}
