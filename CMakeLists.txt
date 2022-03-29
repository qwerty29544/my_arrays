CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MyArrays)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -o 3 -lm -Wall -fopenmp")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -std=c11 -o 3 -lm -Wall -fopenmp")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -std=c11 -o 3 -lm -Wall -fopenmp")

SET(PATH_TO_LIB1 "./src/ArrayAlgos/")
SET(PATH_TO_LIB_CUBE "./src/DiscreteShapes/")

INCLUDE_DIRECTORIES(./src/)



ADD_EXECUTABLE(main src/main.c)

ADD_LIBRARY(my_arrays "${PATH_TO_LIB1}my_arrays.c" "${PATH_TO_LIB1}my_arrays.h")
ADD_LIBRARY(cubed_cube "${PATH_TO_LIB_CUBE}cubed_cube.c" "${PATH_TO_LIB_CUBE}cubed_cube.h")

TARGET_LINK_LIBRARIES(main my_arrays cubed_cube m)