CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MyArrays)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall -fopenmp")
SET(PATH_TO_LIB1 "./src/ArrayAlgos/")


INCLUDE_DIRECTORIES(./src/)


ADD_EXECUTABLE(main src/main.c)

ADD_LIBRARY(my_arrays "${PATH_TO_LIB1}my_arrays.c" "${PATH_TO_LIB1}my_arrays.h")

TARGET_LINK_LIBRARIES(main my_arrays)