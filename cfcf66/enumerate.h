#ifndef ENUMERATE_H
#define ENUMERATE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>


double euclidean_len(double *v, int len);

double len_squared(double *v, int len);

void add_scaled_vector(double *dest, double *src, double scalar, int len);

void enum_recursive(double **basis, int rows, int cols, double *currentVector, double currentLen, int depth, double *shortestVector, double *minLen, int limits);

void enumeration(double **basis, int rows, int cols);

#endif // ENUMERATE_H
