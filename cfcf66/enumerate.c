// include the necessary module
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

// include header file
#include "enumerate.h"


// calculate the euclidean length of a vector
double euclidean_len(double *v, int len) {
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

// find the square length of vector
double len_squared(double *v, int len) {
    double sum = 0;
    for (int i = 0; i < len; i++) {
        sum += v[i] * v[i];
    }
    return sum;
}

// defines scaled vector addition
void add_scaled_vector(double *dest, double *src, double scalar, int len) {
    for (int i = 0; i < len; i++) {
        dest[i] += src[i] * scalar;
    }
}

// recursivly enumerate
void enum_recursive(double **basis, int rows, int cols,
                    double *currentVector, double currentLen,
                    int depth, double *shortestVector,
                    double *minLen, int boundary) {
    if (depth == rows) {
        if (currentLen > 0 && currentLen < *minLen) {
            *minLen = currentLen;
            memcpy(shortestVector, currentVector, cols * sizeof(double));
        }
        return;
    }

    double *newVector = (double *)malloc(cols * sizeof(double));
    for (int i = -boundary; i <= boundary; i++) {
        memcpy(newVector, currentVector, cols * sizeof(double));
        add_scaled_vector(newVector, basis[depth],
                          (double)i, cols);
        double newLen = len_squared(newVector, cols);
        if (newLen < *minLen) {
            enum_recursive(basis, rows, cols, newVector, newLen,
                           depth + 1, shortestVector, minLen, boundary);
        }
    }
    free(newVector);
}

// main enum function
void enumeration(double **basis, int rows, int cols) {
    double minLen = DBL_MAX;
    double *shortestVector = (double *)malloc(cols * sizeof(double));
    double *currentVector = (double *)calloc(cols, sizeof(double));

    int boundary = 10;

    enum_recursive(basis, rows, cols, currentVector, 0, 0,
                   shortestVector, &minLen, boundary);

    printf("Shortest Vector: ");
    for (int i = 0; i < cols; i++) {
        printf("%f ", shortestVector[i]);
    }
    printf("\n");

    double length = euclidean_len(shortestVector, cols);
    printf("Euclidean Length of Shortest Vector: %f\n", length);

    // create result.txt file
    FILE *file = fopen("result.txt", "w");
    if (file != NULL) {
        fprintf(file, "%f", length);
        fclose(file);
    } else {
        printf("Error opening file!\n");
    }

    free(shortestVector);
    free(currentVector);
}
