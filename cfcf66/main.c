// include the necessary modules
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

// inclide headers
#include "lll_reduction.h"
#include "enumerate.h"

// create a structure for parsing vectors
typedef struct {
    double *elements;
    int length;
} Vector;

// declare the function
double euclidean_len(double *v, int len);

// applies an additional processing step to further reduce the basis
void further_processing(double **basis, int dim, int len) {
    int maxIterations = 100;
    for (int iter = 0; iter < maxIterations; iter++) {
        bool madeProgress = false;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (i != j) {
                    double coeff = round(dot_product(basis[i], basis[j], len) /
                                    dot_product(basis[j], basis[j], len));
                    if (coeff != 0.0) {
                        double originalLength = euclidean_len(basis[i], len);
                        for (int k = 0; k < len; k++) {
                            basis[i][k] -= coeff * basis[j][k];
                        }
                        double newLength = euclidean_len(basis[i], len);
                        if (newLength < originalLength) {
                            madeProgress = true;
                        } else {
                            for (int k = 0; k < len; k++) {
                                basis[i][k] += coeff * basis[j][k];
                            }
                        }
                    }
                }
            }
        }
        if (!madeProgress) {
            break;
        }
    }
}

// parsing function
Vector parse_vector(const char* str) {
    // finds length of vector
    int length = 0;
    int inNumber = 0;  // used for checks

    for (int i = 0; str[i]; i++) {
        // finds each vector, if it contains [], it is a vector
        if (str[i] != ' ' && str[i] != '[' && str[i] != ']') {
            if (!inNumber) {
                inNumber = 1;
                length++;
            }
        } else {
            inNumber = 0;
        }
    }
    // define vector struct
    Vector vector;
    vector.length = length;
    vector.elements = (double *)malloc(vector.length * sizeof(double));

    // gets elements from each string
    const char *token = strtok((char *)str, " []");
    int index = 0;
    while (token != NULL && index < length) {
        vector.elements[index++] = atof(token);
        token = strtok(NULL, " []");
    }

    return vector;
}

// function to print matrix in correct format
void print_matrix(double **matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// main function
int main(int argc, char *argv[]) {
    if (argc <= 1) {
        printf("Usage: ./runme [vector1] [vector2] ...\n");
        return 1;
    }

    // Temporary buffer, used for concatenating vector elements
    char temporaryBuffer[1024] = {0};
    int numOfVectors = 0;
    // allocate space for vectors
    Vector *vectors = (Vector *)malloc(argc * sizeof(Vector));

    // processes command line arguments, parses vectors
    for (int i = 1; i < argc; i++) {
        snprintf(temporaryBuffer + strlen(temporaryBuffer),
            sizeof(temporaryBuffer) - strlen(temporaryBuffer),
            "%s ", argv[i]);

        if (strchr(argv[i], ']') != NULL) {
            // used to debug
            printf("Concatenated String: %s\n", temporaryBuffer);
            vectors[numOfVectors] = parse_vector(temporaryBuffer);
            // used to debug
            printf("Parsed Vector: ");
            // printf("num of vectors = %d \n", numberOfVectors);
            for (int j = 0; j < vectors[numOfVectors].length; j++) {
                printf("%f ", vectors[numOfVectors].elements[j]);
            }
            printf("\n");
            numOfVectors++;
            memset(temporaryBuffer, 0, sizeof(temporaryBuffer));
        }
    }

    int vector_len = vectors[0].length;

    // allocate memory to basis vectors
    double **basis = (double **)malloc(numOfVectors * sizeof(double *));
    for (int i = 0; i < numOfVectors; i++) {
        basis[i] = (double *)malloc(vector_len * sizeof(double));
        for (int j = 0; j < vector_len; j++) {
            basis[i][j] = vectors[i].elements[j];
        }
    }

    // Start timing
    clock_t start = clock();

    // print the original basis
    printf("Original Basis:\n");
    print_matrix(basis, numOfVectors, vector_len);

    // set delta parameter (for lll)
    double delta = 0.75;

    // call to LLL reduction function, then print
    lll_reduction(basis, numOfVectors, vector_len, delta);
    printf("Reduced Basis:\n");
    print_matrix(basis, numOfVectors, vector_len);

    // process the basis further to shorten the vectors further, then print
    further_processing(basis, numOfVectors, vector_len);
    printf("Post-Processed Basis:\n");
    print_matrix(basis, numOfVectors, vector_len);

    // finally, apply enumeration function
    printf("Enumerating for the shortest vector in the lattice...\n");
    enumeration(basis, numOfVectors, vector_len);

    // time end of program
    clock_t end = clock();
    // print time taken:
    double time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken to run: %f seconds\n", time_taken);



    // free the allocated memory
    for (int i = 0; i < numOfVectors; i++) {
        free(vectors[i].elements);
    }
    free(vectors);
    for (int i = 0; i < numOfVectors; i++) {
        free(basis[i]);
    }
    free(basis);

    return 0;
}
