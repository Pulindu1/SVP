// include header
#include "lll_reduction.h"
// include the necessary module
#include <stdbool.h>

// calculate dot product
double dot_product(double* v1, double* v2, int len) {
    double sum = 0.0;
    for (int i = 0; i < len; i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

// pre define how to subtract vectors
void subtract_vectors(double* v1, double* v2, double* result, int len) {
    for (int i = 0; i < len; i++) {
        result[i] = v1[i] - v2[i];
    }
}

// vector normalisation
void normalize_vector(double* vector, int len) {
    double norm = sqrt(dot_product(vector, vector, len));
    for (int i = 0; i < len; i++) {
        vector[i] /= norm;
    }
}

//  gram schmidt function
void gram_schmidt(double** basis, double** orthogonal_basis, int n, int len) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < len; j++) {
            orthogonal_basis[i][j] = basis[i][j];
        }
        for (int k = 0; k < i; k++) {
            double mu = dot_product(basis[i], orthogonal_basis[k], len) /
                        dot_product(orthogonal_basis[k],
                                    orthogonal_basis[k], len);
            for (int j = 0; j < len; j++) {
                orthogonal_basis[i][j] -= mu * orthogonal_basis[k][j];
            }
        }
    }
}

// main lll function
void lll_reduction(double** basis, int n, int len, double delta) {
    // allocate enough memory
    double** orthogonal_basis = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        orthogonal_basis[i] = (double*)malloc(len * sizeof(double));
    }

    // perform gram schmidt
    gram_schmidt(basis, orthogonal_basis, n, len);

    // looping LLL algorithm
    for (int k = 1; k < n; k++) {
        // size reduction
        for (int j = 0; j < k; j++) {
            double mu_kj = dot_product(basis[k], orthogonal_basis[j], len) /
                           dot_product(orthogonal_basis[j],
                                       orthogonal_basis[j], len);
            if (fabs(mu_kj) > 0.5) {
                for (int l = 0; l < len; l++) {
                    basis[k][l] -= round(mu_kj) * basis[j][l];
                }
                gram_schmidt(basis, orthogonal_basis, n, len);
            }
        }

        // Exchange Condition
        double mu_kk1 = dot_product(basis[k], orthogonal_basis[k - 1], len) /
                        dot_product(orthogonal_basis[k - 1],
                                    orthogonal_basis[k - 1], len);
        if (dot_product(orthogonal_basis[k], orthogonal_basis[k], len) <
            (delta - mu_kk1 * mu_kk1) *
            dot_product(orthogonal_basis[k - 1],
                        orthogonal_basis[k - 1], len)) {
            double* temp = basis[k];
            basis[k] = basis[k - 1];
            basis[k - 1] = temp;
            // compute Gram-Schmidt again
            gram_schmidt(basis, orthogonal_basis, n, len);
            // keep k above 0
            k = fmax(1, k - 1);
        }
    }

    // Free the allocated memory
    for (int i = 0; i < n; i++) {
        free(orthogonal_basis[i]);
    }
    free(orthogonal_basis);
}
