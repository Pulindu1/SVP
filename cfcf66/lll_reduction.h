#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_basis(const char* msg, double** basis, int n, int len);

double dot_product(double* v1, double* v2, int len);

void subtract_vectors(double* v1, double* v2, double* result, int len);

void gram_schmidt(double** basis, double** orthogonal_basis, int n, int len);

void lll_reduction(double** basis, int n, int len, double delta);


