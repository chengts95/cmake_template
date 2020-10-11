#pragma once

#include "common.h"


int dgetrf(int matrix_layout, int m, int n, Real* a, int lda);
int dgetrs(int matrix_layout, int n, int rhs, Real * a, int lda, Real *b, int ldb);
int dgetrf_pivot(int matrix_layout, int m, int n, Real * a, int lda, int * ipiv);
int dgetrf_pivot2(int matrix_layout, int m, int n, Real * a, int lda, int * ipiv);
int dgetrs_pivot(int matrix_layout, int n, int rhs, Real * a, int lda, int * ipiv, Real *b, int ldb);