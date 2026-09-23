#ifndef HERM3_H_
#define HERM3_H_

#include <complex.h>

// 3x3 复厄米矩阵 Jacobi 对角化（不依赖 LAPACK）
// V 按列存本征矢（正交归一，总体相位任意），传 NULL 则只求本征值
void herm3_jacobi(const double complex A[3][3], double w[3], double complex V[3][3]);
void herm3_eigvals(const double complex A[3][3], double w[3]);   // 只要本征值（升序）

#endif
