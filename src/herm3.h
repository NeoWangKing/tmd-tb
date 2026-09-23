#ifndef HERM3_H_
#define HERM3_H_

#include <complex.h>

// 3x3 复数厄米矩阵 Jacobi 对角化（不依赖 LAPACK）
//   w: 升序本征值；V: 本征矢按列存放（正交归一），传 NULL 则只求本征值。
// 每个 (p,q) 先做一次酉变换把 a_pq 的相位消掉，再用实 Jacobi 旋转把它转掉。
// 本征矢的总体相位任意，做 BSE 时需另行固定。
void herm3_jacobi(const double complex A[3][3], double w[3], double complex V[3][3]);

// 只要本征值（升序）
void herm3_eigvals(const double complex A[3][3], double w[3]);

#endif
