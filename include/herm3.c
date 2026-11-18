#include <stddef.h>
#include <complex.h>
#include <math.h>
#include "herm3.h"

// 3x3 复数厄米矩阵 Jacobi 对角化（不依赖 LAPACK）
//   w: 升序本征值；V: 本征矢按列存放（正交归一），传 NULL 则只求本征值。
// 每个 (p,q) 先做一次酉变换把 a_pq 的相位消掉，再用实 Jacobi 旋转把它转掉。
// 本征矢的总体相位任意，做 BSE 时需另行固定。
void herm3_jacobi(const double complex A[3][3], double w[3], double complex V[3][3])
{
    double complex a[3][3];
    double complex v[3][3];

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            a[i][j] = A[i][j];
            v[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int sweep = 0; sweep < 64; ++sweep) {
        double off = 0.0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (i != j) off += cabs(a[i][j]);
            }
        }
        if (off < 1e-14) break;

        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                double t = cabs(a[p][q]);
                if (t < 1e-300) continue;

                // (1) 消相位：A <- D^dagger A D, D = diag(e^{i phi}, 1) 作用在 (p,q) 上
                double phi = carg(a[p][q]);
                double complex eiph = cos(phi) + I*sin(phi);
                for (int i = 0; i < 3; ++i) {
                    a[p][i] *= conj(eiph);  // 第 p 行 × e^{-i phi}
                    a[i][p] *= eiph;        // 第 p 列 × e^{+i phi}
                    v[i][p] *= eiph;        // 本征矢同步 V <- V D
                }

                // (2) 实 Jacobi 旋转：A <- R^T A R, 消掉 a[p][q]
                double theta = 0.5 * atan2(2.0 * t, creal(a[q][q]) - creal(a[p][p]));
                double c = cos(theta), s = sin(theta);
                for (int i = 0; i < 3; ++i) {           // 右乘 R
                    double complex aip = a[i][p], aiq = a[i][q];
                    a[i][p] = c * aip - s * aiq;
                    a[i][q] = s * aip + c * aiq;
                    double complex vip = v[i][p], viq = v[i][q];
                    v[i][p] = c * vip - s * viq;
                    v[i][q] = s * vip + c * viq;
                }
                for (int j = 0; j < 3; ++j) {           // 左乘 R^T
                    double complex apj = a[p][j], aqj = a[q][j];
                    a[p][j] = c * apj - s * aqj;
                    a[q][j] = s * apj + c * aqj;
                }
            }
        }
    }

    // 取对角元作为本征值，并做升序排序（本征矢的列同步交换）
    for (int i = 0; i < 3; ++i) w[i] = creal(a[i][i]);
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (w[i] > w[j]) {
                double tmp = w[i]; w[i] = w[j]; w[j] = tmp;
                for (int r = 0; r < 3; ++r) {
                    double complex tv = v[r][i]; v[r][i] = v[r][j]; v[r][j] = tv;
                }
            }
        }
    }

    if (V != NULL) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) V[i][j] = v[i][j];
        }
    }
}

// 对角化（只要本征值，升序）
void herm3_eigvals(const double complex H[3][3], double eig[3])
{
    herm3_jacobi(H, eig, NULL);
}
