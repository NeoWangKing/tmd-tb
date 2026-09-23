#ifndef TB_H_
#define TB_H_

#include <complex.h>
#include "mat.h"

// 面内晶格常数 (Å)：由 labelinfo 的三段路径长度反推得 3.182164，用于 k 换算到 Å⁻¹
#define A_ANG 3.182164

// 轨道基：d_z2, d_xy, d_x2-y2
#define TB_NORB 3

// 三个壳层的跳跃矩阵（各 6 个近邻矢量，矩阵实对称部分见 tb.c 的构造）
typedef struct {
    Mat NN[6];
    Mat NNN[6];
    Mat TNN[6];
} TB_Hop;

// 生成跳跃矩阵（内部初始化倒格矢、对称操作、近邻矢量）
void tb_init(TB_Hop *hop);

// 分数坐标 -> 笛卡尔 k（单位 1/a）
void tb_frac_to_cart_k(double f1, double f2, double *kx, double *ky);

// 构造 H(k)，kx/ky 为笛卡尔坐标（单位 1/a）
void tb_build_Hk(const TB_Hop *hop, double kx, double ky, double complex H[TB_NORB][TB_NORB]);

// 含 SOC 的本征值（H0 为无 SOC 的 3x3 块）
void tb_eig_soc(const double complex H0[TB_NORB][TB_NORB], double eig[6]);          // 6 个合并升序
void tb_eig_soc_spin(const double complex H0[TB_NORB][TB_NORB],
                     double eig_up[3], double eig_dn[3]);                           // 分自旋各自升序

#endif
