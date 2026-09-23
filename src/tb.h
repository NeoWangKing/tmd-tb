#ifndef TB_H_
#define TB_H_

#include <complex.h>
#include "mat.h"

// 轨道基：d_z2, d_xy, d_x2-y2
#define TB_NORB 3

// 面内晶格常数 (Å)，由 labelinfo 的三段路径长度反推
#define A_ANG 3.182164

// 模型参数（19 个）。defaults 由 params.h 给出，拟合时可以运行时替换。
typedef struct {
    double e1, e2;
    double t11, t12, t13, t22, t23, t33;
    double r11, r12, r21, r22, r23;
    double u11, u12, u13, u22, u23, u33;
} TB_Params;

// 三个壳层的跳跃矩阵（每个 6 个近邻矢量）
typedef struct {
    Mat NN[6];
    Mat NNN[6];
    Mat TNN[6];
} TB_Hop;

const TB_Params *tb_get_params(void);
void tb_set_params(const TB_Params *p);
void tb_init(TB_Hop *hop);
void tb_update_hopping(TB_Hop *hop);
void tb_frac_to_cart_k(double f1, double f2, double *kx, double *ky);
void tb_build_Hk(const TB_Hop *hop, double kx, double ky, double complex H[TB_NORB][TB_NORB]);
void tb_eig_soc(const double complex H0[TB_NORB][TB_NORB], double eig[6]);      // 合并升序
void tb_eig_soc_spin(const double complex H0[TB_NORB][TB_NORB],
                     double eig_up[3], double eig_dn[3]);                       // 分自旋各自升序

#endif
