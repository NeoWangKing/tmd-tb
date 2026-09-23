#include <math.h>
#include <complex.h>
#include "config.h"
#include "mat.h"
#include "tb.h"
#include "herm3.h"

static double b1x, b1y, b2x, b2y;

static Mat C3; // 绕z轴旋转120° (C3)
static Mat Sh; // 关于xy平面的反射 (S_h)
static Mat Sv; // 关于yz平面的反射 (S_v)

static Mat NN_vec;  // 最近邻 (6个)
#if NNN_MODEL
static Mat NNN_vec; // 次近邻 (6个)
static Mat TNN_vec; // 第三近邻 (6个)
#endif

static void init_reciprocal(void)
{
    b1x = 2*PI; b1y = -2*PI/sqrt(3.0);
    b2x = 0.0;  b2y = 4*PI/sqrt(3.0);
}

static void init_symmetry(void)
{
    // C_3: (d_z^2, d_xy, d_x^2-y^2) 旋转
    C3 = mat_alloc(3, 3);
    MAT_AT(C3, 0, 0) = 1.0;         MAT_AT(C3, 0, 1) =  0.0;         MAT_AT(C3, 0, 2) =  0.0;
    MAT_AT(C3, 1, 0) = 0.0;         MAT_AT(C3, 1, 1) = +cos(4*PI/3); MAT_AT(C3, 1, 2) = -sin(4*PI/3);
    MAT_AT(C3, 2, 0) = 0.0;         MAT_AT(C3, 2, 1) = +sin(4*PI/3); MAT_AT(C3, 2, 2) = +cos(4*PI/3);

    // S_h: (d_z^2, d_xy, d_x^2-y^2) -> (d_z^2, d_xy, d_x^2-y^2)
    Sh = mat_alloc(3, 3);
    MAT_AT(Sh, 0, 0) = 1.0;         MAT_AT(Sh, 0, 1) =  0.0;         MAT_AT(Sh, 0, 2) =  0.0;
    MAT_AT(Sh, 1, 0) = 0.0;         MAT_AT(Sh, 1, 1) =  1.0;         MAT_AT(Sh, 1, 2) =  0.0;
    MAT_AT(Sh, 2, 0) = 0.0;         MAT_AT(Sh, 2, 1) =  0.0;         MAT_AT(Sh, 2, 2) =  1.0;

    // S_v: (d_z^2, d_xy, d_x^2-y^2) -> (d_z^2, -d_xy, d_x^2-y^2)
    Sv = mat_alloc(3, 3);
    MAT_AT(Sv, 0, 0) =  1.0;        MAT_AT(Sv, 0, 1) =  0.0;         MAT_AT(Sv, 0, 2) =  0.0;
    MAT_AT(Sv, 1, 0) =  0.0;        MAT_AT(Sv, 1, 1) = -1.0;         MAT_AT(Sv, 1, 2) =  0.0;
    MAT_AT(Sv, 2, 0) =  0.0;        MAT_AT(Sv, 2, 1) =  0.0;         MAT_AT(Sv, 2, 2) =  1.0;
}

static void init_vectors(void)
{
#if NN_MODEL
    // 最近邻 (NN)
    NN_vec  = mat_alloc(6, 2);
    MAT_AT(NN_vec, 0, 0) =  1.0; MAT_AT(NN_vec, 0, 1) =  0.0;
    MAT_AT(NN_vec, 1, 0) =  0.5; MAT_AT(NN_vec, 1, 1) =  sqrt(3.0)/2.0;
    MAT_AT(NN_vec, 2, 0) = -0.5; MAT_AT(NN_vec, 2, 1) =  sqrt(3.0)/2.0;
    MAT_AT(NN_vec, 3, 0) = -1.0; MAT_AT(NN_vec, 3, 1) =  0.0;
    MAT_AT(NN_vec, 4, 0) = -0.5; MAT_AT(NN_vec, 4, 1) = -sqrt(3.0)/2.0;
    MAT_AT(NN_vec, 5, 0) =  0.5; MAT_AT(NN_vec, 5, 1) = -sqrt(3.0)/2.0;

#if NNN_MODEL
    // 次近邻 (NNN)
    NNN_vec = mat_alloc(6, 2);
    MAT_AT(NNN_vec, 0, 0) =  1.5; MAT_AT(NNN_vec, 0, 1) =  sqrt(3.0)/2.0;
    MAT_AT(NNN_vec, 1, 0) =  0.0; MAT_AT(NNN_vec, 1, 1) =  sqrt(3.0);
    MAT_AT(NNN_vec, 2, 0) = -1.5; MAT_AT(NNN_vec, 2, 1) =  sqrt(3.0)/2.0;
    MAT_AT(NNN_vec, 3, 0) = -1.5; MAT_AT(NNN_vec, 3, 1) = -sqrt(3.0)/2.0;
    MAT_AT(NNN_vec, 4, 0) =  0.0; MAT_AT(NNN_vec, 4, 1) = -sqrt(3.0);
    MAT_AT(NNN_vec, 5, 0) =  1.5; MAT_AT(NNN_vec, 5, 1) = -sqrt(3.0)/2.0;

    // 第三近邻 (TNN)
    TNN_vec = mat_alloc(6, 2);
    MAT_AT(TNN_vec, 0, 0) =  2.0; MAT_AT(TNN_vec, 0, 1) =  0.0;
    MAT_AT(TNN_vec, 1, 0) =  1.0; MAT_AT(TNN_vec, 1, 1) =  sqrt(3.0);
    MAT_AT(TNN_vec, 2, 0) = -1.0; MAT_AT(TNN_vec, 2, 1) =  sqrt(3.0);
    MAT_AT(TNN_vec, 3, 0) = -2.0; MAT_AT(TNN_vec, 3, 1) =  0.0;
    MAT_AT(TNN_vec, 4, 0) = -1.0; MAT_AT(TNN_vec, 4, 1) = -sqrt(3.0);
    MAT_AT(TNN_vec, 5, 0) =  1.0; MAT_AT(TNN_vec, 5, 1) = -sqrt(3.0);
#endif // NNN_MODEL
#endif // NN_MODEL
}

// 对实数矩阵施加对称操作： M' = op * M * op^T
static void apply_symmetry(Mat Mout, Mat op, Mat Min)
{
    Mat tmp = mat_alloc(3, 3); Mat opT = mat_alloc(3, 3);
    mat_transpose(opT, op);
    mat_dot(tmp, op, Min);
    mat_dot(Mout, tmp, opT);
}

#if NN_MODEL
// 最近邻 E(R1) 
static void build_E_NN_R1(Mat E)
{
    MAT_AT(E, 0, 0) =  t11; MAT_AT(E, 0, 1) =  t12; MAT_AT(E, 0, 2) = t13;
    MAT_AT(E, 1, 0) = -t12; MAT_AT(E, 1, 1) =  t22; MAT_AT(E, 1, 2) = t23;
    MAT_AT(E, 2, 0) =  t13; MAT_AT(E, 2, 1) = -t23; MAT_AT(E, 2, 2) = t33;
}

#if NNN_MODEL
// 次近邻 E(~R1)
static void build_E_NNN_R1(Mat E)
{
    MAT_AT(E, 0, 0) =  r11;           MAT_AT(E, 0, 1) = r12; MAT_AT(E, 0, 2) = -r12/sqrt(3.0);
    MAT_AT(E, 1, 0) =  r21;           MAT_AT(E, 1, 1) = r22; MAT_AT(E, 1, 2) = r23;
    MAT_AT(E, 2, 0) = -r21/sqrt(3.0); MAT_AT(E, 2, 1) = r23; MAT_AT(E, 2, 2) = r22 + 2*r23/sqrt(3.0);
}

// 第三近邻 E(2R1)
static void build_E_TNN_R1(Mat E)
{
    MAT_AT(E, 0, 0) =  u11; MAT_AT(E, 0, 1) =  u12; MAT_AT(E, 0, 2) = u13;
    MAT_AT(E, 1, 0) = -u12; MAT_AT(E, 1, 1) =  u22; MAT_AT(E, 1, 2) = u23;
    MAT_AT(E, 2, 0) =  u13; MAT_AT(E, 2, 1) = -u23; MAT_AT(E, 2, 2) = u33;
}
#endif // NNN_NODEL
#endif // NN_NODEL

#if NN_MODEL
static void generate_hopping_matrices_NN(Mat mats[6])
{
    Mat C3sq = mat_alloc(3, 3);
    mat_dot(C3sq, C3, C3);
    Mat sigmaC3 = mat_alloc(3,3);
    mat_dot(sigmaC3, Sv, C3);
    Mat sigmaC3sq = mat_alloc(3,3);
    mat_dot(sigmaC3sq, Sv, C3sq);

    // R1
    build_E_NN_R1(mats[0]);
    // R2 = Sv*C3 -> R1
    apply_symmetry(mats[1], sigmaC3, mats[0]);
    // R3 = C3 -> R1
    apply_symmetry(mats[2], C3, mats[0]);
    // R4 = (R1)T
    apply_symmetry(mats[3], Sv, mats[0]);
    // R5 = C3^2 -> R1
    apply_symmetry(mats[4], C3sq, mats[0]);
    // R6 = Sv*C3*C3 -> R1
    apply_symmetry(mats[5], sigmaC3sq, mats[0]);
}

#if NNN_MODEL
static void generate_hopping_matrices_TNN(Mat mats[6])
{
    Mat C3sq = mat_alloc(3, 3);
    mat_dot(C3sq, C3, C3);
    Mat sigmaC3 = mat_alloc(3,3);   // σv * C3
    mat_dot(sigmaC3, Sv, C3);
    Mat sigmaC3sq = mat_alloc(3,3); // σv * C3²
    mat_dot(sigmaC3sq, Sv, C3sq);

    // R1
    build_E_TNN_R1(mats[0]);
    // R2 = Sv*C3 -> R1
    apply_symmetry(mats[1], sigmaC3, mats[0]);
    // R3 = C3 -> R1
    apply_symmetry(mats[2], C3, mats[0]);
    // R4 = (R1)T
    apply_symmetry(mats[3], Sv, mats[0]);
    // R5 = C3^2 -> R1
    apply_symmetry(mats[4], C3sq, mats[0]);
    // R6 = Sv*C3*C3 -> R1
    apply_symmetry(mats[5], sigmaC3sq, mats[0]);
}

static void generate_hopping_matrices_NNN(Mat mats[6])
{
    Mat C3sq = mat_alloc(3,3);
    mat_dot(C3sq, C3, C3);

    Mat E1 = mat_alloc(3,3);
    build_E_NNN_R1(E1);      // R̃1

    // R̃1
    mat_copy(mats[0], E1);
    // R̃4 = -R̃1 → E^T
    mat_transpose(mats[3], E1);
    // R̃2 = C3² · R̃4
    apply_symmetry(mats[1], C3sq, mats[3]);
    // R̃3 = C3 · R̃1
    apply_symmetry(mats[2], Sv, mats[0]);
    // R̃5 = C3² · R̃1
    apply_symmetry(mats[4], C3sq, mats[0]);
    // R̃6 = C3 · R̃4
    apply_symmetry(mats[5], Sv, mats[3]);
}
#endif // NNN_MODEL
#endif // NN_MODEL

// k 空间哈密顿量构造
void tb_frac_to_cart_k(double f1, double f2, double *kx, double *ky)
{
    *kx = f1 * b1x + f2 * b2x;
    *ky = f1 * b1y + f2 * b2y;
}

// 构建3x3哈密顿矩阵（无SOC）
void tb_build_Hk(const TB_Hop *hop, double kx, double ky, double complex H[3][3])
{
    // 初始化原位能（对角）
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            H[i][j] = 0.0;
        }
    }
    H[0][0] = e1;
    H[1][1] = e2;
    H[2][2] = e2;

    // 最近邻贡献
    for (int n=0; n<6; n++) {
        double phase = kx*MAT_AT(NN_vec, n, 0) + ky*MAT_AT(NN_vec, n, 1);
        double complex eikR = cos(phase) - I*sin(phase);
        for (int i=0;i<3;i++) {
            for (int j=0;j<3;j++) {
                H[i][j] += eikR * MAT_AT(hop->NN[n], i, j);
            }
        }
    }

#if NNN_MODEL
    // 次近邻贡献
    for (int n=0; n<6; n++) {
        double phase = kx*MAT_AT(NNN_vec, n, 0) + ky*MAT_AT(NNN_vec, n, 1);
        double complex eikR = cos(phase) - I*sin(phase);
        for (int i=0;i<3;i++) {
            for (int j=0;j<3;j++) {
                H[i][j] += eikR * MAT_AT(hop->NNN[n], i, j);
            }
        }
    }
    // 第三近邻贡献
    for (int n=0; n<6; n++) {
        double phase = kx*MAT_AT(TNN_vec, n, 0) + ky*MAT_AT(TNN_vec, n, 1);
        double complex eikR = cos(phase) - I*sin(phase);
        for (int i=0;i<3;i++) {
            for (int j=0;j<3;j++) {
                H[i][j] += eikR * MAT_AT(hop->TNN[n], i, j);
            }
        }
    }
#endif
}

void tb_init(TB_Hop *hop)
{
    init_reciprocal();
    init_symmetry();
    init_vectors();

    for (size_t i = 0; i < 6; ++i) {
        hop->NN[i]  = mat_alloc(3, 3);
        hop->NNN[i] = mat_alloc(3, 3);
        hop->TNN[i] = mat_alloc(3, 3);
    }
#if NN_MODEL
    generate_hopping_matrices_NN(hop->NN);
#if NNN_MODEL
    generate_hopping_matrices_NNN(hop->NNN);
    generate_hopping_matrices_TNN(hop->TNN);
#endif
#endif
}

#if SOC_MODEL
// 基于两个 3x3 块对角化，返回升序排列的 6 个本征值
void tb_eig_soc(const double complex H0[3][3], double eig[6])
{
    // Lz 矩阵（轨道顺序：d_z2, d_xy, d_x2-y2）
    static const double complex Lz[3][3] = {
        {0.0,   0.0,     0.0},
        {0.0,   0.0,  -2.0*I},
        {0.0, 2.0*I,     0.0}
    };

    double complex H_up[3][3], H_dn[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            H_up[i][j] = H0[i][j] + (lambda / 2.0) * Lz[i][j];
            H_dn[i][j] = H0[i][j] - (lambda / 2.0) * Lz[i][j];
        }
    }

    double eig_up[3], eig_dn[3];
    herm3_eigvals(H_up, eig_up);
    herm3_eigvals(H_dn, eig_dn);

    // 合并并整体排序（6 个值）
    for (int i = 0; i < 3; i++) eig[i]     = eig_up[i];
    for (int i = 0; i < 3; i++) eig[i + 3] = eig_dn[i];

    for (int i = 0; i < 5; i++) {
        for (int j = i + 1; j < 6; j++) {
            if (eig[i] > eig[j]) {
                double tmp = eig[i];
                eig[i] = eig[j];
                eig[j] = tmp;
            }
        }
    }
}

#if SOC_SPIN
// 分别计算自旋向上和向下的本征值（各自升序），不合并
void tb_eig_soc_spin(const double complex H0[3][3],
                     double eig_up[3], double eig_dn[3])
{
    // Lz 矩阵（轨道顺序：d_z2, d_xy, d_x2-y2）
    static const double complex Lz[3][3] = {
        {0.0, 0.0,      0.0},
        {0.0, 0.0, -2.0*I},
        {0.0, 2.0*I,     0.0}
    };

    double complex H_up[3][3], H_dn[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            H_up[i][j] = H0[i][j] + (lambda / 2.0) * Lz[i][j];
            H_dn[i][j] = H0[i][j] - (lambda / 2.0) * Lz[i][j];
        }
    }

    // 各自对角化（函数内部已做升序排序）
    herm3_eigvals(H_up, eig_up);
    herm3_eigvals(H_dn, eig_dn);
}
#endif // SOC_SPIN
#endif // SOC_MODEL
