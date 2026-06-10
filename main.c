#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <lapacke.h> 
#include "mat.h"

#define NN_MODEL 1
#ifndef NNN_MODEL
#define NNN_MODEL 0
#endif 
#ifndef SOC_MODEL
#define SOC_MODEL 0
#endif
#ifndef SOC_SPIN
#define SOC_SPIN 0
#endif
#include "params.h"

#define PI 3.14159265358979323846

static double b1x, b1y, b2x, b2y;

static Mat C3; // 绕z轴旋转120° (C3)
static Mat Sh; // 关于xy平面的反射 (S_h)
static Mat Sv; // 关于yz平面的反射 (S_v)

static Mat NN_vec;  // 最近邻 (6个)
static Mat NNN_vec; // 次近邻 (6个)
static Mat TNN_vec; // 第三近邻 (6个)

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
void apply_symmetry(Mat Mout, Mat op, Mat Min)
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
    apply_symmetry(mats[3], Sv, mats[0]);   // σv * HR1 * σv^T
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
    apply_symmetry(mats[3], Sv, mats[0]);   // σv * HR1 * σv^T
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
static void frac_to_cart_k(double f1, double f2, double *kx, double *ky)
{
    *kx = f1 * b1x + f2 * b2x;
    *ky = f1 * b1y + f2 * b2y;
}

// 构建3x3哈密顿矩阵（无SOC）
static void build_Hk(double kx, double ky, Mat NN_mats[6], Mat NNN_mats[6], Mat TNN_mats[6], double complex H[3][3])
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
                H[i][j] += eikR * MAT_AT(NN_mats[n], i, j);
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
                H[i][j] += eikR * MAT_AT(NNN_mats[n], i, j);
            }
        }
    }
    // 第三近邻贡献
    for (int n=0; n<6; n++) {
        double phase = kx*MAT_AT(TNN_vec, n, 0) + ky*MAT_AT(TNN_vec, n, 1);
        double complex eikR = cos(phase) - I*sin(phase);
        for (int i=0;i<3;i++) {
            for (int j=0;j<3;j++) {
                H[i][j] += eikR * MAT_AT(TNN_mats[n], i, j);
            }
        }
    }
#endif
}

// 输入升序排列的能带 eig，以及能带总数 nband
// 返回价带顶能量 (VBM)
static double find_vbm_from_sorted_eigenvalues(double *eig, int nband) {
    // 找到相邻能级差最大的位置，作为带隙
    double max_gap = 0.0;
    int gap_idx = -1;
    for (int i = 0; i < nband-1; i++) {
        double gap = eig[i+1] - eig[i];
        if (gap > max_gap) {
            max_gap = gap;
            gap_idx = i;
        }
    }
    // gap_idx 是价带顶的索引，价带是从 0 到 gap_idx (包括)
    if (gap_idx < 0) return eig[0]; // 无带隙时取最低能带？不太可能
    return eig[gap_idx]; // VBM 是价带的最大值
}

// 对角化
static void diagonalize_3x3_hermitian(double complex H[3][3], double eig[3]) {
    char jobz = 'N';
    char uplo = 'U';
    int n = 3;
    double w[3];
    double complex a[3*3];
    // 行优先
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            a[i*3 + j] = H[i][j];
        }
    }
    LAPACKE_zheev(LAPACK_ROW_MAJOR, jobz, uplo, n, a, n, w);
    eig[0] = w[0]; eig[1] = w[1]; eig[2] = w[2];

    for (int i=0; i<2; i++) {
        for (int j=i+1; j<3; j++) {
            if (eig[i] > eig[j]) {
                double tmp = eig[i];
                eig[i] = eig[j];
                eig[j] = tmp;
            }
        }
    }
}

#if SOC_MODEL
// 基于两个 3x3 块对角化，返回升序排列的 6 个本征值
static void diagonalize_6x6_with_soc(double complex H0[3][3], double eig[6])
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
    diagonalize_3x3_hermitian(H_up, eig_up);
    diagonalize_3x3_hermitian(H_dn, eig_dn);

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
static void diagonalize_soc_spin_blocks(double complex H0[3][3],
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
    diagonalize_3x3_hermitian(H_up, eig_up);
    diagonalize_3x3_hermitian(H_dn, eig_dn);
}
#endif // SOC_SPIN
#endif // SOC_MODEL

typedef struct {
    const double *start;
    const double *end;
    int steps;
} Segment;

// 高对称点路径定义
static const double GAMMA[2] = {0.0, 0.0};          // Gamma 点
static const double M[2]     = {0.5, 0.5};          // M 点
static const double K[2]     = {2.0/3.0, 1.0/3.0};  // K 点
static const double Mp[2]    = {0.5, 0.0};          // M' 点
static const double Kp[2]    = {1.0/3.0, -1.0/3.0}; // K' 点

static Segment path[] = {
    {M, GAMMA, 40},
    {GAMMA, K, 40},
    {K, Mp, 40},
    {Mp, Kp, 40}
};
static const int num_segments = sizeof(path)/sizeof(Segment);

int main(void)
{
    init_reciprocal();
    init_symmetry();
    init_vectors();

    Mat NN_mats[6];
    Mat NNN_mats[6];
    Mat TNN_mats[6];
    for (size_t i = 0; i < 6; ++i) {
        NN_mats[i] = mat_alloc(3, 3);
        NNN_mats[i] = mat_alloc(3, 3);
        TNN_mats[i] = mat_alloc(3, 3);
    }
#if NN_MODEL
    generate_hopping_matrices_NN(NN_mats);
#if NNN_MODEL
    generate_hopping_matrices_NNN(NNN_mats);
    generate_hopping_matrices_TNN(TNN_mats);
#endif
#endif

    double total_len = 0.0;
    double seg_len[num_segments];
    for (int s = 0; s < num_segments; s++) {
        double dk1 = path[s].end[0] - path[s].start[0];
        double dk2 = path[s].end[1] - path[s].start[1];
        double dkx = dk1 * b1x + dk2 * b2x;
        double dky = dk1 * b1y + dk2 * b2y;
        seg_len[s] = sqrt(dkx*dkx + dky*dky);
        total_len += seg_len[s];
    }

#if SOC_MODEL
#if SOC_SPIN
    printf("# k-path    Up1    Up2    Up3    Dn1    Dn2    Dn3\n");
#else
    printf("# k-path (frac)    E1    E2    E3    E4    E5    E6\n");
#endif
#else
    printf("# k-path (fractional)    E1 (eV)    E2 (eV)    E3 (eV)\n");
#endif
    printf("# path: M → Γ → K → M' → K'\n");
    printf("# Symmetry points: M=0, Γ=%.6f, K=%.6f, M'=%.6f, K'=%.6f\n",
            (seg_len[0])/(2*PI), (seg_len[0]+seg_len[1])/(2*PI), (seg_len[0]+seg_len[1]+seg_len[2])/(2*PI), (total_len)/(2*PI));

    double cumulative = 0.0;
    double complex Hgamma[3][3];
    double kgx = 0.0, kgy = 0.0; // Γ点分数坐标
    build_Hk(kgx, kgy, NN_mats, NNN_mats, TNN_mats, Hgamma);
    printf("\n");
    printf("# H at Gamma (0,0):\n");
    for (int i=0;i<3;i++) {
        printf("# ");
        for (int j=0;j<3;j++) printf("(%f%+fi) ", creal(Hgamma[i][j]), cimag(Hgamma[i][j]));
        printf("\n");
    }

    double complex HK[3][3];
    double kKx = 2.0/3.0, kKy = 1.0/3.0;
    build_Hk(kKx, kKy, NN_mats, NNN_mats, TNN_mats, HK);
    printf("\n");
    printf("# H at K (2/3,1/3):\n");
    for (int i=0;i<3;i++) {
        printf("# ");
        for (int j=0;j<3;j++) printf("(%f%+fi) ", creal(HK[i][j]), cimag(HK[i][j]));
        printf("\n");
    }

#if SOC_MODEL
    #if SOC_SPIN
    {
        double kx, ky;
        frac_to_cart_k(GAMMA[0], GAMMA[1], &kx, &ky);
        double complex H_gamma[3][3];
        build_Hk(kx, ky, NN_mats, NNN_mats, TNN_mats, H_gamma);

        double eig_up[3], eig_dn[3];
        diagonalize_soc_spin_blocks(H_gamma, eig_up, eig_dn);

        // 合并排序 6 个值
        double eig_all[6];
        for (int i = 0; i < 3; i++) eig_all[i]   = eig_up[i];
        for (int i = 0; i < 3; i++) eig_all[i+3] = eig_dn[i];
        // 排序
        for (int i = 0; i < 5; i++)
            for (int j = i+1; j < 6; j++)
                if (eig_all[i] > eig_all[j]) {
                    double tmp = eig_all[i];
                    eig_all[i] = eig_all[j];
                    eig_all[j] = tmp;
                }

        double vbm = find_vbm_from_sorted_eigenvalues(eig_all, 6);
        printf("\n");
        printf("# TB_Gamma_VBM = %12.6f\n", vbm);
    }
    #else
    {
        double kx, ky;
        frac_to_cart_k(GAMMA[0], GAMMA[1], &kx, &ky);
        double complex H_gamma[3][3];
        build_Hk(kx, ky, NN_mats, NNN_mats, TNN_mats, H_gamma);

        double eig[6];
        diagonalize_6x6_with_soc(H_gamma, eig);  // 已排序升序

        double vbm = find_vbm_from_sorted_eigenvalues(eig, 6);
        printf("\n");
        printf("# TB_Gamma_VBM = %12.6f\n", vbm);
    }
    #endif
#else
    {
        double kx, ky;
        frac_to_cart_k(GAMMA[0], GAMMA[1], &kx, &ky);
        double complex H_gamma[3][3];
        build_Hk(kx, ky, NN_mats, NNN_mats, TNN_mats, H_gamma);
        double eig[3];
        diagonalize_3x3_hermitian(H_gamma, eig); // 已排序
        double vbm = find_vbm_from_sorted_eigenvalues(eig, 3);
        printf("\n");
        printf("# TB_Gamma_VBM = %12.6f\n", vbm);
    }
#endif
#if SOC_MODEL
    // 打印 K 点的 6 个 SOC 本征值，并计算 VBM 劈裂
    {
        double kxK, kyK;
        frac_to_cart_k(K[0], K[1], &kxK, &kyK);
        double complex HK_soc[3][3];
        build_Hk(kxK, kyK, NN_mats, NNN_mats, TNN_mats, HK_soc);

        double eig_K[6];
        diagonalize_6x6_with_soc(HK_soc, eig_K);

        printf("\n# SOC eigenvalues at K (2/3, 1/3):\n");
        for (int b = 0; b < 6; b++) {
            printf("#   E[%d] = %12.6f eV\n", b, eig_K[b]);
        }

        // double gap_threshold = 0.5; // 假设带隙 > 0.5 eV
        double vbm1 = -100.0, vbm2 = -100.0;
        for (int b = 0; b < 6; b++) {
            if (eig_K[b] < 0.5) { // 取低于阈值的为价带
                if (eig_K[b] > vbm1) {
                    vbm2 = vbm1;
                    vbm1 = eig_K[b];
                } else if (eig_K[b] > vbm2) {
                    vbm2 = eig_K[b];
                }
            }
        }
        printf("# VBM1 = %12.6f eV, VBM2 = %12.6f eV\n", vbm1, vbm2);
        printf("# VBM SOC splitting = %12.6f eV  (expected 2*lambda = %12.6f eV)\n",
               vbm1 - vbm2, 2.0 * lambda);
    }
#endif
    for (int s=0; s<num_segments; s++) {
        double start_x = path[s].start[0], start_y = path[s].start[1];
        double end_x   = path[s].end[0],   end_y   = path[s].end[1];
        int steps = path[s].steps;
        for (int i=0; i<=steps; i++) {
            double t = (double)i / steps;
            double fx = start_x + t * (end_x - start_x);
            double fy = start_y + t * (end_y - start_y);
            double kx, ky;
            frac_to_cart_k(fx, fy, &kx, &ky);

            double complex H[3][3];
            build_Hk(kx, ky, NN_mats, NNN_mats, TNN_mats, H);

#if SOC_MODEL
    #if SOC_SPIN
            double eig_up[3], eig_dn[3];
            diagonalize_soc_spin_blocks(H, eig_up, eig_dn);

            double d = t * seg_len[s];
            double xcoord = (cumulative + d) / (2*PI);
            // 输出 7 列：x, up1, up2, up3, dn1, dn2, dn3
            printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                    xcoord,
                    eig_up[0], eig_up[1], eig_up[2],
                    eig_dn[0], eig_dn[1], eig_dn[2]);
    #else
            double eig[6];
            diagonalize_6x6_with_soc(H, eig);

            double d = t * seg_len[s];
            double xcoord = (cumulative + d) / (2*PI);
            printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                    xcoord, eig[0], eig[1], eig[2], eig[3], eig[4], eig[5]);
    #endif
#else
            double eig[3];
            diagonalize_3x3_hermitian(H, eig);

            double d = t * seg_len[s];
            double xcoord = (cumulative + d) / (2*PI);
            printf("%12.6f %12.6f %12.6f %12.6f\n",
                    xcoord, eig[0], eig[1], eig[2]);
#endif

        }
        cumulative += seg_len[s];

    }

    return 0;
}
