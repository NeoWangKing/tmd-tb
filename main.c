#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <lapacke.h> 
#include "mat.h"

#define PI 3.14159265358979323846

static double b1x, b1y, b2x, b2y;

// 模型参数
static const double e1  =  0.683;
static const double e2  =  1.707;

// 最近邻参数 (NN)
static const double t11 = -0.146;
static const double t12 =  0.114;
static const double t13 =  0.506;
static const double t22 =  0.085;
static const double t23 =  0.162;
static const double t33 =  0.073;

// 次近邻参数 (NNN)
static const double r11 =  0.060;
static const double r12 = -0.236;
static const double r21 =  0.067;
static const double r22 =  0.016;
static const double r23 =  0.087;

// 第三近邻参数 (TNN)
static const double u11 = -0.038;
static const double u12 =  0.046;
static const double u13 =  0.001;
static const double u22 =  0.266;
static const double u23 = -0.176;
static const double u33 = -0.150;

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
    // 最近邻 (NN)
    double NN_coords[6][2] = {
        {1.0,0.0}, {0.5,sqrt(3.0)/2.0}, {-0.5,sqrt(3.0)/2.0},
        {-1.0,0.0}, {-0.5,-sqrt(3.0)/2.0}, {0.5,-sqrt(3.0)/2.0}
    };

    // 次近邻 (NNN)
    double NNN_coords[6][2] = {
        {1.5,sqrt(3.0)/2.0}, {0.0,sqrt(3.0)}, {-1.5,sqrt(3.0)/2.0},
        {-1.5,-sqrt(3.0)/2.0}, {0.0,-sqrt(3.0)}, {1.5,-sqrt(3.0)/2.0},
    };

    // 第三近邻 (TNN)
    double TNN_coords[6][2] = {
        {2.0,0.0}, {1.0,sqrt(3.0)}, {-1.0,sqrt(3.0)},
        {-2.0,0.0}, {-1.0,-sqrt(3.0)}, {1.0,-sqrt(3.0)}
    };

    NN_vec  = mat_alloc(6, 2);
    NNN_vec = mat_alloc(6, 2);
    TNN_vec = mat_alloc(6, 2);

    for (size_t i = 0; i < 6; ++i) {
        MAT_AT(NN_vec,  i, 0) = NN_coords[i][0];
        MAT_AT(NN_vec,  i, 1) = NN_coords[i][1];

        MAT_AT(NNN_vec, i, 0) = NNN_coords[i][0];
        MAT_AT(NNN_vec, i, 1) = NNN_coords[i][1];

        MAT_AT(TNN_vec, i, 0) = TNN_coords[i][0];
        MAT_AT(TNN_vec, i, 1) = TNN_coords[i][1];
    }
}

// 对实数矩阵施加对称操作： M' = op * M * op^T
void apply_symmetry(Mat Mout, Mat op, Mat Min)
{
    Mat tmp = mat_alloc(3, 3); Mat opT = mat_alloc(3, 3);
    mat_transpose(opT, op);
    mat_dot(tmp, op, Min);
    mat_dot(Mout, tmp, opT);
}

// 最近邻 E(R1) 
static void build_E_NN_R1(Mat E)
{
    MAT_AT(E, 0, 0) =  t11;             MAT_AT(E, 0, 1) =  t12; MAT_AT(E, 0, 2) = t13;
    MAT_AT(E, 1, 0) = -t12;             MAT_AT(E, 1, 1) =  t22; MAT_AT(E, 1, 2) = t23;
    MAT_AT(E, 2, 0) =  t13;             MAT_AT(E, 2, 1) = -t23; MAT_AT(E, 2, 2) = t33;
}

// 次近邻 E(~R1)
static void build_E_NNN_R1(Mat E)
{
    MAT_AT(E, 0, 0) =  r11;             MAT_AT(E, 0, 1) = r12;  MAT_AT(E, 0, 2) = -r12/sqrt(3.0);
    MAT_AT(E, 1, 0) =  r21;             MAT_AT(E, 1, 1) = r22;  MAT_AT(E, 1, 2) = r23;
    MAT_AT(E, 2, 0) = -r21/sqrt(3.0); MAT_AT(E, 2, 1) = r23;  MAT_AT(E, 2, 2) = r22 + 2*r23/sqrt(3.0);
}

// 第三近邻 E(2R1)
static void build_E_TNN_R1(Mat E)
{
    MAT_AT(E, 0, 0) =  u11;             MAT_AT(E, 0, 1) =  u12; MAT_AT(E, 0, 2) = u13;
    MAT_AT(E, 1, 0) = -u12;             MAT_AT(E, 1, 1) =  u22; MAT_AT(E, 1, 2) = u23;
    MAT_AT(E, 2, 0) =  u13;             MAT_AT(E, 2, 1) = -u23; MAT_AT(E, 2, 2) = u33;
}

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

// k 空间哈密顿量构造
// 将 k 点分数坐标(f1,f2)转换为笛卡尔波矢
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
    
#if 1
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
#endif
    
#if 1
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
    generate_hopping_matrices_NN(NN_mats);
    generate_hopping_matrices_NNN(NNN_mats);
    generate_hopping_matrices_TNN(TNN_mats);

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

    printf("# k-path (fractional)    E1 (eV)    E2 (eV)    E3 (eV)\n");
    printf("# path: M → Γ → K → M' → K'\n");
    printf("# Symmetry points: M=0, Γ=%.6f, K=%.6f, M'=%.6f, K'=%.6f\n",
            (seg_len[0])/(2*PI), (seg_len[0]+seg_len[1])/(2*PI), (seg_len[0]+seg_len[1]+seg_len[2])/(2*PI), (total_len)/(2*PI));

    double cumulative = 0.0;
    double complex Hgamma[3][3];
    double kgx = 0.0, kgy = 0.0; // Γ点分数坐标
    build_Hk(kgx, kgy, NN_mats, NNN_mats, TNN_mats, Hgamma);
    printf("H at Gamma (0,0):\n");
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) printf("(%f%+fi) ", creal(Hgamma[i][j]), cimag(Hgamma[i][j]));
        printf("\n");
    }

    double complex HK[3][3];
    double kKx = 2.0/3.0, kKy = 1.0/3.0;
    build_Hk(kKx, kKy, NN_mats, NNN_mats, TNN_mats, HK);
    printf("H at K (2/3,1/3):\n");
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) printf("(%f%+fi) ", creal(HK[i][j]), cimag(HK[i][j]));
        printf("\n");
    }
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

            double eig[3];
            diagonalize_3x3_hermitian(H, eig);

            double d = t * seg_len[s];
            double xcoord = cumulative + d;
            xcoord /= 2*PI;
            printf("%12.6f %12.6f %12.6f %12.6f\n", xcoord, eig[0], eig[1], eig[2]);

        }
        cumulative += seg_len[s];

    }

    return 0;
}
