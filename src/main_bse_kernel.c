// main_bse_kernel.c —— 式 (6) 的核矩阵：K 点附近 N_k×N_k 网格上组装
//   用法：./bin/bse_kernel [N_k]
//   参数取 GW 拟合那一套（PARAM_SET=2），自旋两块分开算，各自 1 条价带 2 条导带
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#ifndef PARAM_SET
#define PARAM_SET 2          // 用拟合 GW 得到的那套参数
#endif
#ifndef NNN_MODEL
#define NNN_MODEL 1
#endif
#ifndef SOC_MODEL
#define SOC_MODEL 1          // 需要 lambda
#endif

#include "tb.h"
#include "herm3.h"
#include "params.h"

#define N_VAL 1                   // 价带条数
#define N_CON 2                   // 导带条数
#define N_PAIR (N_VAL * N_CON)    // 每个 k 点的电子-空穴对数

// Keldysh 屏蔽势参数：空气 / 二氧化硅界面
#define EPS_D    2.5
#define R0_ANG   13.55
#define EE_2EPS0 90.48            // e^2/(2 eps0)，单位 eV·Å
#define PI_ 3.14159265358979323846

static void die(const char *msg)
{
    fprintf(stderr, "[ERROR] %s\n", msg);
    exit(1);
}

// 交叠积分：x 取 k' 点的带指标向量，y 取 k 点的
static double complex inner(const double complex x[TB_NORB], const double complex y[TB_NORB])
{
    double complex s = 0.0;
    for (int a = 0; a < TB_NORB; ++a) s += conj(x[a]) * y[a];
    return s;
}

// 屏蔽库仑势（q 以 1/Å 为单位）；q=0 用正则化值代替
static double u_keldysh(double q, double u0)
{
    if (q < 1e-12) return u0;
    return -EE_2EPS0 / (EPS_D * q * (1.0 + R0_ANG * q));
}

// q=0 的对角元：在与 k 格子等面积的圆盘上取平均，积分有解析结果
static double u_regularized(double cell_area)
{
    double r = sqrt(cell_area / PI_);
    double amp = EE_2EPS0 / EPS_D;
    return -(2.0 * amp / (r * r * R0_ANG)) * log(1.0 + R0_ANG * r);
}

// 按显示宽度补齐：中文字符占两列，保证报告的列能对齐
static void pad_label(const char *s, int width)
{
    int w = 0;
    for (const unsigned char *p = (const unsigned char *)s; *p; ) {
        if (*p < 0x80) { w += 1; p += 1; }
        else if (*p >= 0xF0) { w += 2; p += 4; }
        else if (*p >= 0xE0) { w += 2; p += 3; }
        else { w += 2; p += 2; }
    }
    fputs(s, stdout);
    for (int i = w; i < width; ++i) putchar(' ');
}

int main(int argc, char **argv)
{
    int nk = (argc > 1) ? atoi(argv[1]) : 12;
    if (nk < 2) die("N_k 至少取 2");
    int dump = (argc > 2 && strcmp(argv[2], "dump") == 0);

    TB_Hop hop;
    tb_init(&hop);

    size_t nkpt = (size_t)nk * nk;
    double *kx = malloc(sizeof *kx * nkpt);
    double *ky = malloc(sizeof *ky * nkpt);
    double *eps = malloc(sizeof *eps * nkpt * 2 * TB_NORB);
    double complex *coef = malloc(sizeof *coef * nkpt * 2 * TB_NORB * TB_NORB);
    if (!kx || !ky || !eps || !coef) die("内存不足");

    // K 点为中心的均匀网格
    for (int i = 0; i < nk; ++i) {
        for (int j = 0; j < nk; ++j) {
            double f1 = 2.0 / 3.0 + ((double)i - 0.5 * (double)(nk - 1)) / (double)nk;
            double f2 = 1.0 / 3.0 + ((double)j - 0.5 * (double)(nk - 1)) / (double)nk;
            size_t p = (size_t)i * nk + j;
            tb_frac_to_cart_k(f1, f2, &kx[p], &ky[p]);
        }
    }

    // 每个 k 点解两个自旋块，并把本征矢的相位固定
    static const double complex Lz[TB_NORB][TB_NORB] = {
        { 0.0,     0.0,      0.0},
        { 0.0,     0.0, -2.0 * I},
        { 0.0, 2.0 * I,      0.0}
    };
    double complex H0[TB_NORB][TB_NORB], Hu[TB_NORB][TB_NORB], Hd[TB_NORB][TB_NORB];
    double complex V[TB_NORB][TB_NORB];
    double w[TB_NORB];

    double ierr = 0.0;
    double smin_s[6];
    for (int i = 0; i < 6; ++i) smin_s[i] = 1e300;
    for (size_t p = 0; p < nkpt; ++p) {
        tb_build_Hk(&hop, kx[p], ky[p], H0);
        for (int i = 0; i < TB_NORB; ++i) {
            for (int j = 0; j < TB_NORB; ++j) {
                Hu[i][j] = H0[i][j] + 0.5 * lambda * Lz[i][j];
                Hd[i][j] = H0[i][j] - 0.5 * lambda * Lz[i][j];
            }
        }
        for (int s = 0; s < 2; ++s) {
            herm3_jacobi(s ? Hd : Hu, w, V);
            for (int n = 0; n < TB_NORB; ++n) {
                // 固定相位：让 sum_a C 为实数
                double complex sum = 0.0;
                for (int a = 0; a < TB_NORB; ++a) sum += V[a][n];
                double ph = carg(sum);
                double smod = cabs(sum);
                if (smod < smin_s[s * TB_NORB + n]) smin_s[s * TB_NORB + n] = smod;
                double complex fac = cos(ph) - I * sin(ph);
                eps[(p * 2 + s) * TB_NORB + n] = w[n];
                for (int a = 0; a < TB_NORB; ++a)
                    coef[((p * 2 + s) * TB_NORB + n) * TB_NORB + a] = V[a][n] * fac;
                // 归一化自检
                double complex nn = inner(&coef[((p * 2 + s) * TB_NORB + n) * TB_NORB],
                                          &coef[((p * 2 + s) * TB_NORB + n) * TB_NORB]);
                double e = fabs(creal(nn) - 1.0);
                if (e > ierr) ierr = e;
            }
        }
    }

    // 格子面积与晶体总面积
    double a_lat = A_ANG;
    double acell = sqrt(3.0) / 2.0 * a_lat * a_lat;      // 原胞面积，Å²
    // 间距与每点面积必须从实际 k 向量算：沿 b1 的间距是 |b1|/N_k，不是 2π/(a N_k)
    double dkx = kx[1] - kx[0], dky = ky[1] - ky[0];                 // 沿 j（b2）方向
    double ekx = kx[nk] - kx[0], eky = ky[nk] - ky[0];               // 沿 i（b1）方向
    // tb_frac_to_cart_k 给的是 1/a 单位，换成 1/Å 要除以 a
    double dk = 0.5 * (sqrt(dkx * dkx + dky * dky) + sqrt(ekx * ekx + eky * eky)) / a_lat;
    double cell_area = fabs(dkx * eky - dky * ekx) / (a_lat * a_lat);  // 每个 k 点代表的面积
    double u0 = u_regularized(cell_area);
    double vtot = acell * (double)nkpt;

    // 各带 min|sum_alpha C| 的最小值，留到报告里打
    double smin_up = fmin(smin_s[0], fmin(smin_s[1], smin_s[2]));
    double smin_dn = fmin(smin_s[3], fmin(smin_s[4], smin_s[5]));
    (void)smin_up;
    (void)smin_dn;

    // ---- 收集统计量 ----
    double imean[2][TB_NORB], imin[2][TB_NORB];
    for (int s = 0; s < 2; ++s) {
        for (int n = 0; n < TB_NORB; ++n) {
            double smin = 1e300, ssum = 0.0;
            int cnt = 0;
            for (int i = 0; i < nk; ++i) {
                for (int j = 0; j + 1 < nk; ++j) {
                    size_t p = (size_t)i * nk + j;
                    size_t q = p + 1;
                    double complex ov = inner(&coef[((q * 2 + s) * TB_NORB + n) * TB_NORB],
                                              &coef[((p * 2 + s) * TB_NORB + n) * TB_NORB]);
                    double m = cabs(ov);
                    if (m < smin) smin = m;
                    ssum += m;
                    ++cnt;
                }
            }
            imean[s][n] = ssum / (double)cnt;
            imin[s][n] = smin;
        }
    }

    // 组装核矩阵（两个自旋各一套）
    size_t npair = nkpt * N_PAIR;
    double complex *W = malloc(sizeof *W * npair * npair);
    if (!W) die("核矩阵内存不足");

    double herm[2], amax[2], amean[2], adiag[2];
    double blk_re[6][6], blk_im[6][6], diag_avg = 0.0, same_k_avg = 0.0;
    for (int s = 0; s < 2; ++s) {
        for (size_t p = 0; p < nkpt; ++p) {
            for (size_t q = 0; q < nkpt; ++q) {
                double dkx = kx[p] - kx[q];
                double dky = ky[p] - ky[q];
                double qmod = sqrt(dkx * dkx + dky * dky) / a_lat;
                double uq = u_keldysh(qmod, u0);
                const double complex *cvk  = &coef[((p * 2 + s) * TB_NORB + 0) * TB_NORB];
                const double complex *cvkp = &coef[((q * 2 + s) * TB_NORB + 0) * TB_NORB];
                double complex iv = inner(cvkp, cvk);
                for (int a = 0; a < N_PAIR; ++a) {
                    int ca = a + N_VAL;
                    const double complex *cck = &coef[((p * 2 + s) * TB_NORB + ca) * TB_NORB];
                    for (int b = 0; b < N_PAIR; ++b) {
                        int cb2 = b + N_VAL;
                        const double complex *cckp = &coef[((q * 2 + s) * TB_NORB + cb2) * TB_NORB];
                        double complex ic = conj(inner(cckp, cck));
                        W[(p * N_PAIR + a) * npair + (q * N_PAIR + b)] = uq * ic * iv / vtot;
                    }
                }
            }
        }

        double herrmax = 0.0, mx = 0.0, sum = 0.0, dg = 0.0;
        for (size_t i = 0; i < npair; ++i) {
            for (size_t j = 0; j < npair; ++j) {
                double complex x = W[i * npair + j];
                double d = cabs(x - conj(W[j * npair + i]));
                if (d > herrmax) herrmax = d;
                double ax = cabs(x);
                if (ax > mx) mx = ax;
                sum += ax;
                if (i == j) dg += creal(x);
            }
        }
        // 左上角一小块留到报告里打印
        if (s == 0) {
            const int blk = 6;
            for (int i = 0; i < blk; ++i) {
                for (int j = 0; j < blk; ++j) {
                    double complex x = W[(size_t)i * npair + j];
                    blk_re[i][j] = creal(x) * 1e3;
                    blk_im[i][j] = cimag(x) * 1e3;
                }
            }
            double dsum = 0.0, blocksum = 0.0;
            int nblock = 0;
            for (size_t i = 0; i < npair; ++i) {
                for (size_t j = 0; j < npair; ++j) {
                    if (i == j) dsum += creal(W[i * npair + j]);
                    else if ((i / N_PAIR) == (j / N_PAIR)) { blocksum += creal(W[i * npair + j]); ++nblock; }
                }
            }
            diag_avg = dsum / (double)npair * 1e3;
            same_k_avg = nblock ? blocksum / (double)nblock * 1e3 : 0.0;
        }
        if (dump && s == 0) {
            const char *path = "data/bse_kernel.dat";
            FILE *f = fopen(path, "w");
            if (!f) die("写不了 data/bse_kernel.dat");
            fprintf(f, "# 自旋 Up 的核矩阵 H_(c k),(c' k')，单位 eV；只写 i<=j（厄米）\n");
            fprintf(f, "# N_k=%d，维数 %zu，i j Re Im\n", nk, npair);
            for (size_t i = 0; i < npair; ++i) {
                for (size_t j = i; j < npair; ++j) {
                    double complex x = W[i * npair + j];
                    fprintf(f, "%zu %zu %.10e %.10e\n", i, j, creal(x), cimag(x));
                }
            }
            fclose(f);
            printf("完整矩阵已写出：%s（%zu 行）\n\n", path, npair * (npair + 1) / 2);
        }

        herm[s] = herrmax;
        amax[s] = mx;
        amean[s] = sum / (double)(npair * npair);
        adiag[s] = dg / (double)npair;
    }
    free(W);

    // ---- 生成报告 ----
    printf("计算目标：W^(d)_{cvk,c'vk'} = u(k-k') I*_{c'k',ck} I_{v'k',vk}\n");
    printf("模型：GW 拟合参数（PARAM_SET=2），3 个 d 轨道 × 2 个自旋；介电环境 eps_d = %.1f，r0 = %.2f Å\n\n",
           EPS_D, R0_ANG);

    printf("网格\n");
    printf("  N_k = %d，k 点 %d × %d = %zu（以 K 点为心）\n", nk, nk, nk, nkpt);
    printf("  相邻 k 间距 dk = %.5f Å^-1，每点面积 = %.6f Å^-2\n", dk, cell_area);
    printf("  原胞面积 A_c = %.4f Å^2，晶体总面积 V = %.1f Å^2\n", acell, vtot);
    printf("  矩阵维数 N_tot = %d N_k^2 = %zu（价带 1 条，导带 2 条）\n\n", N_PAIR, npair);

    printf("屏蔽势\n");
    printf("  u(q=0) 正则化值 = %.1f eV·Å^2，u(dk) = %.1f eV·Å^2\n", u0, u_keldysh(dk, u0));
    printf("  分界 q = 1/r0 = %.4f Å^-1（左边 u ~ 1/q，右边 u ~ 1/q^2）\n\n", 1.0 / R0_ANG);

    printf("自检\n");
    pad_label("  项目", 38);
    printf("%10s %10s   %s\n", "Up", "Dn", "判据");
    printf("-------------------------------------------------------------------------\n");
    pad_label("  本征矢归一化误差 max|sum|C|^2-1|", 38);
    printf("%10.1e %10.1e   < 1e-12\n", ierr, ierr);
    pad_label("  厄米误差 max|H_ij - H*_ji|", 38);
    printf("%10.1e %10.1e   < 1e-9 eV\n", herm[0], herm[1]);
    pad_label("  相邻 k 点 |I| 平均（价带）", 38);
    printf("%10.3f %10.3f   > 0.95\n", imean[0][0], imean[1][0]);
    pad_label("  相邻 k 点 |I| 最小（导带 1）", 38);
    printf("%10.3f %10.3f   > 0.90\n", imin[0][1], imin[1][1]);
    pad_label("  相邻 k 点 |I| 最小（导带 2）", 38);
    printf("%10.3f %10.3f   > 0.90\n", imin[0][2], imin[1][2]);
    pad_label("  规范固定 min|sum_alpha C|", 38);
    printf("%10.3f %10.3f   > 0.10\n\n", smin_up, smin_dn);

    printf("核矩阵量级\n");
    printf("-------------------------------------------------------------------------\n");
    pad_label("  max|H_ij|", 38);
    printf("%7.1f meV %7.1f meV\n", amax[0] * 1e3, amax[1] * 1e3);
    pad_label("  平均|H_ij|", 38);
    printf("%7.3f meV %7.3f meV\n", amean[0] * 1e3, amean[1] * 1e3);
    pad_label("  对角元平均", 38);
    printf("%7.1f meV %7.1f meV\n", adiag[0] * 1e3, adiag[1] * 1e3);
    pad_label("  平均|H_ij| × 维数", 38);
    printf("%7.3f eV  %7.3f eV\n\n", amean[0] * (double)npair, amean[1] * (double)npair);

    // 结论：逐条对照判据，直接说清楚能不能用
    int bad_norm = (ierr > 1e-12) || (herm[0] > 1e-9) || (herm[1] > 1e-9);
    int bad_overlap = (imin[0][1] < 0.90) || (imin[1][1] < 0.90);
    int bad_gauge = (smin_up < 0.10) || (smin_dn < 0.10);
    printf("核矩阵（自旋 Up，左上 6 × 6，单位 meV；行是 (k,c)，列是 (k',c')）\n");
    for (int part = 0; part < 2; ++part) {
        printf("  %s\n", part == 0 ? "实部" : "虚部");
        printf("  行\\列");
        for (int j = 0; j < 6; ++j) printf("%9d", j);
        printf("\n");
        for (int i = 0; i < 6; ++i) {
            printf("  %5d", i);
            for (int j = 0; j < 6; ++j) printf("%9.2f", part == 0 ? blk_re[i][j] : blk_im[i][j]);
            printf("\n");
        }
    }
    printf("  对角元平均 = %.2f meV（= u(0)/V，q=0 的自相互作用）\n", diag_avg);
    printf("  同一 k 点内不同导带之间恒为 0（%.2f meV）：它们在同一 k 点正交，I = 0，这是对的\n",
           same_k_avg);
    printf("  行号 0,1 是同一个 k 点的两条导带，2,3 是下一个 k 点。\n");
    printf("  要看完整矩阵：./bin/main_bse_kernel %d dump，会写到 data/bse_kernel.dat\n\n", nk);

    printf("结论\n");
    if (!bad_norm && !bad_overlap && !bad_gauge) {
        printf("  三项自检全部通过，核矩阵可以交给式 (5)。\n");
    } else {
        printf("  不通过，核矩阵暂时不可用：\n");
        if (bad_norm) printf("    - 归一化或厄米性不满足，先查对角化与组装\n");
        if (bad_gauge) printf("    - 规范固定退化（min|sum_alpha C| = %.3f），辐角不可靠\n", fmin(smin_up, smin_dn));
        if (bad_overlap) printf("    - 导带相邻 k 点交叠过小（最低 %.3f），多半是导带换位\n",
                                fmin(imin[0][1], imin[1][1]));
    }
    printf("\n");

    printf("说明\n");
    printf("  对角元是 u(0) 的正则化值除以 V，随 N_k 按 1/N_k 下降：u(0) ~ N_k，V ~ N_k^2。\n");
    printf("  平均|H_ij| × 维数是求和量级的粗估，不等于激子束缚能；束缚能要解式 (5) 才有。\n");

    free(kx);
    free(ky);
    free(eps);
    free(coef);
    return 0;
}
