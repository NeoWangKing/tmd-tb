// 拟合诊断：GW 参考能带 vs 当前 TB 的残差分解（刚性成分 / 形状成分）
// 读 data/band_gw.dat 与 MoS2-GWBSE-data/wannier90_band.dat，写 data/fit_residual.dat

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "gwdata.h"
#include "tb.h"

#define TB_NBAND 3

// 比较的带对：TB 第 TB_BAND[i] 条  <->  GW 第 GW_BAND[i] 条
static const int TB_BAND[3] = {1, 2, 3};
static const int GW_BAND[3] = {9, 10, 11};

// 候选能级范围（用于重连）：覆盖最低的几条导带
#define CAND_FIRST 10
#define CAND_LAST  14
#define NCAND      (CAND_LAST - CAND_FIRST + 1)
static const char *PAIR[3]  = {"价带  VB ", "导带1 CB1", "导带2 CB2"};

// GW 路径上的高对称点位置 (Å^-1)，取自 wannier90_band.labelinfo.dat
static const double SYM_K[4]  = {0.0, 1.31633, 1.97450, 3.11448};
static const char  *SYM_NM[4] = {"Γ", "K", "M", "Γ(末)"};


// RMS（相对 0）
static double rms0(const double *v, int n)
{
    double s = 0;
    for (int i = 0; i < n; ++i) s += v[i] * v[i];
    return sqrt(s / n);
}

// 均值 / 去均值后的 RMS / 最大绝对偏差
static void stats(const double *v, int n, double *mean, double *rms, double *maxabs)
{
    double m = 0;
    for (int i = 0; i < n; ++i) m += v[i];
    m /= n;
    double s = 0, mx = 0;
    for (int i = 0; i < n; ++i) {
        double d = v[i] - m;
        s += d * d;
        if (fabs(d) > mx) mx = fabs(d);
    }
    *mean = m; *rms = sqrt(s / n); *maxabs = mx;
}

int main(void)
{
    GwData gw;
    int rc = gw_read("MoS2-GWBSE-data/wannier90_band.dat", &gw);
    if (rc) { fprintf(stderr, "[ERROR] 读 GW 能带失败\n"); return 1; }

    int tb_n, nk, nk2;
    double *tb_k, *tb_e;
    if (tb_read("data/band_gw.dat", TB_NBAND, &tb_n, &tb_k, &tb_e)) {
        fprintf(stderr, "[ERROR] 读 TB 能带失败，先跑 ./build.sh 的 gw 目标\n");
        return 1;
    }
    // 重连：wannier90 的能带按能量排序，交叉处排序曲线会互换身份。
    // 价带（第 9 条）本身与上下都分得开，直接用；两条导带用重连后的物理能带。
    static double *gwref[3];
    {
        double *cand = malloc(sizeof(double) * NCAND * gw.nk);
        double *rec  = malloc(sizeof(double) * NCAND * gw.nk);
        for (int c = 0; c < NCAND; ++c)
            memcpy(cand + c * gw.nk, gw.E + (CAND_FIRST - 1 + c) * gw.nk, sizeof(double) * gw.nk);

        int anchor = 0;
        for (int i = 1; i < gw.nk; ++i)
            if (fabs(gw.k[i] - 1.31633) < fabs(gw.k[anchor] - 1.31633)) anchor = i;
        if (band_track(gw.nk, NCAND, cand, anchor, rec)) { fprintf(stderr, "[ERROR] 重连失败\n"); return 1; }

        double r0 = 0, r1 = 0;
        for (int c = 0; c < NCAND; ++c)
            for (int i = 1; i < gw.nk - 1; ++i) {
                double d0 = cand[c*gw.nk+i+1] - 2*cand[c*gw.nk+i] + cand[c*gw.nk+i-1];
                double d1 = rec[c*gw.nk+i+1]  - 2*rec[c*gw.nk+i]  + rec[c*gw.nk+i-1];
                r0 += d0*d0; r1 += d1*d1;
            }
        printf("GW 参考能带重连（斜率外推，锚点 = K）：粗糙度 %.5f -> %.5f\n", r0, r1);

        gwref[0] = gw.E + (GW_BAND[0] - 1) * gw.nk;      // 价带：本身就是干净的
        gwref[1] = rec + 0 * gw.nk;                      // 物理导带 1
        gwref[2] = rec + 1 * gw.nk;                      // 物理导带 2
    }

    double *k  = malloc(sizeof *k  * tb_n);
    double *E  = malloc(sizeof *E  * tb_n * TB_NBAND);
    nk = dedup(tb_n, tb_k, tb_e, TB_NBAND, k, E);
    nk2 = gw.nk;
    if (nk != nk2) {
        fprintf(stderr, "[ERROR] k 点数不一致：TB %d，GW %d\n", nk, nk2);
        return 1;
    }
    double dk = 0;
    for (int i = 0; i < nk; ++i) {
        double d = fabs(k[i] - gw.k[i]);
        if (d > dk) dk = d;
    }

    printf("拟合诊断（Step A）：GW 参考 vs 当前 TB\n");
    printf("  TB : data/band_gw.dat                       %d 点 -> 去重 %d 点\n", tb_n, nk);
    printf("  GW : MoS2-GWBSE-data/wannier90_band.dat     %d 带 x %d 点\n", gw.nband, gw.nk);
    printf("  k 轴最大偏差 = %.2e Å^-1\n\n", dk);

    double *R  = malloc(sizeof *R  * nk * 3);   // 原始残差
    double *RS = malloc(sizeof *RS * nk * 3);   // 刚性平移后
    double *RL = malloc(sizeof *RL * nk * 3);   // 平移 + 线性伸缩后
    double shift[3], scale[3], raw_rms[3], sh_rms[3], sh_max[3], lin_rms[3];
    int    sh_imax[3];

    printf("带对         原始RMS  平移后RMS  max|平移后|   平移量Δ   平移+伸缩RMS   伸缩a\n");
    printf("------------------------------------------------------------------------------\n");
    for (int b = 0; b < 3; ++b) {
        const double *tb  = E + (TB_BAND[b] - 1);
        const double *gwv = gwref[b];
        double *rb = R + b * nk, *rsb = RS + b * nk, *rlb = RL + b * nk;

        double mx = 0, my = 0;
        for (int i = 0; i < nk; ++i) {
            double x = tb[i * TB_NBAND], y = gwv[i];
            rb[i] = y - x;
            mx += x; my += y;
        }
        mx /= nk; my /= nk;

        double mean, rms, mx1;
        stats(rb, nk, &mean, &rms, &mx1);      // 均值即最优刚性平移量
        shift[b] = mean;
        for (int i = 0; i < nk; ++i) rsb[i] = rb[i] - mean;
        stats(rsb, nk, &mean, &rms, &mx1);

        // 平移 + 线性伸缩（最小二乘闭式解）：y ≈ a*(x - <x>) + <y>
        double sxx = 0, sxy = 0;
        for (int i = 0; i < nk; ++i) {
            double dx = tb[i * TB_NBAND] - mx;
            sxx += dx * dx;
            sxy += dx * (gwv[i] - my);
        }
        scale[b] = sxy / sxx;
        for (int i = 0; i < nk; ++i) {
            double dx = tb[i * TB_NBAND] - mx;
            rlb[i] = (gwv[i] - my) - scale[b] * dx;
        }
        stats(rlb, nk, &mean, &rms, &mx1);

        raw_rms[b] = rms0(rb, nk);
        sh_rms[b]  = rms0(rsb, nk);
        sh_max[b]  = mx1;
        lin_rms[b] = rms0(rlb, nk);
        sh_imax[b] = 0;
        for (int i = 1; i < nk; ++i)
            if (fabs(rsb[i]) > fabs(rsb[sh_imax[b]])) sh_imax[b] = i;

        printf("%s   %7.3f   %7.3f     %7.3f     %+8.3f      %7.3f     %6.3f\n",
               PAIR[b], raw_rms[b], sh_rms[b], sh_max[b], shift[b], lin_rms[b], scale[b]);
    }
    printf("\n  单位 eV。'平移后RMS' = 只做最优刚性平移后剩下的形状误差，是后续拟合要打到的下界；\n");
    printf("  '伸缩a' = 线性拟合斜率，a > 1 表示 GW 色散比 TB 宽。a 接近 1 而 RMS 仍大说明\n");
    printf("  误差是曲率/局部结构而非整体展宽。\n");
    printf("\n平移后残差最大的位置：\n");
    for (int b = 0; b < 3; ++b) {
        int i = sh_imax[b], s = 0;
        for (int j = 1; j < 4; ++j)
            if (fabs(SYM_K[j] - k[i]) < fabs(SYM_K[s] - k[i])) s = j;
        printf("  %s  |r| = %.3f eV 在 k = %.4f Å^-1（最近的高对称点 %s，距 %.4f）\n",
               PAIR[b], sh_max[b], k[i], SYM_NM[s], fabs(k[i] - SYM_K[s]));
    }

    printf("\n带隙：GW 需要比 TB 张开 %.3f eV（导带1 平移 %.3f − 价带平移 %.3f）\n",
           shift[1] - shift[0], shift[1], shift[0]);

    printf("\n高对称点残差 (eV)：r 为原始，r* 为各自最优刚性平移后\n");
    printf("%-6s %8s %8s %8s  | %8s %8s %8s\n", "点", "r_VB", "r_CB1", "r_CB2",
           "r_VB*", "r_CB1*", "r_CB2*");
    for (int s = 0; s < 4; ++s) {
        int i = 0;
        for (int j = 0; j < nk; ++j) if (fabs(k[j] - SYM_K[s]) < fabs(k[i] - SYM_K[s])) i = j;
        printf("%-6s %8.3f %8.3f %8.3f  | %8.3f %8.3f %8.3f\n", SYM_NM[s],
               R[0*nk+i], R[1*nk+i], R[2*nk+i],
               RS[0*nk+i], RS[1*nk+i], RS[2*nk+i]);
    }

    FILE *out = fopen("data/fit_residual.dat", "w");
    if (!out) { fprintf(stderr, "[ERROR] cannot write data/fit_residual.dat\n"); return 1; }
    fprintf(out, "# col1 : k (Å^-1)\n");
    fprintf(out, "# col2-4 : 原始残差 r = E_GW - E_TB           [VB, CB1, CB2]\n");
    fprintf(out, "# col5-7 : 各自减去最优刚性平移后的残差\n");
    fprintf(out, "# col8-10: 平移 + 线性伸缩后的残差\n");
    for (int i = 0; i < nk; ++i) {
        fprintf(out, "%10.6f", k[i]);
        for (int b = 0; b < 3; ++b) fprintf(out, " %10.5f", R[b * nk + i]);
        for (int b = 0; b < 3; ++b) fprintf(out, " %10.5f", RS[b * nk + i]);
        for (int b = 0; b < 3; ++b) fprintf(out, " %10.5f", RL[b * nk + i]);
        fprintf(out, "\n");
    }
    fclose(out);
    printf("\n已写出 data/fit_residual.dat\n");
    return 0;
}
