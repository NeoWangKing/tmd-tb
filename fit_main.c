// 拟合诊断（Step A）：把 GW 参考能带与当前 TB 能带逐点比较，
// 把残差分解成"刚性成分"（整体平移 / 带隙张开）与"形状成分"（色散、曲率）。
//
// 读：data/band_gw.dat                   —— TB 在 GW 的 285 个 k 点上的能带
//     MoS2-GWBSE-data/wannier90_band.dat —— GW 参考能带（17 条带，按块存放）
// 写：data/fit_residual.dat              —— 残差曲线，供 gnuplot/plot_fit.gp 使用

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "tb.h"

#define TB_NBAND 3
#define GW_NBAND_MAX 64
#define NK_MAX 100000

// 比较的带对：TB 第 TB_BAND[i] 条  <->  GW 第 GW_BAND[i] 条
static const int TB_BAND[3] = {1, 2, 3};
static const int GW_BAND[3] = {9, 10, 11};
static const char *PAIR[3]  = {"价带  VB ", "导带1 CB1", "导带2 CB2"};

// GW 路径上的高对称点位置 (Å^-1)，取自 wannier90_band.labelinfo.dat
static const double SYM_K[4]  = {0.0, 1.31633, 1.97450, 3.11448};
static const char  *SYM_NM[4] = {"Γ", "K", "M", "Γ(末)"};

typedef struct {
    int nband, nk;
    double *k;      // [nk]
    double *E;      // [nband*nk]，按带分块
} GwData;

static int line_is_blank(const char *s)
{
    for (; *s; ++s) if (*s != ' ' && *s != '\t' && *s != '\n' && *s != '\r') return 0;
    return 1;
}

static int gw_read(const char *path, GwData *g)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "打不开 %s\n", path); return -1; }

    // 第一遍：数每条带的点数，并检查各带点数一致
    char line[512];
    int nband = 0, nk = -1, cnt = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) {
            if (cnt > 0) {
                if (nk < 0) nk = cnt; else if (cnt != nk) return -2;
                nband++; cnt = 0;
            }
            continue;
        }
        if (line[0] == '#') continue;
        cnt++;
    }
    if (cnt > 0) { if (nk < 0) nk = cnt; else if (cnt != nk) return -2; nband++; }
    fclose(f);
    if (nband <= 0 || nk <= 0 || nband > GW_NBAND_MAX || nk > NK_MAX) return -3;

    g->nband = nband; g->nk = nk;
    g->k = malloc(sizeof *g->k * nk);
    g->E = malloc(sizeof *g->E * nk * nband);
    if (!g->k || !g->E) return -4;

    // 第二遍：读入；k 轴以第一块为准，其余块必须一致
    f = fopen(path, "r");
    if (!f) return -1;
    int band = 0, i = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) { if (i > 0) { band++; i = 0; } continue; }
        if (line[0] == '#') continue;
        double kv, ev;
        if (sscanf(line, "%lf %lf", &kv, &ev) != 2) continue;
        if (band >= nband || i >= nk) break;
        if (band == 0) g->k[i] = kv;
        else if (fabs(kv - g->k[i]) > 1e-8) { fclose(f); return -5; }
        g->E[band * nk + i] = ev;
        i++;
    }
    fclose(f);
    return 0;
}

// 读 TB 能带表：跳过 # 注释行，每行 1 + nband 列（k 与 nband 个能量）
static int tb_read(const char *path, int nband, int *nk_out, double **k_out, double **e_out)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "打不开 %s\n", path); return -1; }

    int cap = 1024, n = 0;
    double *k = malloc(sizeof *k * cap);
    double *e = malloc(sizeof *e * cap * nband);
    if (!k || !e) { fclose(f); return -2; }

    char line[512];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line_is_blank(line)) continue;
        double v[8];
        int got = 0;
        char *p = line;
        for (int j = 0; j < nband + 1; ++j) {
            char *end;
            double x = strtod(p, &end);
            if (end == p) break;
            v[got++] = x;
            p = end;
        }
        if (got != nband + 1) continue;         // 列数不对，跳过
        if (n == cap) {
            cap *= 2;
            double *nk2 = realloc(k, sizeof *k * cap);
            double *ne2 = realloc(e, sizeof *e * cap * nband);
            if (!nk2 || !ne2) { fclose(f); return -2; }
            k = nk2; e = ne2;
        }
        k[n] = v[0];
        for (int j = 0; j < nband; ++j) e[n * nband + j] = v[j + 1];
        n++;
    }
    fclose(f);
    *nk_out = n; *k_out = k; *e_out = e;
    return 0;
}

// 去掉相邻重复的 k（段端点被打印两次；GW 文件里每个 k 只出现一次）
static int dedup(int n, const double *k, const double *e, int nband, double *k2, double *e2)
{
    int m = 0;
    for (int i = 0; i < n; ++i) {
        if (m > 0 && fabs(k[i] - k2[m-1]) < 1e-9) continue;
        k2[m] = k[i];
        for (int j = 0; j < nband; ++j) e2[m * nband + j] = e[i * nband + j];
        m++;
    }
    return m;
}

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
    if (rc) { fprintf(stderr, "读 GW 能带失败 (rc=%d)\n", rc); return 1; }

    int tb_n, nk, nk2;
    double *tb_k, *tb_e;
    if (tb_read("data/band_gw.dat", TB_NBAND, &tb_n, &tb_k, &tb_e)) {
        fprintf(stderr, "读 TB 能带失败（先跑 ./build.sh）\n");
        return 1;
    }
    double *k  = malloc(sizeof *k  * tb_n);
    double *E  = malloc(sizeof *E  * tb_n * TB_NBAND);
    nk = dedup(tb_n, tb_k, tb_e, TB_NBAND, k, E);
    nk2 = gw.nk;
    if (nk != nk2) {
        fprintf(stderr, "k 点数不一致：TB %d，GW %d\n", nk, nk2);
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
        const double *gwv = gw.E + (GW_BAND[b] - 1) * gw.nk;
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
    if (!out) { fprintf(stderr, "写 data/fit_residual.dat 失败\n"); return 1; }
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
