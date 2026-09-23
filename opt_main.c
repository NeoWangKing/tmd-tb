// 全路径最小二乘拟合（最朴素的版本：逐参数扫描 + 步长减半）
//
//   min_p  Σ_{n,k} [E_TB(p; n, k) − E_GW(n, k)]²
//   n = 价带（TB 第 1 条 ↔ GW 第 9 条）、第一导带（TB 第 2 条 ↔ GW 第 10 条）
//   k = GW 参考数据里的全部 285 个点
//
// 读 MoS2-GWBSE-data/wannier90_band.kpt（k 点）与 wannier90_band.dat（参考能量）
// 写 data/gw_fit_bands.dat（初值/拟合后的 TB 能带，供 plot_opt.gp 画图）

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "gwdata.h"
#include "tb.h"
#include "herm3.h"

#define NP 19            // e1,e2 + 6 个 t + 5 个 r + 6 个 u

#define GW_VB 9          // GW 参考能带里价带的序号（1-based）
#define GW_CB 10         // 第一导带

// GW 路径上的高对称点位置 (Å^-1)，用于报告分区残差
static const double SYM_K[3]  = {0.0, 1.31633, 1.97450};
static const char  *SYM_NM[3] = {"Γ", "K", "M"};

static int nk;
static double *kx, *ky;        // TB 单位（1/a）的笛卡尔 k
static double *gw_vb, *gw_cb;  // 参考能量

static void unpack(const double *p, TB_Params *q)
{
    q->e1 = p[0];  q->e2 = p[1];
    q->t11 = p[2]; q->t12 = p[3]; q->t13 = p[4];
    q->t22 = p[5]; q->t23 = p[6]; q->t33 = p[7];
    q->r11 = p[8]; q->r12 = p[9]; q->r21 = p[10];
    q->r22 = p[11]; q->r23 = p[12];
    q->u11 = p[13]; q->u12 = p[14]; q->u13 = p[15];
    q->u22 = p[16]; q->u23 = p[17]; q->u33 = p[18];
}

static void pack(const TB_Params *q, double *p)
{
    p[0] = q->e1;  p[1] = q->e2;
    p[2] = q->t11; p[3] = q->t12; p[4] = q->t13;
    p[5] = q->t22; p[6] = q->t23; p[7] = q->t33;
    p[8] = q->r11; p[9] = q->r12; p[10] = q->r21;
    p[11] = q->r22; p[12] = q->r23;
    p[13] = q->u11; p[14] = q->u12; p[15] = q->u13;
    p[16] = q->u22; p[17] = q->u23; p[18] = q->u33;
}

// 目标函数：两条带、全部 k 点的均匀加权最小二乘
static double loss(const double *p, TB_Hop *hop)
{
    TB_Params q;
    unpack(p, &q);
    tb_set_params(&q);
    tb_update_hopping(hop);

    double s = 0;
    for (int i = 0; i < nk; ++i) {
        double complex H[3][3];
        double e[3];
        tb_build_Hk(hop, kx[i], ky[i], H);
        herm3_eigvals(H, e);
        double dv = e[0] - gw_vb[i];
        double dc = e[1] - gw_cb[i];
        s += dv*dv + dc*dc;
    }
    return s;
}

// 把当前参数下的能带算出来（用于输出/画图）
static void eval_bands(const double *p, TB_Hop *hop, double *vb, double *cb)
{
    TB_Params q;
    unpack(p, &q);
    tb_set_params(&q);
    tb_update_hopping(hop);

    for (int i = 0; i < nk; ++i) {
        double complex H[3][3];
        double e[3];
        tb_build_Hk(hop, kx[i], ky[i], H);
        herm3_eigvals(H, e);
        vb[i] = e[0];
        cb[i] = e[1];
    }
}

static double rms(const double *a, const double *b)
{
    double s = 0;
    for (int i = 0; i < nk; ++i) s += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(s / nk);
}

int main(void)
{
    GwData gw;
    if (gw_read("MoS2-GWBSE-data/wannier90_band.dat", &gw)) return 1;
    if (gw.nband < GW_CB) { fprintf(stderr, "[ERROR] GW 数据带数不够\n"); return 1; }

    if (kpt_read("MoS2-GWBSE-data/wannier90_band.kpt", &nk, &kx, &ky)) return 1;
    if (nk != gw.nk) {
        fprintf(stderr, "[ERROR] k 点数不一致：kpt %d，band %d\n", nk, gw.nk);
        return 1;
    }
    gw_vb = gw.E + (GW_VB - 1) * nk;
    gw_cb = gw.E + (GW_CB - 1) * nk;

    // 自检：把 .kpt 重建的 k 点累加成路径长度，应与参考文件的 k 轴一致
    double cum = 0, dk = 0;
    for (int i = 0; i < nk; ++i) {
        if (i > 0) cum += hypot(kx[i] - kx[i-1], ky[i] - ky[i-1]);
        double d = fabs(cum / A_ANG - gw.k[i]);
        if (d > dk) dk = d;
    }
    printf("拟合：GW 参考 vs 三能带模型（全路径，%d 个 k 点，2 条带）\n", nk);
    printf("  k 点自检：由 .kpt 重建的路径长度与参考文件最大差 %.2e Å^-1\n", dk);
    if (dk > 1e-4) { fprintf(stderr, "[ERROR] k 点重建与参考不一致，检查坐标约定\n"); return 1; }

    TB_Hop hop;
    tb_init(&hop);

    double p[NP], p0[NP];
    pack(tb_get_params(), p);
    pack(tb_get_params(), p0);

    double *vb0 = malloc(sizeof *vb0 * nk), *cb0 = malloc(sizeof *cb0 * nk);
    double *vb  = malloc(sizeof *vb  * nk), *cb  = malloc(sizeof *cb  * nk);
    eval_bands(p0, &hop, vb0, cb0);

    // 初值先做最优刚性平移（就是 Step A 里那个 Δ）：给 e1/e2 同时加常数等价于整条
    // 能带上移，所以这个偏移有解析解，不必让优化器慢慢爬。
    double shift = 0;
    for (int i = 0; i < nk; ++i)
        shift += (gw_vb[i] - vb0[i]) + (gw_cb[i] - cb0[i]);
    shift /= 2.0 * nk;
    p0[0] += shift;
    p0[1] += shift;
    pack((const TB_Params *)tb_get_params(), p);   // 占位，下面重算
    eval_bands(p0, &hop, vb0, cb0);
    for (int j = 0; j < NP; ++j) p[j] = p0[j];
    printf("\n初值先整体平移 %.3f eV（e1/e2 同加）\n", shift);

    double best = loss(p, &hop);
    printf("\n初值：loss = %.5f eV², RMS = %.4f eV（价带 %.4f / 导带 %.4f）\n",
           best, sqrt(best / (2.0*nk)), rms(vb0, gw_vb), rms(cb0, gw_cb));

    // ---- 逐参数扫描 + 步长减半 ----
    double step[NP];
    for (int j = 0; j < NP; ++j) step[j] = 0.05;   // 初值 50 meV

    for (int it = 1; it <= 2000; ++it) {
        int improved = 0;
        double before = best;
        for (int j = 0; j < NP; ++j) {
            for (int sgn = +1; sgn >= -1; sgn -= 2) {
                double old = p[j];
                p[j] = old + sgn * step[j];
                double l = loss(p, &hop);
                if (l < best) { best = l; improved = 1; break; }
                p[j] = old;
            }
        }
        double maxstep = 0;
        for (int j = 0; j < NP; ++j) if (step[j] > maxstep) maxstep = step[j];
        // 一整轮都动不了、或改进已经很小（< 0.1%），就把步长减半
        if (!improved || best > before * (1.0 - 1e-3)) {
            for (int j = 0; j < NP; ++j) step[j] *= 0.5;
            maxstep *= 0.5;
            if (maxstep < 1e-6) {
                printf("\n最大步长降到 %.0e，停止（第 %d 轮）\n", maxstep, it);
                break;
            }
        }
        if (it == 1 || it % 20 == 0)
            printf("  第 %3d 轮: loss = %10.5f eV², RMS = %.4f eV, 最大步长 = %.1e\n",
                   it, best, sqrt(best / (2.0*nk)), maxstep);
    }

    // ---- 结果 ----
    eval_bands(p, &hop, vb, cb);
    printf("\n拟合后：loss = %.5f eV², RMS = %.4f eV（价带 %.4f / 导带 %.4f）\n",
           best, sqrt(best / (2.0*nk)), rms(vb, gw_vb), rms(cb, gw_cb));

    printf("\n参数（初值 -> 拟合值，单位 eV）：\n");
    const char *nm[NP] = {"e1","e2","t11","t12","t13","t22","t23","t33",
                          "r11","r12","r21","r22","r23","u11","u12","u13","u22","u23","u33"};
    for (int j = 0; j < NP; ++j)
        printf("  %-4s %+8.4f -> %+8.4f\n", nm[j], p0[j], p[j]);

    printf("\n高对称点残差 (eV)：E_GW − E_TB\n");
    printf("%-4s %9s %9s   | %9s %9s\n", "点", "价带(初)", "导带(初)", "价带(拟合)", "导带(拟合)");
    for (int s = 0; s < 3; ++s) {
        int i = 0;
        for (int j = 0; j < nk; ++j) if (fabs(gw.k[j] - SYM_K[s]) < fabs(gw.k[i] - SYM_K[s])) i = j;
        printf("%-4s %9.3f %9.3f   | %9.3f %9.3f\n", SYM_NM[s],
               gw_vb[i] - vb0[i], gw_cb[i] - cb0[i], gw_vb[i] - vb[i], gw_cb[i] - cb[i]);
    }

    FILE *out = fopen("data/gw_fit_bands.dat", "w");
    if (!out) { fprintf(stderr, "[ERROR] cannot write data/gw_fit_bands.dat\n"); return 1; }
    fprintf(out, "# col1: k (Å^-1)\n");
    fprintf(out, "# col2-3: GW 参考价带/导带\n");
    fprintf(out, "# col4-5: TB 初值（params.h 的文献参数）价带/导带\n");
    fprintf(out, "# col6-7: TB 拟合后 价带/导带\n");
    for (int i = 0; i < nk; ++i)
        fprintf(out, "%10.6f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                gw.k[i], gw_vb[i], gw_cb[i], vb0[i], cb0[i], vb[i], cb[i]);
    fclose(out);
    printf("\n已写出 data/gw_fit_bands.dat\n");
    return 0;
}
