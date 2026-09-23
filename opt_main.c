// 全路径最小二乘拟合（最朴素的版本：逐参数扫描 + 步长减半）
//
//   min_p  Σ_{n,k} [E_TB(p; n, k) − E_GW(n, k)]²
//   n = 价带、第一导带、第二导带（TB 第 1/2/3 条 ↔ GW 第 9/10/11 条，见 NPAIR）
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

// 这个拟合要 19 个参数；不打开 NNN_MODEL 的话 r/u 恒为 0，等于退化成只有 NN 的模型，
// 而且不会报任何错 —— 所以这里直接编译期拦住。
#if !NNN_MODEL
#error "opt_main.c 必须用 -DNNN_MODEL=1 编译（见 build.sh）"
#endif

#define NP 19            // e1,e2 + 6 个 t + 5 个 r + 6 个 u

// 参与拟合的带对：TB 第 TB_PAIR[i] 条 <-> GW 第 GW_PAIR[i] 条
// 2 条 = 只拟合价带+第一导带；3 条 = 再加上第二条导带
#define NPAIR 3
static const int TB_PAIR[3] = {1, 2, 3};
static const int GW_PAIR[3] = {9, 10, 11};

// 重连用的候选能级范围（1-based）：覆盖最低的几条导带
#define CAND_FIRST 10
#define CAND_LAST  14
#define NCAND      (CAND_LAST - CAND_FIRST + 1)

// 参考能带的配对方式：
//   0 = 直接用排序块（"第 10/11 条带"）——交叉处排序曲线互换身份、有折点
//   1 = 先用斜率外推重连成物理能带再拟合（默认）
// wannier90 的能带输出是按能量排序的，两条物理能带交叉处会互换身份；
// 而 TB 模型自己的本征值（herm3_eigvals）同样是排序的 —— 两边都要重连才一致。
#ifndef RECONNECT
#define RECONNECT 1
#endif

// GW 路径上的高对称点位置 (Å^-1)，用于报告分区残差
static const double SYM_K[3]  = {0.0, 1.31633, 1.97450};
static const char  *SYM_NM[3] = {"Γ", "K", "M"};

static int nk;
static double *kx, *ky;        // TB 单位（1/a）的笛卡尔 k
static double *gw_ref[3];      // 各带对的目标能量（重连后的物理能带）
static double *gsorted[3];     // 排序块 9/10/11（仅用于诊断）
static double *recon;          // GW 重连结果 [NCAND*nk]
static double *mb_raw;         // 模型本征值（排序）[3*nk]
static double *mb_trk;         // 模型本征值（重连后）[3*nk]
static int     track_anchor;   // 重连锚点（K 点）

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
// 算模型在全部 k 点的三条带。RECONNECT=1 时先重连：拟合参数下模型自己的两条
// 导带在 M-Γ 之间也会靠近/交叉，而 herm3_eigvals 返回的是排序本征值。
static void model_bands(const double *p, TB_Hop *hop)
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
        for (int b = 0; b < 3; ++b) mb_raw[b*nk + i] = e[b];
    }
#if RECONNECT
    band_track(nk, 3, mb_raw, track_anchor, mb_trk);
#else
    for (int b = 0; b < 3; ++b)
        for (int i = 0; i < nk; ++i) mb_trk[b*nk + i] = mb_raw[b*nk + i];
#endif
}

static double loss(const double *p, TB_Hop *hop)
{
    TB_Params q;
    unpack(p, &q);
    tb_set_params(&q);
    tb_update_hopping(hop);

    model_bands(p, hop);
    double s = 0;
    for (int i = 0; i < nk; ++i)
        for (int b = 0; b < NPAIR; ++b) {
            double d = mb_trk[(TB_PAIR[b] - 1) * nk + i] - gw_ref[b][i];
            s += d * d;
        }
    return s;
}

// 把当前参数下的能带算出来（用于输出/画图）
static void eval_bands(const double *p, TB_Hop *hop, double *out)
{
    model_bands(p, hop);
    for (int b = 0; b < NPAIR; ++b)
        for (int i = 0; i < nk; ++i) out[b*nk + i] = mb_trk[(TB_PAIR[b] - 1) * nk + i];
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
    if (gw.nband < GW_PAIR[NPAIR-1]) { fprintf(stderr, "[ERROR] GW 数据带数不够\n"); return 1; }

    if (kpt_read("MoS2-GWBSE-data/wannier90_band.kpt", &nk, &kx, &ky)) return 1;
    if (nk != gw.nk) {
        fprintf(stderr, "[ERROR] k 点数不一致：kpt %d，band %d\n", nk, gw.nk);
        return 1;
    }
    for (int b = 0; b < 3; ++b) gsorted[b] = gw.E + (GW_PAIR[0] - 1 + b) * nk;

    // ---- 重连：把"每个 k 点排序好的能级"接回物理能带 ----
    {
        double *cand = malloc(sizeof(double) * NCAND * nk);
        for (int c = 0; c < NCAND; ++c)
            memcpy(cand + c * nk, gw.E + (CAND_FIRST - 1 + c) * nk, sizeof(double) * nk);
        recon = malloc(sizeof(double) * NCAND * nk);

        int anchor = 0;
        for (int i = 1; i < nk; ++i)
            if (fabs(gw.k[i] - 1.31633) < fabs(gw.k[anchor] - 1.31633)) anchor = i;
        track_anchor = anchor;
        mb_raw = malloc(sizeof(double) * 3 * nk);
        mb_trk = malloc(sizeof(double) * 3 * nk);
        if (band_track(nk, NCAND, cand, anchor, recon)) {
            fprintf(stderr, "[ERROR] 参考能带重连失败\n");
            return 1;
        }

        double r0 = 0, r1 = 0;
        for (int c = 0; c < NCAND; ++c)
            for (int i = 1; i < nk - 1; ++i) {
                double d0 = cand[c*nk+i+1] - 2*cand[c*nk+i] + cand[c*nk+i-1];
                double d1 = recon[c*nk+i+1] - 2*recon[c*nk+i] + recon[c*nk+i-1];
                r0 += d0*d0;
                r1 += d1*d1;
            }
        printf("参考能带重连（斜率外推，锚点 = K）：\n");
        printf("  总粗糙度 Σ(二阶差分)²：排序块 %.5f -> 重连后 %.5f\n", r0, r1);
        for (int c = 0; c < 3; ++c) {
            int diff = 0;
            for (int i = 0; i < nk; ++i)
                if (fabs(recon[c*nk+i] - cand[c*nk+i]) > 1e-9) diff++;
            printf("    物理带%d：与排序块%d 不同的 k 点 %d/%d（%.0f%%）\n",
                   c + 1, CAND_FIRST + c, diff, nk, 100.0 * diff / nk);
        }
        printf("\n");

        // 写出重连结果，供 plot_reconnect.gp 直观检查
        FILE *fo = fopen("data/band_track.dat", "w");
        if (fo) {
            fprintf(fo, "# col1: k (Å^-1)\n");
            fprintf(fo, "# col2-%d : 排序块 %d..%d（原始数据）\n", 1 + NCAND, CAND_FIRST, CAND_LAST);
            fprintf(fo, "# col%d-%d : 重连后的物理带（编号按 K 点能量顺序）\n", 2 + NCAND, 1 + 2*NCAND);
            for (int i = 0; i < nk; ++i) {
                fprintf(fo, "%10.6f", gw.k[i]);
                for (int c = 0; c < NCAND; ++c) fprintf(fo, " %10.5f", cand[c*nk + i]);
                for (int b = 0; b < NCAND; ++b) fprintf(fo, " %10.5f", recon[b*nk + i]);
                // 再附上：该曲线在这个 k 用的是哪个排序块（用于画"身份图"）
                for (int b = 0; b < NCAND; ++b) {
                    int src = 0;
                    for (int c = 1; c < NCAND; ++c)
                        if (fabs(cand[c*nk+i] - recon[b*nk+i]) < fabs(cand[src*nk+i] - recon[b*nk+i])) src = c;
                    fprintf(fo, " %4d", CAND_FIRST + src);
                }
                fprintf(fo, "\n");
            }
            fclose(fo);
        }
    }

#if RECONNECT
    // 价带用排序块 9（本身很干净），两条导带用重连后的物理能带
    gw_ref[0] = gsorted[0];
    gw_ref[1] = recon + 0 * nk;
    gw_ref[2] = recon + 1 * nk;
#else
    for (int b = 0; b < NPAIR; ++b) gw_ref[b] = gsorted[b];
#endif

    // 自检：把 .kpt 重建的 k 点累加成路径长度，应与参考文件的 k 轴一致
    double cum = 0, dk = 0;
    for (int i = 0; i < nk; ++i) {
        if (i > 0) cum += hypot(kx[i] - kx[i-1], ky[i] - ky[i-1]);
        double d = fabs(cum / A_ANG - gw.k[i]);
        if (d > dk) dk = d;
    }
    printf("拟合：GW 参考 vs 三能带模型（NN + NNN + TNN，%d 个参数）\n", NP);
    printf("      全路径 %d 个 k 点，%d 条带\n", nk, NPAIR);
    printf("  k 点自检：由 .kpt 重建的路径长度与参考文件最大差 %.2e Å^-1\n", dk);
    if (dk > 1e-4) { fprintf(stderr, "[ERROR] k 点重建与参考不一致，检查坐标约定\n"); return 1; }

    TB_Hop hop;
    tb_init(&hop);

    double p[NP], p0[NP];
    pack(tb_get_params(), p);
    pack(tb_get_params(), p0);

    double *tb0 = malloc(sizeof *tb0 * nk * NPAIR);
    double *tb  = malloc(sizeof *tb  * nk * NPAIR);
    eval_bands(p0, &hop, tb0);

    // 初值先做最优刚性平移（就是 Step A 里那个 Δ）：给 e1/e2 同时加常数等价于整条
    // 能带上移，所以这个偏移有解析解，不必让优化器慢慢爬。
    double shift = 0;
    for (int b = 0; b < NPAIR; ++b)
        for (int i = 0; i < nk; ++i)
            shift += gw_ref[b][i] - tb0[b * nk + i];
    shift /= (double)NPAIR * nk;
    p0[0] += shift;
    p0[1] += shift;
    pack((const TB_Params *)tb_get_params(), p);   // 占位，下面重算
    eval_bands(p0, &hop, tb0);
    for (int j = 0; j < NP; ++j) p[j] = p0[j];
    printf("\n初值先整体平移 %.3f eV（e1/e2 同加）\n", shift);

    double best = loss(p, &hop);
    printf("\n初值：loss = %.5f eV², RMS = %.4f eV\n",
           best, sqrt(best / ((double)NPAIR * nk)));
    for (int b = 0; b < NPAIR; ++b)
        printf("        带%d(TB%d↔GW%d) RMS = %.4f eV\n", b+1, TB_PAIR[b], GW_PAIR[b], rms(tb0 + b*nk, gw_ref[b]));

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
                   it, best, sqrt(best / ((double)NPAIR * nk)), maxstep);
    }

    // ---- 结果 ----
    eval_bands(p, &hop, tb);
    printf("\n拟合后：loss = %.5f eV², RMS = %.4f eV\n",
           best, sqrt(best / ((double)NPAIR * nk)));
    for (int b = 0; b < NPAIR; ++b)
        printf("        带%d(TB%d↔GW%d) RMS = %.4f eV\n", b+1, TB_PAIR[b], GW_PAIR[b], rms(tb + b*nk, gw_ref[b]));

    printf("\n参数（初值 -> 拟合值，单位 eV）：\n");
    const char *nm[NP] = {"e1","e2","t11","t12","t13","t22","t23","t33",
                          "r11","r12","r21","r22","r23","u11","u12","u13","u22","u23","u33"};
    for (int j = 0; j < NP; ++j)
        printf("  %-4s %+8.4f -> %+8.4f\n", nm[j], p0[j], p[j]);

    // 歧义诊断：GW 数据是按能量排序的，两条能带交叉处排序曲线会互换身份。
    // 这里统计"候选能级里有两条靠得太近"的 k 点（那里的配对不可信），
    // 并给出只在这些点之外统计的 RMS。
    {
        const double eps = 0.02;   // 近简并判据 (eV)
        int namb = 0;
        double e1 = 0, e2 = 0, es = 0;
        for (int i = 0; i < nk; ++i) {
            int amb = 0;
            for (int a = 0; a < 3 && !amb; ++a)
                for (int b = a + 1; b < 3; ++b) {
                    double da = gsorted[a][i] - gsorted[b][i];
                    if (fabs(da) < eps) amb = 1;
                }
            if (amb) namb++;
            for (int b = 0; b < NPAIR; ++b) {
                double d = gw_ref[b][i] - tb[b * nk + i];
                if (amb) e2 += d * d; else e1 += d * d;
                double ds = gsorted[b][i] - tb[b * nk + i];
                es += ds * ds;
            }
        }
        printf("\n配对歧义诊断（|ΔE| < %.2f eV 视为分不清是哪条物理能带）：\n", eps);
        printf("  歧义 k 点 %d/%d（%.0f%%）\n", namb, nk, 100.0 * namb / nk);
        printf("  无歧义点上的 RMS = %.4f eV（全部点 %.4f eV）\n",
               sqrt(e1 / ((double)NPAIR * (nk - namb))), sqrt((e1 + e2) / ((double)NPAIR * nk)));
        printf("  对照：若直接用排序块当目标，RMS = %.4f eV\n", sqrt(es / ((double)NPAIR * nk)));
    }

    printf("\n高对称点残差 E_GW − E_TB (eV)：初值 -> 拟合后\n");
    printf("%-4s", "点");
    for (int b = 0; b < NPAIR; ++b) printf("   TB%d↔GW%-2d", TB_PAIR[b], GW_PAIR[b]);
    printf("\n");
    for (int s = 0; s < 3; ++s) {
        int i = 0;
        for (int j = 0; j < nk; ++j) if (fabs(gw.k[j] - SYM_K[s]) < fabs(gw.k[i] - SYM_K[s])) i = j;
        printf("%-4s", SYM_NM[s]);
        for (int b = 0; b < NPAIR; ++b)
            printf("  %+.3f->%+.3f", gw_ref[b][i] - tb0[b*nk+i], gw_ref[b][i] - tb[b*nk+i]);
        printf("\n");
    }

    FILE *out = fopen("data/gw_fit_bands.dat", "w");
    if (!out) { fprintf(stderr, "[ERROR] cannot write data/gw_fit_bands.dat\n"); return 1; }
    fprintf(out, "# col1: k (Å^-1)\n");
    fprintf(out, "# col2-%d : GW 参考（带 %d/%d/%d）\n", 1+NPAIR, GW_PAIR[0], GW_PAIR[1 % NPAIR], GW_PAIR[2 % NPAIR]);
    fprintf(out, "# col%d-%d : TB 初值（params.h 文献参数）\n", 2+NPAIR, 1+2*NPAIR);
    fprintf(out, "# col%d-%d : TB 拟合后\n", 2+2*NPAIR, 1+3*NPAIR);
    for (int i = 0; i < nk; ++i) {
        fprintf(out, "%10.6f", gw.k[i]);
        for (int b = 0; b < NPAIR; ++b) fprintf(out, " %10.5f", gw_ref[b][i]);
        for (int b = 0; b < NPAIR; ++b) fprintf(out, " %10.5f", tb0[b*nk+i]);
        for (int b = 0; b < NPAIR; ++b) fprintf(out, " %10.5f", tb[b*nk+i]);
        fprintf(out, "\n");
    }
    fclose(out);
    printf("\n已写出 data/gw_fit_bands.dat\n");
    return 0;
}
