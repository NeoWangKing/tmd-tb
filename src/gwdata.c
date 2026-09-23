#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "gwdata.h"

static int line_is_blank(const char *s)
{
    for (; *s; ++s) if (*s != ' ' && *s != '\t' && *s != '\n' && *s != '\r') return 0;
    return 1;
}

int gw_read(const char *path, GwData *g)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[ERROR] cannot open %s\n", path); return 1; }

    // 第一遍：数每条带的点数，并检查各带点数一致
    char line[512];
    int nband = 0, nk = -1, cnt = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) {
            if (cnt > 0) {
                if (nk < 0) nk = cnt; else if (cnt != nk) return 1;
                nband++; cnt = 0;
            }
            continue;
        }
        if (line[0] == '#') continue;
        cnt++;
    }
    if (cnt > 0) { if (nk < 0) nk = cnt; else if (cnt != nk) return 1; nband++; }
    fclose(f);
    if (nband <= 0 || nk <= 0 || nband > GW_NBAND_MAX || nk > NK_MAX) {
        fprintf(stderr, "[ERROR] %s: nband=%d nk=%d 不合理\n", path, nband, nk);
        return 1;
    }

    g->nband = nband; g->nk = nk;
    g->k = malloc(sizeof *g->k * nk);
    g->E = malloc(sizeof *g->E * nk * nband);
    if (!g->k || !g->E) { fprintf(stderr, "[ERROR] out of memory\n"); return 1; }

    // 第二遍：读入；k 轴以第一块为准，其余块必须一致
    f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[ERROR] cannot open %s\n", path); return 1; }
    int band = 0, i = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) { if (i > 0) { band++; i = 0; } continue; }
        if (line[0] == '#') continue;
        double kv, ev;
        if (sscanf(line, "%lf %lf", &kv, &ev) != 2) continue;
        if (band >= nband || i >= nk) break;
        if (band == 0) g->k[i] = kv;
        else if (fabs(kv - g->k[i]) > 1e-8) {
            fprintf(stderr, "[ERROR] %s: band %d 的 k 轴与第 1 条带不一致\n", path, band + 1);
            fclose(f); return 1;
        }
        g->E[band * nk + i] = ev;
        i++;
    }
    fclose(f);
    return 0;
}

// 读 TB 能带表：跳过 # 注释行，每行 1 + nband 列（k 与 nband 个能量）
int tb_read(const char *path, int nband, int *nk_out, double **k_out, double **e_out)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[ERROR] cannot open %s\n", path); return 1; }

    int cap = 1024, n = 0;
    double *k = malloc(sizeof *k * cap);
    double *e = malloc(sizeof *e * cap * nband);
    if (!k || !e) { fprintf(stderr, "[ERROR] out of memory\n"); fclose(f); return 1; }

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
            if (!nk2 || !ne2) { fprintf(stderr, "[ERROR] out of memory\n"); fclose(f); return 1; }
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
int dedup(int n, const double *k, const double *e, int nband, double *k2, double *e2)
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

// .kpt 里的分数坐标用的是六角惯用基 b1 = (2π/a)(1, 1/√3)、b2 = (2π/a)(1, -1/√3)，
// 与 tb.c 用的基不同但张成同一个倒格子。换算到 TB 单位（1/a）：
//     kx = 2π(f1 + f2),  ky = 2π(f1 - f2)/√3
// （自检：K 点 (1/3,1/3) -> (4π/3, 0)，与 tb.c 里 K=(2/3,1/3) 的结果一致）
int kpt_read(const char *path, int *nk_out, double **kx_out, double **ky_out)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "[ERROR] cannot open %s\n", path); return 1; }

    int n = 0;
    if (fscanf(f, "%d", &n) != 1 || n <= 0 || n > NK_MAX) {
        fprintf(stderr, "[ERROR] %s: 点数异常\n", path);
        fclose(f);
        return 1;
    }
    double *kx = malloc(sizeof *kx * n);
    double *ky = malloc(sizeof *ky * n);
    if (!kx || !ky) { fprintf(stderr, "[ERROR] out of memory\n"); fclose(f); return 1; }

    for (int i = 0; i < n; ++i) {
        double f1, f2, f3, w;
        if (fscanf(f, "%lf %lf %lf %lf", &f1, &f2, &f3, &w) != 4) {
            fprintf(stderr, "[ERROR] %s: 第 %d 行格式不对\n", path, i + 2);
            fclose(f);
            return 1;
        }
        kx[i] = 2*PI*(f1 + f2);
        ky[i] = 2*PI*(f1 - f2)/sqrt(3.0);
    }
    fclose(f);
    *nk_out = n; *kx_out = kx; *ky_out = ky;
    return 0;
}

// 在候选里找一个排列，使各曲线到"预测位置"的平方距离和最小（递归穷举，ncand 很小）
static void best_perm(const double *cnow, const double *pred, int ncand, int *perm, int *used,
                      int depth, double cur, double *best, int *bestperm)
{
    if (cur >= *best) return;
    if (depth == ncand) {
        *best = cur;
        for (int b = 0; b < ncand; ++b) bestperm[b] = perm[b];
        return;
    }
    for (int c = 0; c < ncand; ++c) {
        if (used[c]) continue;
        used[c] = 1;
        perm[depth] = c;
        double d = cnow[c] - pred[depth];
        best_perm(cnow, pred, ncand, perm, used, depth + 1, cur + d * d, best, bestperm);
        used[c] = 0;
    }
}

int gw_reconnect(int nk, int ncand, const double *cand, int anchor, double *out)
{
    if (ncand < 2 || ncand > 8 || anchor <= 0 || anchor >= nk - 1) return 1;

    int *asg = malloc(sizeof(int) * nk * ncand);
    if (!asg) return 1;
    for (int i = 0; i < nk; ++i)
        for (int b = 0; b < ncand; ++b) asg[i*ncand + b] = b;   // 初始 = 排序顺序

    // 从 anchor 向两侧逐点分配：用"斜率外推"预测下一时刻的位置，
    // 再在所有排列里选距离平方和最小的那个 —— 等价于保持斜率连续。
    for (int step = +1; step >= -1; step -= 2) {
        double prev[8], curp[8], cnow[8], pred[8];
        for (int b = 0; b < ncand; ++b) {
            curp[b] = cand[b*nk + anchor];
            prev[b] = curp[b];
        }
        for (int i = anchor + step; i >= 0 && i < nk; i += step) {
            for (int b = 0; b < ncand; ++b) {
                pred[b] = 2*curp[b] - prev[b];
                cnow[b] = cand[b*nk + i];
            }
            int perm[8], used[8] = {0}, bestperm[8];
            double best = 1e300, cur = 0;
            best_perm(cnow, pred, ncand, perm, used, 0, cur, &best, bestperm);
            for (int b = 0; b < ncand; ++b) {
                asg[i*ncand + b] = bestperm[b];
                prev[b] = curp[b];
                curp[b] = cnow[bestperm[b]];
            }
        }
    }

    for (int i = 0; i < nk; ++i)
        for (int b = 0; b < ncand; ++b)
            out[b*nk + i] = cand[asg[i*ncand + b] * nk + i];

    free(asg);
    return 0;
}
