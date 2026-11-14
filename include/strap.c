#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "strap.h"

#define STRAP_LINE 1024

static int line_is_blank(const char *s)
{
    for (; *s; ++s) if (*s != ' ' && *s != '\t' && *s != '\n' && *s != '\r') return 0;
    return 1;
}

// 数一行有几个数，用来区分块式（k E）和列式（k E1 E2 ...）
static int count_fields(const char *s)
{
    int n = 0;
    while (*s) {
        while (*s == ' ' || *s == '\t') ++s;
        if (!*s || *s == '\n' || *s == '\r') break;
        ++n;
        while (*s && *s != ' ' && *s != '\t' && *s != '\n' && *s != '\r') ++s;
    }
    return n;
}

void strap_init(StrapBands *s)
{
    if (!s) return;
    s->nk = 0;
    s->nband = 0;
    s->k = NULL;
    s->E = NULL;
}

int strap_alloc(StrapBands *s, int nk, int nband)
{
    if (!s || nk <= 0 || nband <= 0) return 1;
    strap_free(s);
    s->k = malloc(sizeof *s->k * (size_t)nk);
    s->E = malloc(sizeof *s->E * (size_t)nk * (size_t)nband);
    if (!s->k || !s->E) {
        fprintf(stderr, "[ERROR] strap_alloc: 内存不足 (nk=%d nband=%d)\n", nk, nband);
        strap_free(s);
        return 1;
    }
    s->nk = nk;
    s->nband = nband;
    return 0;
}

void strap_free(StrapBands *s)
{
    if (!s) return;
    free(s->k);
    free(s->E);
    s->k = NULL;
    s->E = NULL;
    s->nk = 0;
    s->nband = 0;
}

// 块式：每块 "k E"，块间空行，各块点数相同且 k 轴逐点一致
static int read_block(const char *path, StrapBands *s)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", path);
        return 1;
    }

    char line[STRAP_LINE];
    int nband = 0, nk = -1, cnt = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) {
            if (cnt > 0) {
                if (nk < 0) nk = cnt;
                else if (cnt != nk) {
                    fprintf(stderr, "[ERROR] %s: 各块点数不一致（%d vs %d）\n", path, cnt, nk);
                    fclose(f);
                    return 1;
                }
                ++nband;
                cnt = 0;
            }
            continue;
        }
        if (line[0] == '#') continue;
        ++cnt;
    }
    if (cnt > 0) {
        if (nk < 0) nk = cnt;
        else if (cnt != nk) {
            fprintf(stderr, "[ERROR] %s: 最后一块点数 %d 与前面 %d 不一致\n", path, cnt, nk);
            fclose(f);
            return 1;
        }
        ++nband;
    }
    fclose(f);

    if (nband <= 0 || nk <= 0) {
        fprintf(stderr, "[ERROR] %s: 没读出数据（nband=%d nk=%d）\n", path, nband, nk);
        return 1;
    }
    if (strap_alloc(s, nk, nband)) return 1;

    f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", path);
        strap_free(s);
        return 1;
    }
    int band = 0, i = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) {
            if (i > 0) {
                ++band;
                i = 0;
            }
            continue;
        }
        if (line[0] == '#') continue;
        double kv = 0, ev = 0;
        if (sscanf(line, "%lf %lf", &kv, &ev) != 2) continue;
        if (band >= nband || i >= nk) break;
        if (band == 0) s->k[i] = kv;
        else if (fabs(kv - s->k[i]) > 1e-8) {
            fprintf(stderr, "[ERROR] %s: 第 %d 块的 k 轴与第 1 块不一致（k[%d]: %.10g vs %.10g）\n",
                    path, band + 1, i, kv, s->k[i]);
            fclose(f);
            strap_free(s);
            return 1;
        }
        s->E[band * nk + i] = ev;
        ++i;
    }
    fclose(f);
    return 0;
}

// 列式：每行 "k E1 E2 ... En"
static int read_column(const char *path, StrapBands *s)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", path);
        return 1;
    }

    char line[STRAP_LINE];
    int ncol = -1, nk = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        int nc = count_fields(line);
        if (ncol < 0) ncol = nc;
        else if (nc != ncol) {
            fprintf(stderr, "[ERROR] %s: 第 %d 行列数 %d 与前面 %d 不一致\n", path, nk + 1, nc, ncol);
            fclose(f);
            return 1;
        }
        ++nk;
    }
    fclose(f);

    if (nk <= 0 || ncol < 2) {
        fprintf(stderr, "[ERROR] %s: 没读出数据（nk=%d ncol=%d）\n", path, nk, ncol);
        return 1;
    }
    if (strap_alloc(s, nk, ncol - 1)) return 1;

    f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", path);
        strap_free(s);
        return 1;
    }
    double *row = malloc(sizeof *row * (size_t)ncol);
    if (!row) {
        fprintf(stderr, "[ERROR] 内存不足\n");
        fclose(f);
        strap_free(s);
        return 1;
    }
    int i = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        if (i >= nk) break;
        // 逐个数手动切，避免 64 列时反复改格式串
        const char *p = line;
        int nc = 0;
        while (nc < ncol) {
            char *end = NULL;
            double v = strtod(p, &end);
            if (end == p) break;
            row[nc++] = v;
            p = end;
        }
        if (nc != ncol) {
            fprintf(stderr, "[ERROR] %s: 第 %d 行只读到 %d 个数\n", path, i + 1, nc);
            free(row);
            fclose(f);
            strap_free(s);
            return 1;
        }
        s->k[i] = row[0];
        for (int b = 0; b < s->nband; ++b) s->E[b * nk + i] = row[b + 1];
        ++i;
    }
    free(row);
    fclose(f);
    return 0;
}

int strap_read(const char *path, StrapBands *s)
{
    if (!path || !s) return 1;
    memset(s, 0, sizeof *s);

    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", path);
        return 1;
    }
    char line[STRAP_LINE];
    int ncol = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        ncol = count_fields(line);
        break;
    }
    fclose(f);
    if (ncol <= 0) {
        fprintf(stderr, "[ERROR] %s: 文件里没有数据行\n", path);
        return 1;
    }
    return (ncol > 2) ? read_column(path, s) : read_block(path, s);
}

int strap_write(const char *path, const StrapBands *s, int layout, int prec, const char *header)
{
    if (!path || !s || !s->E || s->nk <= 0 || s->nband <= 0) return 1;
    if (prec <= 0) prec = 17;

    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "[ERROR] 写不了 %s\n", path);
        return 1;
    }
    if (header) fprintf(f, "# %s\n", header);
    fprintf(f, "# nk = %d, nband = %d\n", s->nk, s->nband);

    if (layout == STRAP_LAYOUT_BLOCK) {
        for (int b = 0; b < s->nband; ++b) {
            if (b) fprintf(f, "\n");
            for (int i = 0; i < s->nk; ++i)
                fprintf(f, "%.10g %.*g\n", s->k[i], prec, s->E[b * s->nk + i]);
        }
    } else {
        for (int i = 0; i < s->nk; ++i) {
            fprintf(f, "%.10g", s->k[i]);
            for (int b = 0; b < s->nband; ++b)
                fprintf(f, " %.*g", prec, s->E[b * s->nk + i]);
            fprintf(f, "\n");
        }
    }
    fclose(f);
    return 0;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

int strap_sort_by_energy(StrapBands *s)
{
    if (!s || !s->E || s->nk <= 0 || s->nband <= 0) return 1;
    double *col = malloc(sizeof *col * (size_t)s->nband);
    if (!col) {
        fprintf(stderr, "[ERROR] 内存不足\n");
        return 1;
    }
    for (int i = 0; i < s->nk; ++i) {
        for (int b = 0; b < s->nband; ++b) col[b] = s->E[b * s->nk + i];
        qsort(col, (size_t)s->nband, sizeof *col, cmp_double);
        for (int b = 0; b < s->nband; ++b) s->E[b * s->nk + i] = col[b];
    }
    free(col);
    return 0;
}

int strap_find_anchor(const StrapBands *s, int *anchor)
{
    if (!s || !s->E || !anchor || s->nk < 3 || s->nband < 2) return 1;
    double *col = malloc(sizeof *col * (size_t)s->nband);
    if (!col) {
        fprintf(stderr, "[ERROR] 内存不足\n");
        return 1;
    }
    int best_i = s->nk / 2;
    double best_gap = -1;
    for (int i = 0; i < s->nk; ++i) {
        for (int b = 0; b < s->nband; ++b) col[b] = s->E[b * s->nk + i];
        qsort(col, (size_t)s->nband, sizeof *col, cmp_double);
        double gap = col[1] - col[0];
        for (int b = 2; b < s->nband; ++b) {
            double d = col[b] - col[b - 1];
            if (d < gap) gap = d;
        }
        if (gap > best_gap) {
            best_gap = gap;
            best_i = i;
        }
    }
    free(col);
    // 锚点两侧都要有邻居，端点处没法外推
    if (best_i == 0) best_i = 1;
    if (best_i >= s->nk - 1) best_i = s->nk - 2;
    *anchor = best_i;
    return 0;
}

// 匈牙利算法（O(n³)）：求 n×n 代价矩阵的最优指派，assign[i] = 第 i 行分到的列
static int hungarian(int n, const double *cost, int *assign)
{
    double *u = calloc((size_t)n + 1, sizeof *u);
    double *v = calloc((size_t)n + 1, sizeof *v);
    int *p = calloc((size_t)n + 1, sizeof *p);
    int *way = calloc((size_t)n + 1, sizeof *way);
    double *minv = malloc(sizeof *minv * ((size_t)n + 1));
    char *used = malloc((size_t)n + 1);
    if (!u || !v || !p || !way || !minv || !used) {
        fprintf(stderr, "[ERROR] hungarian: 内存不足\n");
        free(u); free(v); free(p); free(way); free(minv); free(used);
        return 1;
    }

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        for (int j = 0; j <= n; ++j) {
            minv[j] = 1e300;
            used[j] = 0;
        }
        do {
            used[j0] = 1;
            int i0 = p[j0], j1 = 0;
            double delta = 1e300;
            for (int j = 1; j <= n; ++j) {
                if (used[j]) continue;
                double cur = cost[(size_t)(i0 - 1) * n + (j - 1)] - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    for (int j = 1; j <= n; ++j) assign[p[j] - 1] = j - 1;
    free(u); free(v); free(p); free(way); free(minv); free(used);
    return 0;
}

int strap_reconnect(StrapBands *s, int anchor)
{
    if (!s || !s->E || s->nk < 3 || s->nband < 2) return 1;
    if (anchor <= 0 || anchor >= s->nk - 1) {
        fprintf(stderr, "[ERROR] strap_reconnect: 锚点 %d 越界（需要 0 < anchor < %d）\n",
                anchor, s->nk - 1);
        return 1;
    }
    int nk = s->nk, nb = s->nband;

    int *asg = malloc(sizeof *asg * (size_t)nk * (size_t)nb);
    double *cost = malloc(sizeof *cost * (size_t)nb * (size_t)nb);
    int *perm = malloc(sizeof *perm * (size_t)nb);
    double *cur = malloc(sizeof *cur * (size_t)nb);
    double *prev = malloc(sizeof *prev * (size_t)nb);
    double *pred = malloc(sizeof *pred * (size_t)nb);
    double *now = malloc(sizeof *now * (size_t)nb);
    double *out = malloc(sizeof *out * (size_t)nk * (size_t)nb);
    if (!asg || !cost || !perm || !cur || !prev || !pred || !now || !out) {
        fprintf(stderr, "[ERROR] strap_reconnect: 内存不足\n");
        free(asg); free(cost); free(perm); free(cur); free(prev); free(pred); free(now); free(out);
        return 1;
    }

    // 初始身份 = 各 k 点内部的排列顺序
    for (int i = 0; i < nk; ++i)
        for (int b = 0; b < nb; ++b) asg[(size_t)i * nb + b] = b;

    for (int dir = 1; dir >= -1; dir -= 2) {
        for (int b = 0; b < nb; ++b) {
            cur[b] = s->E[b * nk + anchor];
            prev[b] = cur[b];
        }
        for (int i = anchor + dir; i >= 0 && i < nk; i += dir) {
            for (int b = 0; b < nb; ++b) {
                pred[b] = 2.0 * cur[b] - prev[b];
                now[b] = s->E[b * nk + i];
            }
            // 代价 = 该带接到该能级后，偏离线性外推预测多少
            for (int b = 0; b < nb; ++b)
                for (int c = 0; c < nb; ++c) {
                    double d = now[c] - pred[b];
                    cost[(size_t)b * nb + c] = d * d;
                }
            if (hungarian(nb, cost, perm)) {
                free(asg); free(cost); free(perm); free(cur);
                free(prev); free(pred); free(now); free(out);
                return 1;
            }
            for (int b = 0; b < nb; ++b) {
                asg[(size_t)i * nb + b] = perm[b];
                prev[b] = cur[b];
                cur[b] = now[perm[b]];
            }
        }
    }

    for (int i = 0; i < nk; ++i)
        for (int b = 0; b < nb; ++b)
            out[b * nk + i] = s->E[asg[(size_t)i * nb + b] * nk + i];
    memcpy(s->E, out, sizeof *s->E * (size_t)nk * (size_t)nb);

    free(asg); free(cost); free(perm); free(cur);
    free(prev); free(pred); free(now); free(out);
    return 0;
}

double strap_roughness(const StrapBands *s)
{
    if (!s || !s->E || s->nk < 3) return 0.0;
    double sum = 0;
    for (int b = 0; b < s->nband; ++b)
        for (int i = 1; i + 1 < s->nk; ++i) {
            double d2 = s->E[b * s->nk + i + 1] - 2.0 * s->E[b * s->nk + i]
                        + s->E[b * s->nk + i - 1];
            sum += d2 * d2;
        }
    return sum;
}
