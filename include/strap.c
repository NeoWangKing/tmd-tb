#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "strap.h"

#define BUF_SIZE       1024
#define BAND_INIT_CAP  256
#define STRAP_INIT_CAP 8

static int line_is_blank(const char *s)
{
    for (; *s; ++s) if (*s != ' ' && *s != '\t' && *s != '\n' && *s != '\r') return 0;
    return 1;
}

// 数一行有几个数，用来区分块式（k E）和列式（k E1 E2 ...）
static int count_datas(const char *s)
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

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

//=====================================
// 一条带
//=====================================

void band_init(Band *b)
{
    if (!b) return;
    b->keys = NULL;
    b->values = NULL;
    b->count = 0;
    b->capacity = 0;
}

void band_free(Band *b)
{
    if (!b) return;
    free(b->keys);
    free(b->values);
    band_init(b);
}

int band_reserve(Band *b, size_t cap)
{
    if (!b) return 1;
    if (cap <= b->capacity) return 0;

    size_t nc = b->capacity ? b->capacity : BAND_INIT_CAP;
    while (nc < cap) nc *= 2;

    // 先接临时指针：realloc 失败时原指针不能丢
    double *tmp = realloc(b->keys, nc * sizeof *tmp);
    if (!tmp) {
        fprintf(stderr, "[ERROR] band_reserve: keys 扩容到 %zu 失败\n", nc);
        return 1;
    }
    b->keys = tmp;
    tmp = realloc(b->values, nc * sizeof *tmp);
    if (!tmp) {
        fprintf(stderr, "[ERROR] band_reserve: values 扩容到 %zu 失败\n", nc);
        return 1;
    }
    b->values = tmp;
    b->capacity = nc;
    return 0;
}

int band_push(Band *b, double key, double value)
{
    if (!b) return 1;
    if (band_reserve(b, b->count + 1)) return 1;
    b->keys[b->count] = key;
    b->values[b->count] = value;
    ++b->count;
    return 0;
}

//=====================================
// 一组带
//=====================================

void strap_init(StrapBands *s)
{
    if (!s) return;
    s->bands = NULL;
    s->count = 0;
    s->capacity = 0;
}

void strap_free(StrapBands *s)
{
    if (!s) return;
    for (size_t b = 0; b < s->count; ++b) band_free(&s->bands[b]);
    free(s->bands);
    strap_init(s);
}

size_t strap_nb(const StrapBands *s)
{
    return s ? s->count : 0;
}

size_t strap_nk(const StrapBands *s)
{
    return (s && s->count > 0) ? s->bands[0].count : 0;
}

static int strap_reserve(StrapBands *s, size_t cap)
{
    if (cap <= s->capacity) return 0;

    size_t nc = s->capacity ? s->capacity : STRAP_INIT_CAP;
    while (nc < cap) nc *= 2;

    Band *tmp = realloc(s->bands, nc * sizeof *tmp);
    if (!tmp) {
        fprintf(stderr, "[ERROR] strap_reserve: bands 扩容到 %zu 失败\n", nc);
        return 1;
    }
    s->bands = tmp;
    s->capacity = nc;
    return 0;
}

Band *strap_push_band(StrapBands *s)
{
    if (!s) return NULL;
    if (strap_reserve(s, s->count + 1)) return NULL;
    band_init(&s->bands[s->count]);
    return &s->bands[s->count++];
}

int strap_check_uniform(const StrapBands *s)
{
    if (!s || s->count == 0) {
        fprintf(stderr, "[ERROR] check_uniform: 一条带都没有\n");
        return 1;
    }
    const Band *ref = &s->bands[0];
    if (ref->count == 0) {
        fprintf(stderr, "[ERROR] check_uniform: 第 1 条带没有点\n");
        return 1;
    }
    for (size_t b = 1; b < s->count; ++b) {
        const Band *q = &s->bands[b];
        if (q->count != ref->count) {
            fprintf(stderr, "[ERROR] check_uniform: 第 %zu 条带有 %zu 个点，第 1 条有 %zu 个\n",
                    b + 1, q->count, ref->count);
            return 1;
        }
        for (size_t i = 0; i < ref->count; ++i) {
            double tol = 1e-8 * (1.0 + fabs(ref->keys[i]));
            if (fabs(q->keys[i] - ref->keys[i]) > tol) {
                fprintf(stderr,
                        "[ERROR] check_uniform: 第 %zu 条带第 %zu 个横坐标 %.10g 与第 1 条的 %.10g 不一致\n",
                        b + 1, i, q->keys[i], ref->keys[i]);
                return 1;
            }
        }
    }
    return 0;
}

//=====================================
// 读
//=====================================

// 块式：每块 "k E"，块与块之间用空行隔开，一块就是一条带
static int read_block(const char *file_path, StrapBands *s)
{
    FILE *f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", file_path);
        return 1;
    }

    char line[BUF_SIZE];
    Band *cur = NULL;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) {
            cur = NULL;     // 空行 = 一条带的结束
            continue;
        }
        if (line[0] == '#') continue;

        double kv = 0, ev = 0;
        if (sscanf(line, "%lf %lf", &kv, &ev) != 2) continue;
        if (!cur) {
            cur = strap_push_band(s);
            if (!cur) {
                fclose(f);
                return 1;
            }
        }
        if (band_push(cur, kv, ev)) {
            fclose(f);
            return 1;
        }
    }
    fclose(f);

    if (s->count == 0) {
        fprintf(stderr, "[ERROR] %s: 没读出数据\n", file_path);
        return 1;
    }
    return strap_check_uniform(s);
}

// 列式：每行 "k E1 E2 ... En"，第 c 个能量攒成第 c 条带
static int read_column(const char *file_path, StrapBands *s)
{
    FILE *f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", file_path);
        return 1;
    }

    char line[BUF_SIZE];
    int ncol = -1;
    size_t nk = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        int nc = count_datas(line);
        if (ncol < 0) ncol = nc;
        else if (nc != ncol) {
            fprintf(stderr, "[ERROR] %s: 第 %zu 行有 %d 个数，前面是 %d 个\n",
                    file_path, nk + 1, nc, ncol);
            fclose(f);
            return 1;
        }
        ++nk;
    }
    fclose(f);

    if (nk == 0 || ncol < 2) {
        fprintf(stderr, "[ERROR] %s: 没读出数据（nk=%zu ncol=%d）\n", file_path, nk, ncol);
        return 1;
    }

    f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", file_path);
        return 1;
    }
    for (int c = 1; c < ncol; ++c)
        if (!strap_push_band(s)) {
            fclose(f);
            return 1;
        }
    for (size_t b = 0; b < s->count; ++b)
        if (band_reserve(&s->bands[b], nk)) {
            fclose(f);
            return 1;
        }

    size_t i = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        if (i >= nk) break;

        const char *p = line;
        char *end = NULL;
        double kv = strtod(p, &end);
        if (end == p) continue;
        p = end;
        for (size_t b = 0; b < s->count; ++b) {
            double ev = strtod(p, &end);
            if (end == p) {
                fprintf(stderr, "[ERROR] %s: 第 %zu 行只读到 %zu 个能量\n", file_path, i + 1, b);
                fclose(f);
                return 1;
            }
            p = end;
            if (band_push(&s->bands[b], kv, ev)) {
                fclose(f);
                return 1;
            }
        }
        ++i;
    }
    fclose(f);
    return strap_check_uniform(s);
}

int strap_read(const char *file_path, StrapBands *s)
{
    if (!file_path || !s) return 1;
    strap_free(s);

    FILE *f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] 打不开 %s\n", file_path);
        return 1;
    }
    char line[BUF_SIZE];
    int ncol = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        ncol = count_datas(line);
        break;
    }
    fclose(f);

    if (ncol <= 0) {
        fprintf(stderr, "[ERROR] %s: 文件里没有数据行\n", file_path);
        return 1;
    }
    return (ncol > 2) ? read_column(file_path, s) : read_block(file_path, s);
}

//=====================================
// 写
//=====================================

int strap_write(const char *file_path, const StrapBands *s, int layout, int prec,
                const char *header)
{
    if (!file_path || !s || s->count == 0) return 1;
    if (prec <= 0) prec = 17;

    size_t nb = s->count;
    size_t nk = strap_nk(s);
    if (layout == STRAP_LAYOUT_COLUMN && strap_check_uniform(s)) return 1;

    FILE *f = fopen(file_path, "w");
    if (!f) {
        fprintf(stderr, "[ERROR] 写不了 %s\n", file_path);
        return 1;
    }
    if (header) fprintf(f, "# %s\n", header);
    fprintf(f, "# nk = %zu, nb = %zu\n", nk, nb);

    if (layout == STRAP_LAYOUT_BLOCK) {
        for (size_t b = 0; b < nb; ++b) {
            if (b) fprintf(f, "\n");
            const Band *p = &s->bands[b];
            for (size_t i = 0; i < p->count; ++i)
                fprintf(f, "%.10g %.*g\n", p->keys[i], prec, p->values[i]);
        }
    } else {
        const Band *ref = &s->bands[0];
        for (size_t i = 0; i < nk; ++i) {
            fprintf(f, "%.10g", ref->keys[i]);
            for (size_t b = 0; b < nb; ++b)
                fprintf(f, " %.*g", prec, s->bands[b].values[i]);
            fprintf(f, "\n");
        }
    }
    fclose(f);
    return 0;
}

//=====================================
// 处理
//=====================================

int strap_sort_by_energy(StrapBands *s)
{
    if (!s || s->count == 0) return 1;
    if (strap_check_uniform(s)) return 1;

    size_t nb = s->count, nk = strap_nk(s);
    double *col = malloc(sizeof *col * nb);
    if (!col) {
        fprintf(stderr, "[ERROR] 内存不足\n");
        return 1;
    }
    for (size_t i = 0; i < nk; ++i) {
        for (size_t b = 0; b < nb; ++b) col[b] = s->bands[b].values[i];
        qsort(col, nb, sizeof *col, cmp_double);
        for (size_t b = 0; b < nb; ++b) s->bands[b].values[i] = col[b];
    }
    free(col);
    return 0;
}

int strap_find_anchor(const StrapBands *s, size_t *anchor)
{
    if (!s || !anchor || s->count < 2) return 1;
    if (strap_check_uniform(s)) return 1;

    size_t nb = s->count, nk = strap_nk(s);
    if (nk < 3) {
        fprintf(stderr, "[ERROR] strap_find_anchor: 只有 %zu 个 k 点，没法外推\n", nk);
        return 1;
    }
    double *col = malloc(sizeof *col * nb);
    if (!col) {
        fprintf(stderr, "[ERROR] 内存不足\n");
        return 1;
    }
    size_t best_i = nk / 2;
    double best_gap = -1;
    for (size_t i = 0; i < nk; ++i) {
        for (size_t b = 0; b < nb; ++b) col[b] = s->bands[b].values[i];
        qsort(col, nb, sizeof *col, cmp_double);
        double gap = col[1] - col[0];
        for (size_t b = 2; b < nb; ++b) {
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
    if (best_i + 1 >= nk) best_i = nk - 2;
    *anchor = best_i;
    return 0;
}

// 匈牙利算法（O(n³)）：求 n×n 代价矩阵的最优指派，assign[i] = 第 i 行分到的列
static int hungarian(size_t n, const double *cost, size_t *assign)
{
    double *u = calloc(n + 1, sizeof *u);
    double *v = calloc(n + 1, sizeof *v);
    size_t *p = calloc(n + 1, sizeof *p);
    size_t *way = calloc(n + 1, sizeof *way);
    double *minv = malloc(sizeof *minv * (n + 1));
    char *used = malloc(n + 1);
    if (!u || !v || !p || !way || !minv || !used) {
        fprintf(stderr, "[ERROR] hungarian: 内存不足\n");
        free(u); free(v); free(p); free(way); free(minv); free(used);
        return 1;
    }

    for (size_t i = 1; i <= n; ++i) {
        p[0] = i;
        size_t j0 = 0;
        for (size_t j = 0; j <= n; ++j) {
            minv[j] = 1e300;
            used[j] = 0;
        }
        do {
            used[j0] = 1;
            size_t i0 = p[j0], j1 = 0;
            double delta = 1e300;
            for (size_t j = 1; j <= n; ++j) {
                if (used[j]) continue;
                double cur = cost[(i0 - 1) * n + (j - 1)] - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }
            for (size_t j = 0; j <= n; ++j) {
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
            size_t j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }

    for (size_t j = 1; j <= n; ++j) assign[p[j] - 1] = j - 1;
    free(u); free(v); free(p); free(way); free(minv); free(used);
    return 0;
}

int strap_reconnect(StrapBands *s, size_t anchor)
{
    if (!s || s->count < 2) return 1;
    if (strap_check_uniform(s)) return 1;

    size_t nb = s->count, nk = strap_nk(s);
    if (nk < 3 || anchor == 0 || anchor + 1 >= nk) {
        fprintf(stderr, "[ERROR] strap_reconnect: 锚点 %zu 越界（需要 0 < anchor < %zu）\n",
                anchor, nk - 1);
        return 1;
    }

    size_t *asg = malloc(sizeof *asg * nk * nb);
    double *cost = malloc(sizeof *cost * nb * nb);
    size_t *perm = malloc(sizeof *perm * nb);
    double *cur = malloc(sizeof *cur * nb);
    double *prev = malloc(sizeof *prev * nb);
    double *pred = malloc(sizeof *pred * nb);
    double *now = malloc(sizeof *now * nb);
    double *out = malloc(sizeof *out * nk * nb);
    if (!asg || !cost || !perm || !cur || !prev || !pred || !now || !out) {
        fprintf(stderr, "[ERROR] strap_reconnect: 内存不足\n");
        free(asg); free(cost); free(perm); free(cur);
        free(prev); free(pred); free(now); free(out);
        return 1;
    }

    // 初始身份 = 各 k 点内部的排列顺序
    for (size_t i = 0; i < nk; ++i)
        for (size_t b = 0; b < nb; ++b) asg[i * nb + b] = b;

    // 自锚点向两侧逐点走：用前两步线性外推预测位置，再取最优指派
    for (int dir = 1; dir >= -1; dir -= 2) {
        for (size_t b = 0; b < nb; ++b) {
            cur[b] = s->bands[b].values[anchor];
            prev[b] = cur[b];
        }
        for (long i = (long)anchor + dir; i >= 0 && i < (long)nk; i += dir) {
            for (size_t b = 0; b < nb; ++b) {
                pred[b] = 2.0 * cur[b] - prev[b];
                now[b] = s->bands[b].values[i];
            }
            for (size_t b = 0; b < nb; ++b)
                for (size_t c = 0; c < nb; ++c) {
                    double d = now[c] - pred[b];
                    cost[b * nb + c] = d * d;
                }
            if (hungarian(nb, cost, perm)) {
                free(asg); free(cost); free(perm); free(cur);
                free(prev); free(pred); free(now); free(out);
                return 1;
            }
            for (size_t b = 0; b < nb; ++b) {
                asg[(size_t)i * nb + b] = perm[b];
                prev[b] = cur[b];
                cur[b] = now[perm[b]];
            }
        }
    }

    for (size_t i = 0; i < nk; ++i)
        for (size_t b = 0; b < nb; ++b)
            out[b * nk + i] = s->bands[asg[i * nb + b]].values[i];
    for (size_t i = 0; i < nk; ++i)
        for (size_t b = 0; b < nb; ++b)
            s->bands[b].values[i] = out[b * nk + i];

    free(asg); free(cost); free(perm); free(cur);
    free(prev); free(pred); free(now); free(out);
    return 0;
}

double strap_roughness(const StrapBands *s)
{
    if (!s) return 0.0;
    double sum = 0;
    for (size_t b = 0; b < s->count; ++b) {
        const Band *p = &s->bands[b];
        for (size_t i = 1; i + 1 < p->count; ++i) {
            double d2 = p->values[i + 1] - 2.0 * p->values[i] + p->values[i - 1];
            sum += d2 * d2;
        }
    }
    return sum;
}
