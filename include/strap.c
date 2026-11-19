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

int band_reserve(Band *band, size_t nk)
{
    if (!band) return 1;
    if (nk <= band->capacity) return 0;

    size_t capacity = band->capacity ? band->capacity : BAND_INIT_CAP;
    while (capacity < nk) capacity *= 2;

    // 先接临时指针：realloc 失败时原指针不能丢
    double *tmp = realloc(band->keys, capacity * sizeof *tmp);
    if (!tmp) {
        fprintf(stderr, "[ERROR] band_reserve: could not expand key to %zu\n", capacity);
        return 1;
    }
    band->keys = tmp;
    tmp = realloc(band->values, capacity * sizeof *tmp);
    if (!tmp) {
        fprintf(stderr, "[ERROR] band_reserve: could not expand values to %zu\n", capacity);
        return 1;
    }
    band->values = tmp;
    band->capacity = capacity;
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

size_t strap_band_count(const StrapBands *s)
{
    return s ? s->count : 0;
}

size_t strap_sample_count(const StrapBands *s)
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
        fprintf(stderr, "[ERROR] strap_reserve: could not expand bands to %zu\n", nc);
        return 1;
    }
    s->bands = tmp;
    s->capacity = nc;
    return 0;
}

Band *strap_new_band(StrapBands *s)
{
    if (!s) return NULL;
    if (strap_reserve(s, s->count + 1)) return NULL;
    band_init(&s->bands[s->count]);
    return &s->bands[s->count++];
}

int strap_check_uniform(const StrapBands *strap)
{
    if (!strap || strap->count == 0) {
        fprintf(stderr, "[ERROR] check_uniform: there is not any band\n");
        return 1;
    }
    const Band *ref = &strap->bands[0];
    if (ref->count == 0) {
        fprintf(stderr, "[ERROR] check_uniform: the band 1 has no point\n");
        return 1;
    }
    for (size_t i = 1; i < strap->count; ++i) {
        const Band *band = &strap->bands[i];
        if (band->count != ref->count) {
            fprintf(stderr, "[ERROR] check_uniform: the band %zu has %zu points (%zu are expected)\n",
                    i + 1, band->count, ref->count);
            return 1;
        }
        for (size_t j = 0; j < ref->count; ++j) {
            double tol = 1e-8 * (1.0 + fabs(ref->keys[j]));
            if (fabs(band->keys[j] - ref->keys[j]) > tol) {
                fprintf(stderr,
                        "[ERROR] check_uniform: the point %zu of the band %zu (%.10g) has a different key than the band 1 (%.10g)\n",
                        j, i + 1, band->keys[j], ref->keys[j]);
                return 1;
            }
        }
    }
    return 0;
}

// 块式：每块 "k E"，块与块之间用空行隔开，一块就是一条带
static int read_block(const char *file_path, StrapBands *strap)
{
    FILE *f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] could not open %s\n", file_path);
        return 1;
    }

    char line[BUF_SIZE];
    Band *band = NULL;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line)) { band = NULL; continue; }
        if (line[0] == '#') continue;

        double key = 0, value = 0;
        if (sscanf(line, "%lf %lf", &key, &value) != 2) continue;
        if (!band) {
            band = strap_new_band(strap);
            if (!band) {
                fclose(f);
                return 1;
            }
        }
        if (band_push(band, key, value)) {
            fclose(f);
            return 1;
        }
    }
    fclose(f);

    if (strap->count == 0) {
        fprintf(stderr, "[ERROR] %s: there is no data\n", file_path);
        return 1;
    }
    return strap_check_uniform(strap);
}

// 列式：每行 "k E1 E2 ... En"，第 c 个能量攒成第 c 条带
static int read_column(const char *file_path, StrapBands *strap)
{
    FILE *f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] could not open %s\n", file_path);
        return 1;
    }

    char line[BUF_SIZE];
    int column_count = -1;
    size_t key_count = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        int data_count = count_datas(line);
        if (column_count < 0) column_count = data_count;
        else if (data_count != column_count) {
            fprintf(stderr, "[ERROR] %s: there are %d datas in the line %zu (%d datas are expected)\n",
                    file_path, data_count, key_count + 1, column_count);
            fclose(f);
            return 1;
        }
        ++key_count;
    }
    fclose(f);

    if (key_count == 0 || column_count < 2) {
        fprintf(stderr, "[ERROR] %s: there is no data to read (nk=%zu ncol=%d)\n", file_path, key_count, column_count);
        return 1;
    }

    f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] could not open %s\n", file_path);
        return 1;
    }
    for (int i = 1; i < column_count; ++i) {
        if (!strap_new_band(strap)) {
            fclose(f);
            return 1;
        }
    }
    for (size_t i = 0; i < strap->count; ++i) {
        if (band_reserve(&strap->bands[i], key_count)) {
            fclose(f);
            return 1;
        }
    }

    size_t i = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        if (i >= key_count) break;

        const char *p = line;
        char *end = NULL;
        double key = strtod(p, &end);
        if (end == p) continue;
        p = end;
        for (size_t j = 0; j < strap->count; ++j) {
            double value = strtod(p, &end);
            if (end == p) {
                fprintf(stderr, "[ERROR] %s: there are only %zu values in the line %zu\n", file_path, j, i + 1);
                fclose(f);
                return 1;
            }
            p = end;
            if (band_push(&strap->bands[j], key, value)) {
                fclose(f);
                return 1;
            }
        }
        ++i;
    }
    fclose(f);
    return strap_check_uniform(strap);
}

int strap_read(const char *file_path, StrapBands *strap)
{
    if (!file_path || !strap) return 1;
    strap_free(strap);

    FILE *f = fopen(file_path, "r");
    if (!f) {
        fprintf(stderr, "[ERROR] could not open %s\n", file_path);
        return 1;
    }
    char line[BUF_SIZE];
    int column_count = 0;
    while (fgets(line, sizeof line, f)) {
        if (line_is_blank(line) || line[0] == '#') continue;
        column_count = count_datas(line);
        break;
    }
    fclose(f);

    if (column_count <= 0) {
        fprintf(stderr, "[ERROR] there is no data in %s\n", file_path);
        return 1;
    }
    return (column_count > 2) ? read_column(file_path, strap) : read_block(file_path, strap);
}

int strap_write(const char *file_path, const StrapBands *strap, int layout,
                int sig, const char *header)
{
    if (!file_path || !strap || strap->count == 0) return 1;
    if (sig <= 0) sig = 17;

    size_t band_count = strap->count;
    size_t sample_count = strap_sample_count(strap);
    if (layout == STRAP_LAYOUT_COLUMN && strap_check_uniform(strap)) return 1;

    FILE *f = fopen(file_path, "w");
    if (!f) {
        fprintf(stderr, "[ERROR] could not write %s\n", file_path);
        return 1;
    }
    if (header) fprintf(f, "# %s\n", header);
    fprintf(f, "# sample_count = %zu, band_count = %zu\n", sample_count, band_count);

    switch (layout) {
        case STRAP_LAYOUT_BLOCK : 
            {
                for (size_t b = 0; b < band_count; ++b) {
                    if (b) fprintf(f, "\n");
                    const Band *p = &strap->bands[b];
                    for (size_t i = 0; i < p->count; ++i)
                        fprintf(f, "%.*g %.*g\n", sig, p->keys[i], sig, p->values[i]);
                }
                break;
            }
        case STRAP_LAYOUT_COLUMN :
            {
                const Band *ref = &strap->bands[0];
                for (size_t i = 0; i < sample_count; ++i) {
                    fprintf(f, "%.*g", sig, ref->keys[i]);
                    for (size_t b = 0; b < band_count; ++b)
                        fprintf(f, " %.*g", sig, strap->bands[b].values[i]);
                    fprintf(f, "\n");
                }
                break;
            }
        default : 
            {
                fprintf(stderr, "[ERROR] unknown layout %d\n", layout);
                fclose(f);
                return 1;
            }
    }
    fclose(f);
    return 0;
}

int strap_sort_by_energy(StrapBands *strap)
{
    if (!strap || strap->count == 0) return 1;
    if (strap_check_uniform(strap)) return 1;

    size_t band_count = strap->count;
    size_t sample_count = strap_sample_count(strap);
    double *column = malloc(sizeof *column * band_count);
    if (!column) {
        fprintf(stderr, "[ERROR] buy more ram lol!\n");
        return 1;
    }
    for (size_t i = 0; i < sample_count; ++i) {
        for (size_t j = 0; j < band_count; ++j) column[j] = strap->bands[j].values[i];
        qsort(column, band_count, sizeof *column, cmp_double);
        for (size_t j = 0; j < band_count; ++j) strap->bands[j].values[i] = column[j];
    }
    free(column);
    return 0;
}

int strap_find_anchor(const StrapBands *strap, size_t *anchor)
{
    if (!strap || !anchor || strap->count < 2) return 1;
    if (strap_check_uniform(strap)) return 1;

    size_t band_count = strap->count, sample_count = strap_sample_count(strap);
    if (sample_count < 3) {
        fprintf(stderr, "[ERROR] strap_find_anchor: there are only %zu sample point\n", sample_count);
        return 1;
    }
    double *column = malloc(sizeof *column * band_count);
    if (!column) {
        fprintf(stderr, "[ERROR] buy more ram lol!\n");
        return 1;
    }
    size_t best_i = sample_count / 2;
    double best_gap = -1;
    for (size_t i = 0; i < sample_count; ++i) {
        for (size_t j = 0; j < band_count; ++j) column[j] = strap->bands[j].values[i];
        qsort(column, band_count, sizeof *column, cmp_double);
        double gap = column[1] - column[0];
        for (size_t b = 2; b < band_count; ++b) {
            double d = column[b] - column[b - 1];
            if (d < gap) gap = d;
        }
        if (gap > best_gap) {
            best_gap = gap;
            best_i = i;
        }
    }
    free(column);
    // 锚点两侧都要有邻居，端点处没法外推
    if (best_i == 0) best_i = 1;
    if (best_i + 1 >= sample_count) best_i = sample_count - 2;
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
        fprintf(stderr, "[ERROR] hungarian: buy more ram lol!\n");
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

int strap_reconnect(StrapBands *strap, size_t anchor)
{
    if (!strap || strap->count < 2) return 1;
    if (strap_check_uniform(strap)) return 1;

    size_t nb = strap->count, nk = strap_sample_count(strap);
    if (nk < 3 || anchor == 0 || anchor + 1 >= nk) {
        fprintf(stderr, "[ERROR] strap_reconnect: anchor out of edge %zu (expect 0 < anchor < %zu)\n",
                anchor, nk - 1);
        return 1;
    }

    size_t *asg  = malloc(sizeof *asg  * nb * nk);
    double *cost = malloc(sizeof *cost * nb * nb);
    size_t *perm = malloc(sizeof *perm * nb);
    double *cur  = malloc(sizeof *cur  * nb);
    double *prev = malloc(sizeof *prev * nb);
    double *pred = malloc(sizeof *pred * nb);
    double *now  = malloc(sizeof *now  * nb);
    double *out  = malloc(sizeof *out  * nb * nk);
    if (!asg || !cost || !perm || !cur || !prev || !pred || !now || !out) {
        fprintf(stderr, "[ERROR] strap_reconnect: buy more ram lol!\n");
        free(asg);  free(cost); free(perm); free(cur);
        free(prev); free(pred); free(now);  free(out);
        return 1;
    }

    // 初始身份 = 各 k 点内部的排列顺序
    for (size_t i = 0; i < nk; ++i)
        for (size_t b = 0; b < nb; ++b) asg[i * nb + b] = b;

    // 自锚点向两侧逐点走：用前两步线性外推预测位置，再取最优指派
    for (int dir = 1; dir >= -1; dir -= 2) {
        for (size_t b = 0; b < nb; ++b) {
            cur[b] = strap->bands[b].values[anchor];
            prev[b] = cur[b];
        }
        for (long i = (long)anchor + dir; i >= 0 && i < (long)nk; i += dir) {
            for (size_t b = 0; b < nb; ++b) {
                pred[b] = 2.0 * cur[b] - prev[b];
                now[b] = strap->bands[b].values[i];
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
            out[b * nk + i] = strap->bands[asg[i * nb + b]].values[i];
    for (size_t i = 0; i < nk; ++i)
        for (size_t b = 0; b < nb; ++b)
            strap->bands[b].values[i] = out[b * nk + i];

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

// 裸数组版入口：内部装一遍 StrapBands，用完就释放
int strap_reconnect_array(size_t nk, size_t nb, const double *values,
                          size_t anchor, double *out)
{
    if (!values || !out) {
        fprintf(stderr, "[ERROR] reconnect_array: null pointer\n");
        return 1;
    }
    if (nb < 2 || nk < 3) {
        fprintf(stderr, "[ERROR] reconnect_array: need nb >= 2 and nk >= 3 (got %zu, %zu)\n", nb, nk);
        return 1;
    }

    StrapBands s;
    strap_init(&s);
    for (size_t b = 0; b < nb; ++b) {
        Band *band = strap_new_band(&s);
        if (!band) {
            fprintf(stderr, "[ERROR] reconnect_array: cannot allocate band %zu\n", b);
            strap_free(&s);
            return 1;
        }
        for (size_t i = 0; i < nk; ++i) {
            if (band_push(band, (double)i, values[b * nk + i])) {
                fprintf(stderr, "[ERROR] reconnect_array: cannot push point %zu of band %zu\n", i, b);
                strap_free(&s);
                return 1;
            }
        }
    }

    if (strap_reconnect(&s, anchor)) {
        strap_free(&s);
        return 1;
    }

    for (size_t b = 0; b < nb; ++b)
        for (size_t i = 0; i < nk; ++i)
            out[b * nk + i] = s.bands[b].values[i];

    strap_free(&s);
    return 0;
}
