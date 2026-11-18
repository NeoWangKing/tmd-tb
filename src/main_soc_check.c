// 用 strap 库拼接"按能量排序"的能带，与真实物理通道对比，考察简并时的表现
//
//   data/band_soc.dat       ：含 SOC，每个 k 点按能量从小到大输出（身份已丢）
//   data/band_soc_spin.dat  ：同一次计算，按物理通道排（Up1-3 / Dn1-3），当参考
//
// 判据不看能量数值，看身份标签：把参考通道的劈裂按 α 缩放后排序，每个位置记下
// 它来自哪条通道；拼接完回看每条输出带的标签序列，标签中途变了说明身份接错。
// 切换点若正好落在简并处（间距 < 1e-4 eV），只是退化点上的任意挑选，不算错。
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "strap.h"

#define K_POINT   1.244017      // 阶段一路径上的 K 点
#define AMBIG_TOL 1e-4          // 小于这个间距视为简并，切换不算错
#define SAME_TOL  1e-9          // 判定"同一个值"的容差

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static int cmp_pair(const void *a, const void *b)
{
    double x = ((const double *)a)[0], y = ((const double *)b)[0];
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

// 建 nb 条空带，每条预留 nk 个点（strap 是 push 式接口）
static int make_bands(StrapBands *s, size_t nk, size_t nb)
{
    for (size_t b = 0; b < nb; ++b) {
        Band *band = strap_new_band(s);
        if (!band) return 1;
        if (band_reserve(band, nk)) return 1;
        band->count = nk;
    }
    return 0;
}

// 对一档劈裂强度做完整评估；返回认错身份的 k 点数
static size_t evaluate(const StrapBands *chan, size_t nb, size_t nk, size_t anchor,
                       double alpha, int verbose, int detail)
{
    double *val    = malloc(sizeof *val * nk * nb);
    double *sorted = malloc(sizeof *sorted * nk * nb);
    size_t *label  = malloc(sizeof *label * nk * nb);
    double *pair   = malloc(sizeof *pair * nb * 2);
    double *gap    = malloc(sizeof *gap * nk);
    double *gs     = malloc(sizeof *gs * nk);
    if (!val || !sorted || !label || !pair || !gap || !gs) {
        fprintf(stderr, "[ERROR] 内存不足\n");
        exit(1);
    }

    // 按 α 缩放劈裂（围绕每个 k 点上各通道的平均能量）
    for (size_t i = 0; i < nk; ++i) {
        double mean = 0;
        for (size_t c = 0; c < nb; ++c) mean += chan->bands[c].values[i];
        mean /= (double)nb;
        for (size_t c = 0; c < nb; ++c)
            val[c * nk + i] = mean + alpha * (chan->bands[c].values[i] - mean);
    }

    // 每个 k 点排序，记录每个排序位置来自哪条通道
    for (size_t i = 0; i < nk; ++i) {
        for (size_t c = 0; c < nb; ++c) {
            pair[c * 2 + 0] = val[c * nk + i];
            pair[c * 2 + 1] = (double)c;
        }
        qsort(pair, nb, sizeof *pair * 2, cmp_pair);
        for (size_t b = 0; b < nb; ++b) {
            sorted[i * nb + b] = pair[b * 2 + 0];
            label[i * nb + b]  = (size_t)pair[b * 2 + 1];
        }
        gap[i] = sorted[i * nb + 1] - sorted[i * nb + 0];
        for (size_t b = 2; b < nb; ++b) {
            double d = sorted[i * nb + b] - sorted[i * nb + b - 1];
            if (d < gap[i]) gap[i] = d;
        }
    }

    // 装进 StrapBands，按"已排序"的数据拼接
    StrapBands t;
    strap_init(&t);
    if (make_bands(&t, nk, nb)) exit(1);
    for (size_t b = 0; b < nb; ++b)
        for (size_t i = 0; i < nk; ++i) {
            t.bands[b].keys[i]   = chan->bands[0].keys[i];
            t.bands[b].values[i] = sorted[i * nb + b];
        }
    double rough0 = strap_roughness(&t);
    if (strap_reconnect(&t, anchor)) exit(1);
    double rough1 = strap_roughness(&t);

    // 逐带追踪标签
    size_t switches = 0, real_err = 0;
    double err_gap_max = 0;
    size_t *track = malloc(sizeof *track * nb);
    for (size_t b = 0; b < nb; ++b) {
        size_t prev = (size_t)-1, major = 0, major_cnt = 0;
        for (size_t i = 0; i < nk; ++i) {
            // 在排好序的列里找到该值的位置，取它的通道号
            size_t lab = (size_t)-1;
            for (size_t j = 0; j < nb; ++j) {
                if (fabs(sorted[i * nb + j] - t.bands[b].values[i]) < SAME_TOL) {
                    if (lab == (size_t)-1) lab = label[i * nb + j];
                    // 简并处多个候选：优先沿用上一点的标签，避免假切换
                    if (label[i * nb + j] == prev) { lab = prev; break; }
                }
            }
            if (prev != (size_t)-1 && lab != prev) {
                ++switches;
                if (gap[i] > AMBIG_TOL) {
                    ++real_err;
                    if (gap[i] > err_gap_max) err_gap_max = gap[i];
                    if (detail)
                        printf("      k=%3zu (k=%.4f)：间距 %.3f meV，通道 %zu -> %zu\n",
                               i, chan->bands[0].keys[i], gap[i] * 1000.0, prev + 1, lab + 1);
                }
            }
            prev = lab;
        }
        (void)major; (void)major_cnt;
        track[b] = prev;      // 记录一下最后落在哪条通道（仅用于打印）
    }

    memcpy(gs, gap, sizeof *gs * nk);
    qsort(gs, nk, sizeof *gs, cmp_double);
    if (verbose) {
        printf("  %-6.3f %14.3f %12zu %12zu %16.3f %10.3f\n",
               alpha, gs[nk / 2] * 1000.0, switches, real_err,
               err_gap_max * 1000.0, rough1 / (rough0 > 0 ? rough0 : 1.0));
    } else {
        printf("α=%.3f：标签切换 %zu 次，其中真错 %zu 次，粗糙度 %.6g -> %.6g\n",
               alpha, switches, real_err, rough0, rough1);
    }
    free(track);
    free(val); free(sorted); free(label); free(pair); free(gap); free(gs);
    strap_free(&t);
    return real_err;
}


// 把"参考通道 + 拼接结果 + 每条拼接带在各 k 点跟的是哪条通道"导出成一张表
static int export_compare(const char *path, const StrapBands *ref, const StrapBands *work,
                          size_t nb, size_t nk, size_t anchor)
{
    StrapBands t;
    strap_init(&t);
    if (make_bands(&t, nk, nb)) return 1;
    for (size_t b = 0; b < nb; ++b)
        for (size_t i = 0; i < nk; ++i) {
            t.bands[b].keys[i]   = work->bands[0].keys[i];
            t.bands[b].values[i] = work->bands[b].values[i];
        }
    strap_sort_by_energy(&t);
    if (strap_reconnect(&t, anchor)) return 1;

    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "[ERROR] 写不了 %s\n", path);
        strap_free(&t);
        return 1;
    }
    // 按锚点处的通道归属给拼接带排序，这样第 c 列对应第 c 条参考通道，
    // 两个面板里同一种颜色就是同一条物理通道
    // 锚点处可能有两条通道恰好简并，所以必须做一对一指派，不能"值相等就配"
    size_t *order = malloc(sizeof *order * nb);
    char *used = calloc(nb, 1);
    if (!order || !used) { fclose(f); free(order); free(used); strap_free(&t); return 1; }
    for (size_t c = 0; c < nb; ++c) {
        size_t best = (size_t)-1;
        double best_d = 1e300;
        for (size_t b = 0; b < nb; ++b) {
            if (used[b]) continue;
            double d = fabs(t.bands[b].values[anchor] - ref->bands[c].values[anchor]);
            if (d < best_d) { best_d = d; best = b; }
        }
        order[c] = best;
        used[best] = 1;
    }
    free(used);

    fprintf(f, "# col1: k\n");
    fprintf(f, "# col2-7: 参考通道 Up1 Up2 Up3 Dn1 Dn2 Dn3\n");
    fprintf(f, "# col8-13: 拼接结果，按锚点处的通道归属重排，第 c 列对应第 c 条参考通道\n");
    fprintf(f, "# col14-19: 每条拼接带在该 k 跟的是哪条参考通道（1-6）\n");
    for (size_t i = 0; i < nk; ++i) {
        fprintf(f, "%.6f", ref->bands[0].keys[i]);
        for (size_t c = 0; c < nb; ++c) fprintf(f, " %10.5f", ref->bands[c].values[i]);
        for (size_t c = 0; c < nb; ++c) fprintf(f, " %10.5f", t.bands[order[c]].values[i]);
        for (size_t b = 0; b < nb; ++b) {
            size_t lab = 1;
            for (size_t c = 0; c < nb; ++c)
                if (fabs(ref->bands[c].values[i] - t.bands[b].values[i]) < 1e-9) { lab = c + 1; break; }
            fprintf(f, " %2zu", lab);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    free(order);
    strap_free(&t);
    printf("对比数据已写出：%s\n", path);
    return 0;
}

int main(int argc, char **argv)
{
    const char *f_sorted = (argc > 1) ? argv[1] : "data/band_soc.dat";
    const char *f_ref    = (argc > 2) ? argv[2] : "data/band_soc_spin.dat";

    StrapBands ref, work;
    strap_init(&ref);
    strap_init(&work);
    if (strap_read(f_ref, &ref)) return 1;
    if (strap_read(f_sorted, &work)) return 1;

    size_t nb = strap_band_count(&ref), nk = strap_sample_count(&ref);
    if (strap_band_count(&work) != nb || strap_sample_count(&work) != nk) {
        fprintf(stderr, "[ERROR] 两个文件结构不一致\n");
        return 1;
    }

    size_t anchor = 0;
    for (size_t i = 1; i < nk; ++i)
        if (fabs(ref.bands[0].keys[i] - K_POINT) < fabs(ref.bands[0].keys[anchor] - K_POINT))
            anchor = i;

    printf("参考文件：%s（%zu 条通道，%zu 个 k 点）\n", f_ref, nb, nk);
    printf("待拼接：  %s，锚点取 K 点（第 %zu 个 k，k=%.6f）\n\n",
           f_sorted, anchor, ref.bands[0].keys[anchor]);

    // 先看输入数据本身的简并程度
    double *col = malloc(sizeof *col * nb);
    size_t n_exact = 0, n_tiny = 0;
    for (size_t i = 0; i < nk; ++i) {
        for (size_t b = 0; b < nb; ++b) col[b] = work.bands[b].values[i];
        qsort(col, nb, sizeof *col, cmp_double);
        double g = col[1] - col[0];
        for (size_t b = 2; b < nb; ++b) {
            double d = col[b] - col[b - 1];
            if (d < g) g = d;
        }
        if (g < 1e-6) ++n_exact;
        if (g < 1e-3) ++n_tiny;
    }
    printf("输入数据的简并情况（每个 k 点相邻能级的最小间距）：\n");
    printf("  严格简并 < 1e-6 eV：%zu / %zu 个 k 点\n", n_exact, nk);
    printf("  近似简并 < 1e-3 eV：%zu / %zu 个 k 点\n\n", n_tiny, nk);

    printf("劈裂扫描（每档都重新排序 + 拼接，再看身份标签有没有接错）：\n");
    printf("  α     最小间距中位数(meV)   标签切换   真错   真错处最大间距(meV)  粗糙度比\n");
    const double alphas[] = {1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.0};
    for (size_t ia = 0; ia < sizeof alphas / sizeof alphas[0]; ++ia)
        evaluate(&ref, nb, nk, anchor, alphas[ia], 1, 0);

    if (export_compare("data/soc_stitch_compare.dat", &ref, &work, nb, nk, anchor)) return 1;

    printf("\nα=1（即真实劈裂，等价于直接拼 data/band_soc.dat）里认错身份的位置：\n");
    evaluate(&ref, nb, nk, anchor, 1.0, 0, 1);

    free(col);
    strap_free(&ref);
    strap_free(&work);
    return 0;
}
