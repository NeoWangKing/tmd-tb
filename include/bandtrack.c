// 最小曲率能带拼接，见 bandtrack.h
#include <stdlib.h>
#include "bandtrack.h"

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

int band_track(int nk, int ncand, const double *cand, int anchor, double *out)
{
    if (ncand < 2 || ncand > BANDTRACK_NCAND_MAX || anchor <= 0 || anchor >= nk - 1) return 1;

    int *asg = malloc(sizeof(int) * nk * ncand);
    if (!asg) return 1;
    for (int i = 0; i < nk; ++i)
        for (int b = 0; b < ncand; ++b) asg[i*ncand + b] = b;   // 初始 = 排序顺序

    // 从 anchor 向两侧逐点分配：用"斜率外推"预测下一时刻的位置，
    // 再在所有排列里选距离平方和最小的那个 —— 等价于保持斜率连续。
    for (int step = +1; step >= -1; step -= 2) {
        double prev[BANDTRACK_NCAND_MAX], curp[BANDTRACK_NCAND_MAX];
        double cnow[BANDTRACK_NCAND_MAX], pred[BANDTRACK_NCAND_MAX];
        for (int b = 0; b < ncand; ++b) {
            curp[b] = cand[b*nk + anchor];
            prev[b] = curp[b];
        }
        for (int i = anchor + step; i >= 0 && i < nk; i += step) {
            for (int b = 0; b < ncand; ++b) {
                pred[b] = 2*curp[b] - prev[b];
                cnow[b] = cand[b*nk + i];
            }
            int perm[BANDTRACK_NCAND_MAX], used[BANDTRACK_NCAND_MAX] = {0};
            int bestperm[BANDTRACK_NCAND_MAX];
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
