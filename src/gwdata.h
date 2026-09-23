#ifndef GWDATA_H_
#define GWDATA_H_

// GW 参考能带（wannier90 格式：按"带"分块，块间空行）
// 与 TB 能带表（每行 k + nband 个能量）的读取

#define GW_NBAND_MAX 64
#define NK_MAX 100000

typedef struct {
    int nband;
    int nk;
    double *k;      // [nk]
    double *E;      // [nband*nk]，按带分块
} GwData;

int gw_read(const char *path, GwData *g);
int tb_read(const char *path, int nband, int *nk_out, double **k_out, double **e_out);

// 读 wannier90 的 .kpt（第一行是点数，之后每行 f1 f2 f3 weight），
// 输出 TB 单位（1/a）的笛卡尔 k
int kpt_read(const char *path, int *nk_out, double **kx_out, double **ky_out);

// 把"按能量排序"的能级重连成物理能带（GW 参考和 TB 模型两边都用）。
// 判据：从 anchor 点（各带分得开、顺序无歧义）向两侧逐点走，用斜率外推预测下一时刻
// 的位置，再在所有排列里选距离平方和最小的 —— 等价于"保持斜率连续"，也就是让
// 重连后的曲线尽量光滑（交叉点处排序曲线会互换身份产生折点，这样能把它们接回去）。
//   cand[c*nk+i]: 第 c 个候选能级在 k_i 的能量
//   out[b*nk+i] : 重连后第 b 条曲线（编号按 anchor 点的能量顺序）
// 返回 0 成功。
int band_track(int nk, int ncand, const double *cand, int anchor, double *out);

// 去掉相邻重复的 k（TB 输出会把段端点打印两次）
int dedup(int n, const double *k, const double *e, int nband, double *k2, double *e2);

#endif
