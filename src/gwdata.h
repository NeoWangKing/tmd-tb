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

// 去掉相邻重复的 k（TB 输出会把段端点打印两次）
int dedup(int n, const double *k, const double *e, int nband, double *k2, double *e2);

#endif
