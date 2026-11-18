#include <math.h>
#include "kpath.h"

// 高对称点路径定义
const double KPATH_GAMMA[2] = {0.0, 0.0};          // Gamma 点
const double KPATH_M[2]     = {0.5, 0.5};          // M 点
const double KPATH_K[2]     = {2.0/3.0, 1.0/3.0};  // K 点
#if PATH_MODE == 0
static const double Mp[2]    = {0.5, 0.0};          // M' 点
static const double Kp[2]    = {1.0/3.0, -1.0/3.0}; // K' 点
#endif

static Segment path[] = {
#if PATH_MODE == 0
    {KPATH_M, KPATH_GAMMA, 40},
    {KPATH_GAMMA, KPATH_K, 40},
    {KPATH_K, Mp, 40},
    {Mp, Kp, 40}
#else
    {KPATH_GAMMA, KPATH_K, 120},
    {KPATH_K, KPATH_M, 60},
    {KPATH_M, KPATH_GAMMA, 104}
#endif
};
static const size_t num_segments = sizeof(path)/sizeof(Segment);

const Segment *kpath_segments(size_t *n)
{
    *n = num_segments;
    return path;
}

double kpath_seg_len(const Segment *seg)
{
    // 与原来的算法一致：先在分数坐标下求差，再换算成笛卡尔坐标
    double df1 = seg->end[0] - seg->start[0];
    double df2 = seg->end[1] - seg->start[1];
    double kx, ky;
    tb_frac_to_cart_k(df1, df2, &kx, &ky);
    return sqrt(kx*kx + ky*ky);
}
