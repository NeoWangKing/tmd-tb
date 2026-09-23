#ifndef KPATH_H_
#define KPATH_H_

#include <stddef.h>
#include "config.h"
#include "tb.h"

// 一段直线路径（分数坐标）与采样步数
typedef struct {
    const double *start;
    const double *end;
    int steps;
} Segment;

// 横坐标单位与标签：PATH_MODE=0 用 2π/a，=1 用 Å⁻¹
#if PATH_MODE == 0
  #define K_AXIS_UNIT    (2*PI)
  #define K_AXIS_LABEL   "k-path (fractional)"
  #define K_AXIS_LABEL_S "k-path (frac)"
#else
  #define K_AXIS_UNIT    (A_ANG)
  #define K_AXIS_LABEL   "k-path (Å^-1)"
  #define K_AXIS_LABEL_S "k-path (Å^-1)"
#endif

// 高对称点（分数坐标）
extern const double KPATH_GAMMA[2];
extern const double KPATH_M[2];
extern const double KPATH_K[2];

// 当前 PATH_MODE 对应的路径段列表
const Segment *kpath_segments(size_t *n);

// 一段路径的 |Δk|（单位 1/a）
double kpath_seg_len(const Segment *seg);

#endif
