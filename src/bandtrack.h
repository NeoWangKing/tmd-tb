#ifndef BANDTRACK_H_
#define BANDTRACK_H_

// 最小曲率能带拼接（band stitching）
//
// 要解决的问题：wannier90 的能带文件、以及 TB 模型自己的本征值，每个 k 点都是按
// 能量从小到大排好序输出的。两条物理能带一旦交叉，排序结果里它们就互换了身份 ——
// 画出来是两条互相"弹开"的折线，而不是交叉的两条直线。要把它们接回物理能带，
// 就得定出相邻 k 点之间"排序能级"和"物理能带"的对应关系。
//
// 判据（最小曲率）：从 anchor 点（各带分得开、顺序无歧义的地方）出发向两侧逐点走。
// 用前两步的位置线性外推（保持二阶差分最小 ⇔ 曲线最不弯）预测下一列各带的位置，
// 再在候选能级的全部排列里挑一个让 Σ(位置 − 预测)² 最小的分配。对光滑能带来说，
// 交叉点处自然会把身份接回去。
//
// 用法（cand/out 均为"带优先"排布：第 b 条带在 k_i 的值 = x[b*nk + i]）：
//
//     double *cand = malloc(sizeof(double) * ncand * nk);   // 每个 k 点内部已排序
//     double *out  = malloc(sizeof(double) * ncand * nk);
//     ...
//     if (band_track(nk, ncand, cand, anchor, out)) { /* 失败 */ }
//
// 只依赖 <stdlib.h>，可以整个抄到别的项目里用。

#define BANDTRACK_NCAND_MAX 8

// 把 ncand 条"按能量排序"的能级接成 ncand 条物理能带。
//   cand[c*nk+i] : 第 c 个候选能级在 k_i 的能量
//   anchor       : 锚点下标（不能取端点，需要两侧都有邻居）
//   out[b*nk+i]  : 拼接结果，带编号按锚点处的能量顺序
// 返回 0 成功，1 参数不合法或内存不足。
int band_track(int nk, int ncand, const double *cand, int anchor, double *out);

#endif
