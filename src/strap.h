#ifndef STRAP_H_
#define STRAP_H_

// 能带数据表的读写、排序与最小曲率拼接。
//
// 要解决的问题：出能带的程序（wannier90、紧束缚对角化等）在每个 k 点都把
// 本征值按能量排序输出，两条物理能带一旦交叉，排序结果里它们就互换了身份，
// 画出来是两条互相"弹开"的折线，而不是交叉的两条直线。
// 本库只做纯数据处理，不含物理：读表、按能量重排、把身份接回去、写表。
//
// 排布约定：E 为"带优先"，第 b 条带在第 i 个 k 点的值 = E[b*nk + i]。
// 候选能级在每个 k 点内部的排列可以任意，但各 k 点之间必须一致
// （同一列位置对应同一来源，wannier90 输出的这个顺序就是能量升序）。
//
// 用法：
//     StrapBands d;
//     int anchor;
//     strap_init(&d);
//     if (strap_read("wannier90_band.dat", &d)) return 1;
//     if (strap_find_anchor(&d, &anchor)) return 1;
//     if (strap_reconnect(&d, anchor)) return 1;
//     strap_write("out.dat", &d, STRAP_LAYOUT_BLOCK, 6, NULL);
//     strap_free(&d);
//
// 只依赖 libc，不引用本项目任何头文件，可以整对拷到别的项目里用。

#define STRAP_LAYOUT_BLOCK  0   // 每块 "k E"，块间空行（wannier90 band.dat）
#define STRAP_LAYOUT_COLUMN 1   // 每行 "k E1 E2 ... En"（VASP / 紧束缚输出）

typedef struct
{
    int     nk;     // k 点个数
    int     nband;  // 能带条数
    double *k;      // [nk]
    double *E;      // [nband*nk]
} StrapBands;

// 清零。StrapBands 必须先经 strap_init 或写成 = {0}，
// 否则 strap_alloc / strap_free 会去 free 未初始化的指针。
void strap_init(StrapBands *s);

// 分配 nk×nband 的表（旧内容会被释放）
int  strap_alloc(StrapBands *s, int nk, int nband);
void strap_free(StrapBands *s);

// 读表：按第一条数据行的列数自动区分块式（2 列）与列式（>2 列）
int  strap_read(const char *path, StrapBands *s);

// 写表：header 可为 NULL；prec 是 %g 的有效位数（<=0 时取 17，可无损往返）
int  strap_write(const char *path, const StrapBands *s, int layout, int prec,
                 const char *header);

// 每个 k 点内部按能量升序重排（原地）
int  strap_sort_by_energy(StrapBands *s);

// 找一个"各带分得最开"的 k 点当锚点：取相邻能级最小间距最大的那一列
int  strap_find_anchor(const StrapBands *s, int *anchor);

// 最小曲率拼接：从锚点向两侧逐点走，用前两步线性外推预测位置，
// 再在全部排列里选让 Σ(位置-预测)² 最小的分配（匈牙利算法，O(nband³)）。
// 必须 0 < anchor < nk-1。
int  strap_reconnect(StrapBands *s, int anchor);

// 总粗糙度 Σ(二阶差分)²，用来判断拼接前后哪一份更像物理能带
double strap_roughness(const StrapBands *s);

#endif
