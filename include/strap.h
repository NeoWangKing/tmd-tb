#ifndef STRAP_H_
#define STRAP_H_

// 能带数据表的读写、排序与最小曲率拼接。
//
// 数据模型：一条带就是一串 (横坐标, 纵坐标) 点。能带表里横坐标是 k、纵坐标是能量，
// 但这个库本身不关心物理含义，任何 x-y 序列都能装。
//
// 于是 sort / reconnect 这种"逐 k 点跨带"的操作隐含一个前提：所有带共用同一条
// 横坐标轴。strap_check_uniform() 检查这一点，相关函数内部也会先查一遍，
// 不一致就报错，不会静默算错。
//
// 用法：
//     StrapBands d;
//     size_t anchor;
//     strap_init(&d);
//     if (strap_read("wannier90_band.dat", &d)) return 1;
//     if (strap_find_anchor(&d, &anchor)) return 1;
//     if (strap_reconnect(&d, anchor)) return 1;
//     strap_write("out.dat", &d, STRAP_LAYOUT_BLOCK, 6, NULL);
//     strap_free(&d);
//
// 只依赖 libc，不引用本项目任何头文件，可以整对拷到别的项目里用。

#define STRAP_LAYOUT_BLOCK  0   // 每条带一块 "k E"，块间空行（wannier90 band.dat）
#define STRAP_LAYOUT_COLUMN 1   // 每行一个 k："k E1 E2 ... En"（VASP / 紧束缚输出）

typedef struct
{
    double *keys;       // [capacity] 横坐标
    double *values;     // [capacity] 纵坐标
    size_t  count;      // 已有多少个点
    size_t  capacity;   // keys / values 各能装多少个
} Band;

typedef struct
{
    Band  *bands;       // [capacity]
    size_t count;       // 已有多少条带
    size_t capacity;    // bands 数组能装多少条
} StrapBands;

// --- 一条带 ---
void band_init(Band *b);
void band_free(Band *b);
int  band_reserve(Band *b, size_t cap);
int  band_push(Band *b, double key, double value);

// --- 一组带 ---
// StrapBands 必须先 strap_init（或写成 = {0}），否则 strap_free 会 free 垃圾指针。
void   strap_init(StrapBands *s);
void   strap_free(StrapBands *s);
size_t strap_nb(const StrapBands *s);
size_t strap_nk(const StrapBands *s);
// 追加一条空带并返回它的指针。该指针在下一次 strap_push_band 之前有效
// （bands 数组可能被 realloc 搬走）。
Band  *strap_push_band(StrapBands *s);

// 各带点数相同、横坐标逐个相同才返回 0
int strap_check_uniform(const StrapBands *s);

// --- 读写 ---
// 按首行列数自动区分：2 列 → 块式，超过 2 列 → 列式
int strap_read(const char *file_path, StrapBands *s);
// header 可为 NULL；prec 是 %g 的有效位数（<=0 取 17，可无损往返）
int strap_write(const char *file_path, const StrapBands *s, int layout, int prec,
                const char *header);

// --- 处理 ---
// 每个 k 点内部跨带按能量升序重排（原地）
int strap_sort_by_energy(StrapBands *s);
// 找"各带分得最开"的 k 点当锚点：相邻能级最小间距最大的那一列
int strap_find_anchor(const StrapBands *s, size_t *anchor);
// 最小曲率拼接：自锚点向两侧逐点走，用前两步线性外推预测位置，
// 再在全部排列里选让 Σ(位置-预测)² 最小的分配（匈牙利算法，O(nb³)）。
// 要求 0 < anchor < nk-1。
int strap_reconnect(StrapBands *s, size_t anchor);
// 总粗糙度 Σ(二阶差分)²，用来判断拼接前后哪一份更像物理能带
double strap_roughness(const StrapBands *s);

#endif
