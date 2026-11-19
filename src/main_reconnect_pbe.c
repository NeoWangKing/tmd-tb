// 把 MoS2-pbe-ncl.txt 里按能量排序的能带做一次最小曲率拼接，输出同格式的列式表。
//
// 输入是 M-Γ-K-M'-K' 路径上 160 个 k 点、32 条能带的 ncl（含 SOC）计算结果：
// 第一列是 k 路径（单位 2π），后面 32 列是能量，每个 k 点内部按能量从小到大写。
// 两条能带交叉时这个顺序会互换，画出来的曲线在交叉处变成两条互相弹开的折线。
// 拼接（include/strap.h 的 strap_reconnect）从锚点向两侧逐点走，用前两步线性
// 外推预测下一列各带的位置，再在全部排列里选让总位移平方和最小的那个分配，
// 把被排序打乱的身份接回去。
//
// 用法：
//   ./bin/main_reconnect_pbe [输入] [输出] [锚点序号]
// 三个参数都可以不给：
//   输入默认 MoS2-pbe-ncl.txt，输出默认 MoS2-pbe-ncl-reconnected.txt，
//   锚点默认交给 strap_find_anchor 自动找（相邻能级分得最开的那个 k 点）。
// 想换锚点就传第 3 个参数：0 基编号（文件里第一行数据是第 0 个），
// 只能取 1 .. nk-2 —— 两端没有前两步可以外推。

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "strap.h"

// 输出精度：输入是 6 位小数，10 位有效数字足够原样往返，又不会拖出 17 位尾数
#define OUT_SIG 10

int main(int argc, char **argv)
{
    const char *f_in  = (argc > 1) ? argv[1] : "MoS2-pbe-ncl.txt";
    const char *f_out = (argc > 2) ? argv[2] : "MoS2-pbe-ncl-reconnected.txt";

    StrapBands d;
    strap_init(&d);
    if (strap_read(f_in, &d)) return 1;

    size_t nb = strap_band_count(&d);
    size_t nk = strap_sample_count(&d);
    if (nb < 2 || nk < 3) {
        fprintf(stderr, "[ERROR] 至少要 2 条带、3 个 k 点（实际 %zu 条带、%zu 个点）\n", nb, nk);
        strap_free(&d);
        return 1;
    }

    // 输入本来就每个 k 点按能量升序；这里再排一次，把"拼接的起点一定是排序块"
    // 这件事写死，以后换一份顺序不同的数据也不会得到别的结果。
    if (strap_sort_by_energy(&d)) { strap_free(&d); return 1; }

    double rough_before = strap_roughness(&d);

    size_t anchor = 0;
    if (argc > 3) {
        char *end = NULL;
        errno = 0;
        long v = strtol(argv[3], &end, 10);
        if (errno || end == argv[3] || *end != '\0' || v < 1 || (size_t)v + 1 >= nk) {
            fprintf(stderr, "[ERROR] 锚点序号（0 基）要在 1..%zu 之间（给了 \"%s\"）\n", nk - 2, argv[3]);
            strap_free(&d);
            return 1;
        }
        anchor = (size_t)v;
    } else {
        if (strap_find_anchor(&d, &anchor)) { strap_free(&d); return 1; }
    }

    // 拼接前的值留一份，用来统计有多少个位置真的被改派了
    double *before = malloc(sizeof *before * nb * nk);
    if (!before) {
        fprintf(stderr, "[ERROR] buy more ram lol!\n");
        strap_free(&d);
        return 1;
    }
    for (size_t b = 0; b < nb; ++b)
        for (size_t i = 0; i < nk; ++i)
            before[b * nk + i] = d.bands[b].values[i];

    if (strap_reconnect(&d, anchor)) {
        free(before);
        strap_free(&d);
        return 1;
    }

    double rough_after = strap_roughness(&d);

    // strap_reconnect 是把别的带的值搬到这条带来，值本身没变，所以直接比相等即可
    size_t moved = 0;
    for (size_t b = 0; b < nb; ++b)
        for (size_t i = 0; i < nk; ++i)
            if (before[b * nk + i] != d.bands[b].values[i]) ++moved;
    free(before);

    char header[256];
    snprintf(header, sizeof header,
             "%s minimum-curvature reconnected (anchor k=%.6f); "
             "kpath(in_2pi) band-levels(nkpoints_x_nbands)",
             f_in, d.bands[0].keys[anchor]);

    if (strap_write(f_out, &d, STRAP_LAYOUT_COLUMN, OUT_SIG, header)) {
        strap_free(&d);
        return 1;
    }

    printf("最小曲率拼接：%s -> %s\n", f_in, f_out);
    printf("  规模：%zu 条带 × %zu 个 k 点\n", nb, nk);
    printf("  锚点：第 %zu 个 k 点（k = %.6f，%s）\n", anchor, d.bands[0].keys[anchor],
           (argc > 3) ? "命令行指定" : "strap_find_anchor 自动选取");
    printf("  粗糙度 Σ(二阶差分)²：%.8g -> %.8g\n", rough_before, rough_after);
    printf("  被改派的位置：%zu / %zu（%.1f%%）\n",
           moved, nb * nk, 100.0 * (double)moved / (double)(nb * nk));

    strap_free(&d);
    return 0;
}
