#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>
#include "strap.h"

#define OUT_SIG 10

char *args_shift(int *argc, char ***argv)
{
    assert(*argc > 0);
    char *result = **argv;
    (*argc) -= 1;
    (*argv) += 1;
    return result;
}

int main(int argc, char **argv)
{
    const char *program = args_shift(&argc, &argv);

    if (argc <= 0) {
        fprintf(stderr, "[Usage] %s <input_file> <output_file>\n", program);
        fprintf(stderr, "[ERROR] no input file is provided\n");
        return 1;
    }
    const char *input_file  = args_shift(&argc, &argv);

    if (argc <= 0) {
        fprintf(stderr, "[Usage] %s <input_file> <output_file>\n", program);
        fprintf(stderr, "[ERROR] no output file is provided\n");
        return 1;
    }
    const char *output_file  = args_shift(&argc, &argv);

    StrapBands strap;
    strap_init(&strap);
    if (strap_read(input_file, &strap)) return 1;

    size_t band_count = strap_band_count(&strap);
    size_t sample_count = strap_sample_count(&strap);
    if (band_count < 2 || sample_count < 3) {
        fprintf(stderr, "[ERROR] there are only %zu band(s), %zu point(s) (more than 2 bands, 3 sample pints are required)\n", band_count, sample_count);
        strap_free(&strap);
        return 1;
    }

    if (strap_sort_by_energy(&strap)) { strap_free(&strap); return 1; }

    double rough_before = strap_roughness(&strap);

    size_t anchor = 0;
    if (argc <= 0) {
        if (strap_find_anchor(&strap, &anchor)) { strap_free(&strap); return 1; }
    } else {
        char *end = NULL;
        errno = 0;
        anchor = (size_t)(strtol(args_shift(&argc, &argv), &end, 10));
        if (errno || end == argv[3] || *end != '\0' || anchor < 1 || (size_t)anchor + 1 >= sample_count) {
            fprintf(stderr, "[ERROR] the anchor located at %zu (number between 1..%zu is required)\n", anchor, sample_count - 2);
            strap_free(&strap);
            return 1;
        }
    }

    StrapBands before;
    strap_init(&before);

    for (size_t b = 0; b < band_count; ++b) {
        for (size_t i = 0; i < sample_count; ++i) {
            before.bands[b].values[i] = strap.bands[b].values[i];
        }
            before.bands[b].count = strap.bands[b].count;
    }

    if (strap_reconnect(&strap, anchor)) {
        strap_free(&before);
        strap_free(&strap);
        return 1;
    }

    double rough_after = strap_roughness(&strap);

    // strap_reconnect 是把别的带的值搬到这条带来，值本身没变，所以直接比相等即可
    size_t moved = 0;
    for (size_t b = 0; b < band_count; ++b)
        for (size_t i = 0; i < sample_count; ++i)
            if (before[b * sample_count + i] != strap.bands[b].values[i]) ++moved;
    free(before);

    char header[256];
    snprintf(header, sizeof header,
             "%s minimum-curvature reconnected (anchor k=%.6f); "
             "kpath(in_2pi) band-levels(nkpoints_x_nbands)",
             input_file, strap.bands[0].keys[anchor]);

    if (strap_write(output_file, &strap, STRAP_LAYOUT_COLUMN, OUT_SIG, header)) {
        strap_free(&strap);
        return 1;
    }

    printf("最小曲率拼接：%s -> %s\n", input_file, output_file);
    printf("  规模：%zu 条带 × %zu 个 k 点\n", band_count, sample_count);
    printf("  锚点：第 %zu 个 k 点（k = %.6f，%s）\n", anchor, strap.bands[0].keys[anchor],
           (argc > 3) ? "命令行指定" : "strap_find_anchor 自动选取");
    printf("  粗糙度 Σ(二阶差分)²：%.8g -> %.8g\n", rough_before, rough_after);
    printf("  被改派的位置：%zu / %zu（%.1f%%）\n",
           moved, band_count * sample_count, 100.0 * (double)moved / (double)(band_count * sample_count));

    strap_free(&strap);
    return 0;
}
