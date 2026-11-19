#!/bin/sh
# 四个模型 + 两个 GW 路径模型：编译、算能带、出图。
#
# 作图统一放在最后、按依赖顺序做，不在算完某条能带后"顺带"画 ——
# 有的作图脚本要读后面几步才生成的表（plot_gwfit.gp 要读 opt_main 写的
# data/gw_fit_bands.dat），顺带画的时候那张表还不存在，gnuplot 的 stats 会报
# "undefined variable: STATS_max" 并只留下一个 0 字节的 PNG。

set -xe

# 允许从任意目录调用这个脚本，所有相对路径都按脚本所在目录算
cd "$(dirname "$0")"

CFLAGS="-Wall -Wextra -std=c99 -O2 -I. -Iinclude"
LIBSRC="include/mat.c include/tb.c include/herm3.c include/kpath.c include/gwdata.c include/strap.c"
SRC="src/main.c $LIBSRC"

mkdir -p bin data img

# 系统兼容层
CC=${CC:-}
if [ -z "$CC" ]; then
    if command -v clang >/dev/null 2>&1; then CC=clang; else CC=cc; fi
fi

GTERM=pngcairo
HAVE_GP=1
if command -v gnuplot >/dev/null 2>&1; then
    if ! gnuplot -e 'set print "-"; print GPVAL_TERMINALS' 2>/dev/null | grep -q pngcairo; then
        GTERM=png
    fi
else
    HAVE_GP=0
fi
SKIPPED=

# $1 = gnuplot 脚本；$2 = 该脚本应该写出的 PNG（用来检查图真的画出来了）
plot() {
    if [ "$HAVE_GP" = 0 ]; then
        echo "[WARN] 没有 gnuplot，跳过 $1（本应生成 $2）" >&2
        SKIPPED="$SKIPPED $2"
        return 0
    fi
    # 写成 if ! ... 是为了不让 set -e 抢先退出，好把是哪个脚本出的问题打出来
    if ! gnuplot -e "GTERM='$GTERM'" "$1"; then
        echo "[ERROR] $1 作图失败（多半是它要读的 data/ 表还没生成，或脚本里写错了文件名）" >&2
        exit 1
    fi
    # 0 字节的 PNG 基本都是脚本中途报错留下的，直接报出来
    if [ ! -s "$2" ]; then
        echo "[ERROR] $1 没画出 $2（文件不存在或是空的）" >&2
        exit 1
    fi
}

# $1 = 名字：对应 bin/main_$1、data/band_$1.dat
# 其余参数原样传给编译器（用来切换 -DNNN_MODEL / -DSOC_MODEL / -DSOC_SPIN / -DPATH_MODE）
build_run() {
    name=$1
    shift
    "$CC" $CFLAGS -o "bin/main_$name" $SRC "$@" -lm
    "./bin/main_$name" > "data/band_$name.dat"
}

# ---- 算能带 ----
# 对 VASP PBE 能带（PATH_MODE=0）
# M → Γ → K → M' → K'
build_run nn          -DNNN_MODEL=0 -DSOC_MODEL=0 -DSOC_SPIN=0
build_run tnn         -DNNN_MODEL=1 -DSOC_MODEL=0 -DSOC_SPIN=0
build_run soc         -DNNN_MODEL=1 -DSOC_MODEL=1 -DSOC_SPIN=0
build_run soc_spin    -DNNN_MODEL=1 -DSOC_MODEL=1 -DSOC_SPIN=1

# 对 wannier90 的 GW 能带（PATH_MODE=1）
# Γ → K → M → Γ
# TODO: gw_nn 只出数据（它没有自己的作图脚本），被 plot_gw.gp 当参考读取
build_run nn_gw_path  -DNNN_MODEL=0 -DPATH_MODE=1
build_run tnn_gw_path -DNNN_MODEL=1 -DPATH_MODE=1

# 用落地到 params.h 的 GW 参数算能带（GW 路径）
build_run tnn_gw_fit  -DNNN_MODEL=1 -DPARAM_SET=2 -DPATH_MODE=1

# 拟合诊断（Step A）：GW 参考能带 vs 当前 TB 的残差；读上面生成的 data/band_tnn_gw_path.dat
"$CC" $CFLAGS -o bin/main_residual $LIBSRC src/main_residual.c -lm
./bin/main_residual > data/report_residual.txt

# 全路径拟合（Step B）：直接读 GW 参考能带，调 tb_set_params 试参数
"$CC" $CFLAGS -o bin/main_fitting $LIBSRC src/main_fitting.c -DNNN_MODEL=1 -lm
./bin/main_fitting > data/report_fitting.txt

# 对照组：不拼接、直接按排序块配对（CB_PICK 0,1），同样的拟合再跑一遍。
# 输出前缀换成 data/sorted_，不会覆盖上面的主结果。
"$CC" $CFLAGS -o bin/main_fitting_unconnected $LIBSRC src/main_fitting.c -DNNN_MODEL=1 \
    -DRECONNECT=0 -DCB_PICK_A=0 -DCB_PICK_B=1 -DOUT_PREFIX='"data/unconnected_"' -lm
./bin/main_fitting_unconnected > data/report_fitting_unconnected.txt

# 式 (6) 的核矩阵：K 点附近 N_k×N_k 网格上组装，参数用 GW 拟合那一套
# 网格与自旋两块都在程序里定，不需要额外的宏
"$CC" $CFLAGS -o bin/main_bse_kernel src/main_bse_kernel.c $LIBSRC -lm
./bin/main_bse_kernel 16 > data/report_bse_kernel.txt

# 文档插图的数据：两个拟合表按行并排（gnuplot 只能在一个文件里算列间差）
awk -f doc/fig/pair.awk data/gw_bands.dat data/unconnected_gw_bands.dat > data/fig_stitch.dat

# ---- 作图 ----
# 阶段一
mkdir -p img/TB
plot gnuplot/TB/plot_nn.gp       img/TB/NN.png
plot gnuplot/TB/plot_tnn.gp      img/TB/NN-NNN-TNN.png
plot gnuplot/TB/plot_soc.gp      img/TB/NN-NNN-TNN-SOC.png
plot gnuplot/TB/plot_soc_spin.gp img/TB/NN-NNN-TNN-SOC-SPIN.png
# 阶段二
mkdir -p img/GW
plot gnuplot/GW/plot_gw_vs_pbe.gp   img/GW/gw_vs_pbe.png
plot gnuplot/GW/plot_residual.gp   img/GW/residual.png
plot gnuplot/GW/plot_fit_result.gp img/GW/fit_result.png
plot gnuplot/GW/plot_reconnect.gp  img/GW/reconnect.png
plot gnuplot/GW/plot_result.gp     img/GW/result.png
# 文档插图（doc/theory.tex）
plot doc/fig/fig_stitch.gp  doc/fig/fig_stitch.png
plot doc/fig/fig_fit.gp     doc/fig/fig_fit.png
plot doc/fig/fig_keldysh.gp doc/fig/fig_keldysh.png

if [ -n "$SKIPPED" ]; then
    echo "[WARN] 没装 gnuplot，这些图没生成：$SKIPPED" >&2
fi


