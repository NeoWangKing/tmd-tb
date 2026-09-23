#!/bin/sh
# 四个模型 + 两个 GW 路径模型：编译、算能带、出图。
# 不依赖 LAPACK/BLAS，Linux / macOS 通用；CC 可用环境变量覆盖。

set -xe

CC=${CC:-}
if [ -z "$CC" ]; then
    if command -v clang >/dev/null 2>&1; then CC=clang; else CC=cc; fi
fi
CFLAGS="-Wall -Wextra -std=c99 -O2"

mkdir -p bin data img

# $1 = 名字：对应 bin/main_$1、data/band_$1.dat；若存在 gnuplot/plot_$1.gp 则顺带作图
# 其余参数原样传给编译器（用来切换 -DNNN_MODEL / -DSOC_MODEL / -DSOC_SPIN / -DPATH_MODE）
build_run() {
    name=$1
    shift
    "$CC" $CFLAGS -o "bin/main_$name" main.c "$@" -lm
    "./bin/main_$name" > "data/band_$name.dat"
    if [ -f "gnuplot/plot_$name.gp" ]; then
        gnuplot "gnuplot/plot_$name.gp"
    fi
}

# 阶段一：对 VASP PBE 能带（PATH_MODE=0，默认）
build_run nn
build_run tnn      -DNNN_MODEL=1
build_run soc      -DNNN_MODEL=1 -DSOC_MODEL=1
build_run soc_spin -DNNN_MODEL=1 -DSOC_MODEL=1 -DSOC_SPIN=1

# 阶段二：对 wannier90 的 GW 能带（PATH_MODE=1）
# gw_nn 只出数据，被 plot_gw.gp 读取
build_run gw_nn    -DPATH_MODE=1
build_run gw       -DNNN_MODEL=1 -DPATH_MODE=1
