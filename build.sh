#!/bin/sh
# 阶段一：四个 TB 模型的能带计算 + gnuplot 作图
#
# 不再依赖 LAPACK/BLAS（3x3 厄米本征问题由 main.c 内置的 Jacobi 求解器完成），
# 也没有 macOS/Homebrew 的硬编码路径，Linux / macOS 都可以直接跑。
#
# 用法：./build.sh           （CC 可以用环境变量覆盖，如 CC=gcc ./build.sh）

set -xe

CC=${CC:-}
if [ -z "$CC" ]; then
    if command -v clang >/dev/null 2>&1; then CC=clang; else CC=cc; fi
fi
CFLAGS="-Wall -Wextra -std=c99 -O2"

mkdir -p bin data img

# $1 = 模型名：对应 bin/main_$1、data/band_$1.dat、gnuplot/plot_$1.gp
# 其余参数原样传给编译器（用来切换 -DNNN_MODEL / -DSOC_MODEL / -DSOC_SPIN）
build_run() {
    name=$1
    shift
    "$CC" $CFLAGS -o "bin/main_$name" main.c "$@" -lm
    "./bin/main_$name" > "data/band_$name.dat"
    gnuplot "gnuplot/plot_$name.gp"
}

build_run nn
build_run tnn      -DNNN_MODEL=1
build_run soc      -DNNN_MODEL=1 -DSOC_MODEL=1
build_run soc_spin -DNNN_MODEL=1 -DSOC_MODEL=1 -DSOC_SPIN=1
