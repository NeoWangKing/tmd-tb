#!/bin/sh
# 能带计算 + gnuplot 作图
#
#   阶段一（PATH_MODE=0，路径 M-Γ-K-M'-K'，横坐标 2π/a）：
#           与 VASP 的 MoS2-pbe*.txt 对比（plot_nn / plot_tnn / plot_soc / plot_soc_spin）
#   阶段二（PATH_MODE=1，路径 Γ-K-M-Γ，横坐标 Å⁻¹）：
#           与 wannier90 的 GW 能带逐点对比（plot_gw，k 点与 GW 的 285 个点重合）
#
# 不依赖 LAPACK/BLAS（3x3 厄米本征问题由 main.c 内置的 Jacobi 求解器完成），
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

# --- 阶段一：与 VASP PBE 能带对比（PATH_MODE 默认为 0）---
build_run nn
build_run tnn      -DNNN_MODEL=1
build_run soc      -DNNN_MODEL=1 -DSOC_MODEL=1
build_run soc_spin -DNNN_MODEL=1 -DSOC_MODEL=1 -DSOC_SPIN=1

# --- 阶段二：与 wannier90 的 GW 能带对比 ---
# gw_nn 没有对应的作图脚本，只出 data/band_gw_nn.dat（被 plot_gw.gp 读取）
build_run gw_nn    -DPATH_MODE=1
build_run gw       -DNNN_MODEL=1 -DPATH_MODE=1
