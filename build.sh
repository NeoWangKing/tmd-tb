#!/bin/sh
# 四个模型 + 两个 GW 路径模型：编译、算能带、出图。
# 不依赖 LAPACK/BLAS，Linux / macOS 通用；CC 可用环境变量覆盖。
#
# 作图统一放在最后、按依赖顺序做，不在算完某条能带后"顺带"画 ——
# 有的作图脚本要读后面几步才生成的表（plot_gwfit.gp 要读 opt_main 写的
# data/gw_fit_bands.dat），顺带画的时候那张表还不存在，gnuplot 的 stats 会报
# "undefined variable: STATS_max" 并只留下一个 0 字节的 PNG。

set -xe

# 允许从任意目录调用这个脚本，所有相对路径都按脚本所在目录算
cd "$(dirname "$0")"

CC=${CC:-}
if [ -z "$CC" ]; then
    if command -v clang >/dev/null 2>&1; then CC=clang; else CC=cc; fi
fi
CFLAGS="-Wall -Wextra -std=c99 -O2 -I. -Isrc"
LIBSRC="src/mat.c src/tb.c src/herm3.c src/kpath.c src/bandtrack.c src/gwdata.c"
SRC="main.c $LIBSRC"

mkdir -p bin data img

# ---- 作图环境 ----
# pngcairo 需要 libcairo，有的机器上没编进去，那就退回老的 png 终端（能出图，
# 字体糙一点）；连 gnuplot 都没有就只警告不中断：能带数据照样生成，只是没有图。
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
# 阶段一：对 VASP PBE 能带（PATH_MODE=0，默认）
build_run nn
build_run tnn      -DNNN_MODEL=1
build_run soc      -DNNN_MODEL=1 -DSOC_MODEL=1
build_run soc_spin -DNNN_MODEL=1 -DSOC_MODEL=1 -DSOC_SPIN=1

# 阶段二：对 wannier90 的 GW 能带（PATH_MODE=1）
# gw_nn 只出数据（它没有自己的作图脚本），被 plot_gw.gp 当参考读取
build_run gw_nn    -DPATH_MODE=1
build_run gw       -DNNN_MODEL=1 -DPATH_MODE=1

# 任务①：用落地到 params.h 的 GW 参数算能带（GW 路径）
build_run gwfit    -DNNN_MODEL=1 -DPARAM_SET=2 -DPATH_MODE=1

# 拟合诊断（Step A）：GW 参考能带 vs 当前 TB 的残差；读上面生成的 data/band_gw.dat
"$CC" $CFLAGS -o bin/fit_main $LIBSRC fit_main.c -lm
./bin/fit_main > data/fit_report.txt

# 全路径拟合（Step B）：直接读 GW 参考能带，调 tb_set_params 试参数
"$CC" $CFLAGS -o bin/opt_main $LIBSRC opt_main.c -DNNN_MODEL=1 -lm
./bin/opt_main > data/opt_report.txt

# ---- 作图 ----
# 阶段一
plot gnuplot/plot_nn.gp       img/TB-NN.png
plot gnuplot/plot_tnn.gp      img/TB-NN-NNN-TNN.png
plot gnuplot/plot_soc.gp      img/TB-NN-NNN-TNN-SOC.png
plot gnuplot/plot_soc_spin.gp img/TB-NN-NNN-TNN-SOC-SPIN.png
# 阶段二
plot gnuplot/plot_gw.gp        img/GW-vs-TB.png
plot gnuplot/plot_fit.gp       img/GW-residual.png
plot gnuplot/plot_opt.gp       img/GW-fit.png
plot gnuplot/plot_reconnect.gp img/GW-reconnect.png
plot gnuplot/plot_gwfit.gp     img/GW-params-landed.png

if [ -n "$SKIPPED" ]; then
    echo "[WARN] 没装 gnuplot，这些图没生成：$SKIPPED" >&2
fi
