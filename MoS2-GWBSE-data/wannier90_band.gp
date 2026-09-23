# ---------------------------------------------------------------------------
# 输出模式：0 = 弹出 Qt 窗口；1 = 写成 PNG 文件
# 可以用命令行覆盖默认值，不必改文件：gnuplot -e "MODE=1" wannier90_band.gp
# ---------------------------------------------------------------------------
if (!exists("MODE")) MODE = 1        # 建议默认 1：本机 Qt 不可用时至少留下文件

if (MODE == 1) {
  set terminal pngcairo font "Arial,12" size 800,600 enhanced
  set output "wannier90_band.png"
  print "MODE = 1 -> 输出 wannier90_band.png"
} else {
  set terminal qt persist font "Arial,12" size 800,600 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}


set style data dots
set nokey

GAMMA1 = 0
K = 1.31633
M = 1.97450
GAMMA2 = 3.11448

set xlabel "k-path"
set ylabel "Energy (eV)"
# set xrange [GAMMA1: GAMMA2]
# set yrange [-15.81492 : 13.62639]
# set yrange [-3:5]

set arrow from K, graph 0 to K, graph 1 nohead
set arrow from M, graph 0 to M, graph 1 nohead

set xtics ("GAMMA" GAMMA1,"K" K,"M" M,"GAMMA" GAMMA2)

plot "wannier90_band.dat" using 1:2 title "wannier_band"
