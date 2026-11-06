# 文档 §6 用图：Keldysh 屏蔽势与裸库仑势的对比（双对数）
# 纯公式作图，不依赖数据文件
if (!exists("MODE")) MODE = 1
OUTFILE = "doc/fig/fig_keldysh.png"
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,10" size 760,520 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,10" size 760,520 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口"
}

r0   = 13.55
epsd = 2.5
kq(q) = 1 + r0*q
u(q)    = 90.47 / (epsd * q * kq(q))
bare(q) = 90.47 / (epsd * q)

set logscale xy
set xlabel "q (Å^{-1})"
set ylabel "|u(q)| (eV Å^2)"
set xrange [1e-3:2]
set yrange [0.4:1e5]
set grid
set key top right
set arrow from 1/r0, graph 0 to 1/r0, graph 1 nohead dt 3 lc rgb "#606060"
set label 1 "q = 1/r_0" at 1/r0*1.3, 1.2e4 font "Arial,9"
set label 2 "1/q" at 2e-3, 120 font "Arial,10" tc rgb "#808080"
set label 3 "1/(r_0 q^2)" at 1.15, 2.2 font "Arial,10" tc rgb "#1f77b4"

plot bare(x) w l lc rgb "#808080" lw 1.8 dt 2 title "裸库仑（r_0 = 0）", \
     u(x)    w l lc rgb "#1f77b4" lw 2.6 title "Keldysh（r_0 = 13.55 Å，ε_d = 2.5）"
