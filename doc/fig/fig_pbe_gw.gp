# 文档 §3 用图：PBE 与 GW 能带
#   (a) PBE：MoS2-pbe.txt，阶段一那条路径 M-Γ-K-M'-K'（x 以 2π/l 为单位）
#   (b) GW ：MoS2-GWBSE-data/wannier90_band.dat，路径 Γ-K-M-Γ（x 以 Å⁻¹ 为单位）
# 两份数据不在同一条路径上，没法叠在一个坐标里，所以分成两栏；
# 能量都按各自的价带顶归零，带隙可以直接比。
if (!exists("MODE")) MODE = 1
OUTFILE = "doc/fig/fig_pbe_gw.png"
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,19" size 1180,520 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,19" size 1180,520 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口"
}

NV_PBE = 9
NV_GW  = 9
stats "MoS2-pbe.txt" using (column(1+NV_PBE)) nooutput; pbe_vbm = STATS_max
stats "MoS2-GWBSE-data/wannier90_band.dat" every ::0:(NV_GW-1):(284):(NV_GW-1) using 2 nooutput
gw_vbm = STATS_max

PBE = "#1f77b4"
GW_C = "#e8890c"

set multiplot layout 1,2 spacing 0.95,1.0
set yrange [-1.5:3.4]
set ylabel "Energy (eV, VBM 归零)"
set grid

# (a) PBE
set title "(a) PBE" font "Arial,18"
set xrange [0:0.598960]
set xtics ("M" 0, "Γ" 0.180988, "K" 0.389974, "M'" 0.494467, "K'" 0.598960)
set arrow from 0.180988, graph 0 to 0.180988, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from 0.389974, graph 0 to 0.389974, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from 0.494467, graph 0 to 0.494467, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from 0.402, 0 to 0.402, 1.675 nohead lc rgb "#d62728" dt 2 lw 2 front
set label 1 "E_g = 1.675 eV" at 0.415, 0.75 font "Arial,15"
set key top right font "Arial,16"
plot for [c=2:17] "MoS2-pbe.txt" using 1:(column(c)-pbe_vbm) w l lc rgb PBE lw 1.4 \
       title (c == 2 ? "PBE" : "")

# (b) GW
unset label
unset arrow
set ylabel ""
set title "(b) GW" font "Arial,18"
set xrange [0:3.11448]
set xtics ("Γ" 0, "K" 1.31633, "M" 1.97450, "Γ" 3.11448)
set arrow from 1.31633, graph 0 to 1.31633, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from 1.97450, graph 0 to 1.97450, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from 1.38, 0 to 1.38, 2.642 nohead lc rgb "#d62728" dt 2 lw 2 front
set label 2 "E_g = 2.642 eV" at 1.40, 1.20 font "Arial,15"
set key top right font "Arial,16"
plot "MoS2-GWBSE-data/wannier90_band.dat" using 1:($2-gw_vbm) w l lc rgb GW_C lw 1.4 title "GW"
unset multiplot
