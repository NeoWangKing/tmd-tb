# 文档 §4.2 用图：按排序能带拟合 vs 按物理能带拟合
#   data/fig_stitch.dat  = 并排(gw_fit_bands.dat, sorted_gw_fit_bands.dat)
#      col1-10 正确口径：k, GW 物理参考(3), TB 初值(3), TB 拟合(3)
#      col11-20 错误口径：k, GW 排序参考(3), TB 初值(3), TB 拟合(3)
#   灰色虚线来自 Wannier90 原始文件，即所有按能量排序的能带
if (!exists("MODE")) MODE = 1
OUTFILE = "doc/fig/fig_stitch.png"
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,10" size 880,760 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,10" size 880,760 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口"
}

KX = 1.31633
MX = 1.97450
GX = 3.11448
CB1 = "#1f77b4"
CB2 = "#2ca02c"
WB1 = "#d62728"
WB2 = "#ff9f40"
DATA = "data/fig_stitch.dat"

set multiplot layout 2,1 title "按排序能带拟合与按物理能带拟合" font "Arial,11"
set grid
set xrange [0:GX]
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.8 front

# 上：能带。灰虚线 = 原始按能量排序的全部能带；粗线 = 两种口径拟合出的导带
set ylabel "Energy (eV)"
set yrange [-3.0:3.5]
set key top right
plot "MoS2-GWBSE-data/wannier90_band.dat" using 1:2 w l dt 2 lc rgb "#a8a8a8" lw 1.0 \
       title "排序能带（Wannier90 原始，共 17 条）", \
     DATA using 1:9  w l lc rgb CB1 lw 2.2 title "正确口径拟合", \
     DATA using 1:10 w l lc rgb CB1 lw 2.2 notitle, \
     DATA using 1:8  w l lc rgb CB1 lw 2.2 notitle, \
     DATA using 1:19 w l lc rgb WB1 lw 2.2 title "错误口径拟合", \
     DATA using 1:20 w l lc rgb WB1 lw 2.2 notitle, \
     DATA using 1:18 w l lc rgb WB1 lw 2.2 notitle

# 下：两种口径对"物理参考"的残差
set ylabel "E_GW − E_TB (eV)"
set yrange [-0.45:0.45]
set key top right
plot DATA using 1:($3-$19) w l lc rgb WB1 lw 1.8 title "错误口径 CB1：RMS 0.129", \
     DATA using 1:($4-$20) w l lc rgb WB2 lw 1.8 title "错误口径 CB2：RMS 0.173", \
     DATA using 1:($3-$9)  w l lc rgb CB1 lw 2.4 title "正确口径 CB1：RMS 0.038", \
     DATA using 1:($4-$10) w l lc rgb CB2 lw 2.4 title "正确口径 CB2：RMS 0.024"
unset multiplot
