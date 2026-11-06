# 文档 §4.2 用图：两种口径各占一栏，每栏把"拟合目标"和"拟合结果"用不同线型分开
#   (a) 正确口径：拼接后的物理能带（目标，实线） vs 按它拟合的结果（虚线）
#   (b) 错误口径：按能量排序的能带（目标，实线） vs 按它拟合的结果（虚线）
# 数据：data/fig_stitch.dat = 并排(gw_fit_bands.dat, sorted_gw_fit_bands.dat)
#   col1-10  正确口径：k, 物理参考(3), TB 初值(3), TB 拟合(3)
#   col11-20 错误口径：k, 排序参考(3), TB 初值(3), TB 拟合(3)
if (!exists("MODE")) MODE = 1
OUTFILE = "doc/fig/fig_stitch.png"
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,16" size 880,720 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,16" size 880,720 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口"
}

KX = 1.31633
MX = 1.97450
GX = 3.11448
DATA = "data/fig_stitch.dat"
TGT_A = "#e8890c"
FIT_A = "#1f77b4"
TGT_B = "#505050"
FIT_B = "#d62728"

set multiplot layout 2,1 spacing 1.15,0.95
set grid
set xrange [0:GX]
set yrange [0.9:3.6]
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set ylabel "Energy (eV)"
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.8 front
set key top right font "Arial,13"
set dashtype 2 (5,3)

# (a) 正确口径：目标 = 拼接后的物理能带
set title "(a) 拼接后的物理能带 → 拟合　（残差 CB1 0.038 / CB2 0.024 eV）" font "Arial,14"
plot DATA using 1:3 w l lc rgb TGT_A lw 3.4 title "拟合目标：物理带 1、3（拼接后）", \
     DATA using 1:4 w l lc rgb TGT_A lw 3.4 notitle, \
     DATA using 1:9  w l lc rgb FIT_A lw 2.4 dt 2 title "拟合结果", \
     DATA using 1:10 w l lc rgb FIT_A lw 2.4 dt 2 notitle

# (b) 错误口径：目标 = 按能量排序的能带
set title "(b) 按能量排序的能带 → 拟合　（残差 CB1 0.129 / CB2 0.173 eV）" font "Arial,14"
plot "MoS2-GWBSE-data/wannier90_band.dat" using 1:2 w l dt 3 lc rgb "#dcdcdc" lw 1.0 \
       title "其它排序能带", \
     DATA using 1:13 w l lc rgb TGT_B lw 3.0 title "拟合目标：排序块 10、11", \
     DATA using 1:14 w l lc rgb TGT_B lw 3.0 notitle, \
     DATA using 1:19 w l lc rgb FIT_B lw 2.4 dt 2 title "拟合结果", \
     DATA using 1:20 w l lc rgb FIT_B lw 2.4 dt 2 notitle
unset multiplot
