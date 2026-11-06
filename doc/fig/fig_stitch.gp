# 文档 §4.2 用图：按能量排序的能级 vs 拼接后的物理能带，以及用错参考时的残差
# 数据：data/fig_stitch.dat（build.sh 里由 doc/fig/pair.awk 把下面两个文件并排拼成）
#   col1-16  data/gw_reconnect.dat：k, 排序块 10..14, 物理带 1..5, 来源块号
#   col17-26 data/gw_fit_bands.dat：k, GW 参考(3), TB 初值(3), TB 拟合(3)
if (!exists("MODE")) MODE = 1
OUTFILE = "doc/fig/fig_stitch.png"
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,10" size 880,740 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,10" size 880,740 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口"
}

KX = 1.31633
MX = 1.97450
GX = 3.11448
C1 = "#1f77b4"
C3 = "#2ca02c"
DATA = "data/fig_stitch.dat"

set multiplot layout 2,1 title "排序本征值与拼接后的物理能带" font "Arial,11"
set grid
set xrange [0:GX]
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.8 front
set key top right

# 上：能带
set ylabel "Energy (eV)"
set yrange [0.9:3.5]
plot DATA using 1:2 w l lc rgb "black" lw 1.5 dt 2 title "排序块 10、11（按能量）", \
     DATA using 1:3 w l lc rgb "black" lw 1.5 dt 2 notitle, \
     DATA using 1:7 w l lc rgb C1 lw 2.6 title "物理带 1（拼接后）", \
     DATA using 1:9 w l lc rgb C3 lw 2.6 title "物理带 3（拼接后）"

# 下：残差。细线 = 拿排序块当参考；粗线 = 拿拼接后的物理带当参考
set ylabel "E_GW − E_TB (eV)"
set yrange [-0.65:0.65]
set key top right
plot DATA using 1:($2-$25)  w l lc rgb "#e8a0a0" lw 1.4 title "CB1，参考用排序块 10", \
     DATA using 1:($3-$26)  w l lc rgb "#f0c070" lw 1.4 title "CB2，参考用排序块 11", \
     DATA using 1:($19-$25) w l lc rgb C1       lw 2.4 title "CB1，参考用物理带 1", \
     DATA using 1:($20-$26) w l lc rgb C3       lw 2.4 title "CB2，参考用物理带 3"
unset multiplot
