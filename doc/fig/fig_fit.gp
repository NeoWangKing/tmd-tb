# 文档 §5.3 用图：GW 参考能带与拟合后的三能带模型
# 数据：data/gw_fit_bands.dat（col2-4 GW 参考，col5-7 TB 初值，col8-10 TB 拟合）
if (!exists("MODE")) MODE = 1
OUTFILE = "doc/fig/fig_fit.png"
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,18" size 880,560 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,18" size 880,560 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口"
}

stats "data/gw_fit_bands.dat" using 2 nooutput; gw_vbm   = STATS_max
stats "data/gw_fit_bands.dat" using 5 nooutput; init_vbm = STATS_max
stats "data/gw_fit_bands.dat" using 8 nooutput; fit_vbm  = STATS_max

KX = 1.31633
MX = 1.97450
GX = 3.11448

set xlabel "k-path (Å^{-1})"
set ylabel "Energy (eV, VBM 归零)"
set xrange [0:GX]
set yrange [-1.9:5.7]
set grid
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.8 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.8 front
set key top right font "Arial,12"
set dashtype 2 (5,3)

plot "data/gw_fit_bands.dat" using 1:($2-gw_vbm) w l lc "orange" lw 3   title "GW 参考", \
     "data/gw_fit_bands.dat" using 1:($3-gw_vbm) w l lc "orange" lw 3   notitle, \
     "data/gw_fit_bands.dat" using 1:($4-gw_vbm) w l lc "orange" lw 3   notitle, \
     "data/gw_fit_bands.dat" using 1:($8-fit_vbm) w l lc "#1f77b4" lw 2.2 dt 2 title "拟合后（RMS 33 meV）", \
     "data/gw_fit_bands.dat" using 1:($9-fit_vbm) w l lc "#1f77b4" lw 2.2 dt 2 notitle, \
     "data/gw_fit_bands.dat" using 1:($10-fit_vbm) w l lc "#1f77b4" lw 2.2 dt 2 notitle, \
     "data/gw_fit_bands.dat" using 1:($5-init_vbm) w l lc rgb "#606060" lw 1.8 dt 3 title "上一阶段的模型（初值）", \
     "data/gw_fit_bands.dat" using 1:($6-init_vbm) w l lc rgb "#606060" lw 1.8 dt 3 notitle, \
     "data/gw_fit_bands.dat" using 1:($7-init_vbm) w l lc rgb "#606060" lw 1.8 dt 3 notitle
