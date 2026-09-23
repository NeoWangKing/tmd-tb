# 任务①收尾：落地到 params.h 的 GW 参数算出的 TB 能带 vs GW 参考
#   GW 参考（3 条拟合目标）: data/gw_fit_bands.dat 第 2~4 列
#   TB（文献参数）        : data/band_gw.dat
#   TB（GW 参数，PARAM_SET=2）: data/band_gwfit.dat
if (!exists("MODE")) MODE = 1
OUTFILE = "img/GW-params-landed.png"
# pngcairo 需要 libcairo，缺了可以用 -e "GTERM='png'" 换成老终端（见 build.sh）
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,11" size 950,620 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,11" size 950,620 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}

# 各自的价带顶（零点）
stats "data/gw_fit_bands.dat" using 2 nooutput; gw_vbm   = STATS_max
stats "data/gw_fit_bands.dat" using 5 nooutput; init_vbm = STATS_max
stats "data/band_gw.dat"      using 2 nooutput; lit_vbm = STATS_max
stats "data/band_gwfit.dat"   using 2 nooutput; fit_vbm = STATS_max

# "去零点后"的 RMS（都在 data/gw_fit_bands.dat 里：第 2~4 列是 GW 目标，
# 第 5~7 列是文献参数，第 8~10 列是 GW 拟合参数）
stats "data/gw_fit_bands.dat" using ((($8-fit_vbm) - ($2-gw_vbm))**2) nooutput; r_f1 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($9-fit_vbm) - ($3-gw_vbm))**2) nooutput; r_f2 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($10-fit_vbm) - ($4-gw_vbm))**2) nooutput; r_f3 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($5-init_vbm) - ($2-gw_vbm))**2) nooutput; r_l1 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($6-init_vbm) - ($3-gw_vbm))**2) nooutput; r_l2 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($7-init_vbm) - ($4-gw_vbm))**2) nooutput; r_l3 = sqrt(STATS_mean)

KX = 1.31633
MX = 1.97450
GX = 3.11448

set xlabel "k-path (Å^{-1})"
set ylabel "Energy (eV, VBM 归零)"
set xrange [0:GX]
set yrange [-1.9:5.2]
set grid
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.5 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.5 front
set key top right
set label 1 sprintf("GW 参数 RMS: VB %.3f  CB1 %.3f  CB2 %.3f eV", r_f1, r_f2, r_f3) \
      at graph 0.015, 0.95 font "Arial,9"
set label 2 sprintf("文献参数 RMS: VB %.3f  CB1 %.3f  CB2 %.3f eV", r_l1, r_l2, r_l3) \
      at graph 0.015, 0.90 font "Arial,9"
set dashtype 2 (5,3)

plot "MoS2-GWBSE-data/wannier90_band.dat" using 1:($2-gw_vbm) w l lc rgb "#f2b070" lw 1 title "GW 其他能带", \
     "data/gw_fit_bands.dat" using 1:($2-gw_vbm) w l lc "orange" lw 3 title "GW 参考（拟合目标）", \
     "data/gw_fit_bands.dat" using 1:($3-gw_vbm) w l lc "orange" lw 3 notitle, \
     "data/gw_fit_bands.dat" using 1:($4-gw_vbm) w l lc "orange" lw 3 notitle, \
     "data/band_gwfit.dat"   using 1:($2-fit_vbm) w l lc "blue" lw 2.5 dt 2 title "TB（GW 参数，已落地）", \
     "data/band_gwfit.dat"   using 1:($3-fit_vbm) w l lc "blue" lw 2.5 dt 2 notitle, \
     "data/band_gwfit.dat"   using 1:($4-fit_vbm) w l lc "blue" lw 2.5 dt 2 notitle, \
     "data/band_gw.dat"      using 1:($2-lit_vbm) w l lc rgb "#404040" lw 2.4 dt 3 title "TB（文献参数）", \
     "data/band_gw.dat"      using 1:($3-lit_vbm) w l lc rgb "#404040" lw 2.4 dt 3 notitle, \
     "data/band_gw.dat"      using 1:($4-lit_vbm) w l lc rgb "#404040" lw 2.4 dt 3 notitle
