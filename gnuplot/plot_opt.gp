# 全路径拟合结果：GW 参考 vs 拟合后的 TB（各自 VBM 归零）
# 数据由 bin/opt_main 生成（data/gw_fit_bands.dat，列见文件头）
if (!exists("MODE")) MODE = 1
OUTFILE = "img/GW-fit.png"

if (MODE == 1) {
  set terminal pngcairo font "Arial,12" size 900,600 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,12" size 900,600 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}

# 三条带各自的价带顶（零点）与"去掉零点后"的 RMS，全部从数据里算
stats "data/gw_fit_bands.dat" using 2 nooutput; gw_vbm   = STATS_max
stats "data/gw_fit_bands.dat" using 5 nooutput; init_vbm = STATS_max
stats "data/gw_fit_bands.dat" using 8 nooutput; fit_vbm  = STATS_max

stats "data/gw_fit_bands.dat" using ((($8-fit_vbm) - ($2-gw_vbm))**2) nooutput; rf1 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($9-fit_vbm) - ($3-gw_vbm))**2) nooutput; rf2 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($10-fit_vbm) - ($4-gw_vbm))**2) nooutput; rf3 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($5-init_vbm) - ($2-gw_vbm))**2) nooutput; ri1 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($6-init_vbm) - ($3-gw_vbm))**2) nooutput; ri2 = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($7-init_vbm) - ($4-gw_vbm))**2) nooutput; ri3 = sqrt(STATS_mean)

KX = 1.31633
MX = 1.97450
GX = 3.11448

set xlabel "k-path (Å^{-1})"
set ylabel "Energy (eV, VBM 归零)"
set xrange [0:GX]
set yrange [-1.8:5.2]
set grid
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.5 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.5 front
set key top right
set label 1 sprintf("拟合 RMS: VB %.3f  CB1 %.3f  CB2 %.3f eV", rf1, rf2, rf3) \
      at graph 0.015, 0.95 font "Arial,9"
set label 2 sprintf("初值 RMS: VB %.3f  CB1 %.3f  CB2 %.3f eV", ri1, ri2, ri3) \
      at graph 0.015, 0.90 font "Arial,9"

set dashtype 2 (5,3)
# 第一行：GW 的全部能带（原始文件里 17 条按块存放，一条 plot 命令即可全部画出）
# —— 淡橙色细线，只作为背景参考（窗口外的会被裁掉）
plot "MoS2-GWBSE-data/wannier90_band.dat" using 1:($2-gw_vbm) w l lc rgb "#f2b070" lw 1 title "GW 其他能带", \
     "data/gw_fit_bands.dat" using 1:($2-gw_vbm)   w l lc "orange" lw 3   title "GW 参考（拟合用的 3 条）", \
     "data/gw_fit_bands.dat" using 1:($3-gw_vbm)   w l lc "orange" lw 3   notitle, \
     "data/gw_fit_bands.dat" using 1:($4-gw_vbm)   w l lc "orange" lw 3   notitle, \
     "data/gw_fit_bands.dat" using 1:($8-fit_vbm)  w l lc "blue"   lw 2.5 dt 2 title "拟合 TB", \
     "data/gw_fit_bands.dat" using 1:($9-fit_vbm)  w l lc "blue"   lw 2.5 dt 2 notitle, \
     "data/gw_fit_bands.dat" using 1:($10-fit_vbm) w l lc "blue"   lw 2.5 dt 2 notitle, \
     "data/gw_fit_bands.dat" using 1:($5-init_vbm) w l lc rgb "#404040" lw 2.4 dt 3 title "初值 TB", \
     "data/gw_fit_bands.dat" using 1:($6-init_vbm) w l lc rgb "#404040" lw 2.4 dt 3 notitle, \
     "data/gw_fit_bands.dat" using 1:($7-init_vbm) w l lc rgb "#404040" lw 2.4 dt 3 notitle
