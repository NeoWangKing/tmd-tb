# 全路径拟合结果：GW 参考 vs 拟合后的 TB（各自 VBM 归零）
# 数据由 bin/opt_main 生成（data/gw_fit_bands.dat）
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

# 各自的价带顶（零点），以及"去掉零点后"的 RMS —— 都从数据里算，不写死
stats "data/gw_fit_bands.dat" using 2 nooutput; gw_vbm   = STATS_max
stats "data/gw_fit_bands.dat" using 4 nooutput; init_vbm = STATS_max
stats "data/gw_fit_bands.dat" using 6 nooutput; fit_vbm  = STATS_max

stats "data/gw_fit_bands.dat" using ((($6-fit_vbm) - ($2-gw_vbm))**2) nooutput; r_fit_vb = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($7-fit_vbm) - ($3-gw_vbm))**2) nooutput; r_fit_cb = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($4-init_vbm) - ($2-gw_vbm))**2) nooutput; r_ini_vb = sqrt(STATS_mean)
stats "data/gw_fit_bands.dat" using ((($5-init_vbm) - ($3-gw_vbm))**2) nooutput; r_ini_cb = sqrt(STATS_mean)

KX = 1.31633
MX = 1.97450
GX = 3.11448

set xlabel "k-path (Å^{-1})"
set ylabel "Energy (eV, VBM 归零)"
set xrange [0:GX]
set yrange [-1.6:4.6]
set grid
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.5 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.5 front
set key top right
set label 1 sprintf("RMS(去零点): 拟合 VB %.3f / CB %.3f eV   |   初值 VB %.3f / CB %.3f eV", \
                    r_fit_vb, r_fit_cb, r_ini_vb, r_ini_cb) at graph 0.015, 0.95 font "Arial,9"

set dashtype 2 (5,3)
plot "data/gw_fit_bands.dat" using 1:($2-gw_vbm)   w l lc "orange" lw 3   title "GW 参考 (VB)", \
     "data/gw_fit_bands.dat" using 1:($3-gw_vbm)   w l lc "orange" lw 3   notitle, \
     "data/gw_fit_bands.dat" using 1:($6-fit_vbm)  w l lc "blue"   lw 2.5 dt 2 title "拟合 TB (VB)", \
     "data/gw_fit_bands.dat" using 1:($7-fit_vbm)  w l lc "blue"   lw 2.5 dt 2 notitle, \
     "data/gw_fit_bands.dat" using 1:($4-init_vbm) w l lc "gray50" lw 1.5 dt 3 title "初值 TB (VB)", \
     "data/gw_fit_bands.dat" using 1:($5-init_vbm) w l lc "gray50" lw 1.5 dt 3 notitle
