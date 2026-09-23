# ---------------------------------------------------------------------------
# 输出模式：0 = 弹出 Qt 窗口；1 = 写成 PNG 文件
# ---------------------------------------------------------------------------
if (!exists("MODE")) MODE = 1
OUTFILE = "img/GW-vs-TB.png"
# pngcairo 需要 libcairo，缺了可以用 -e "GTERM='png'" 换成老终端（见 build.sh）
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,12" size 900,600 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,12" size 900,600 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}

# ---------------------------------------------------------------------------
# GW (wannier90) vs TB：同一路径 Γ-K-M-Γ、同一横坐标单位 Å⁻¹、两侧各自 VBM 归零
#
#   数据来源：
#     MoS2-GWBSE-data/wannier90_band.dat
#     data/band_gw.dat
#     data/band_gw_nn.dat
#   后两个由 ./build.sh 的 gw / gw_nn 目标生成。
# ---------------------------------------------------------------------------
NV_GW = 9
NK    = 285
NV_TB = 1

# GW 文件按带分块（每块 NK 行），取第 NV_GW 条带要用 every 的起止块号（从 0 起）
stats "MoS2-GWBSE-data/wannier90_band.dat" every ::0:(NV_GW-1):(NK-1):(NV_GW-1) using 2 nooutput
gw_vbm = STATS_max
stats "data/band_gw.dat"    using (column(1+NV_TB)) nooutput
tb_vbm    = STATS_max
stats "data/band_gw_nn.dat" using (column(1+NV_TB)) nooutput
tb_vbm_nn = STATS_max

# 高对称点
GX1 = 0.0
KX  = 1.31633
MX  = 1.97450
GX2 = 3.11448

set xlabel "k-path (Å^{-1})"
set ylabel "Energy (eV, VBM 归零)"
set xrange [0:GX2]
set yrange [-2.5:5]
set grid
set xtics ("Γ" GX1, "K" KX, "M" MX, "Γ" GX2)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.5 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.5 front
set key top right
set label 1 sprintf("zero: GW VBM=%.4f / TB VBM=%.4f / TB-NN VBM=%.4f eV", gw_vbm, tb_vbm, tb_vbm_nn) \
      at graph 0.012, 0.955 font "Arial,9"

plot "MoS2-GWBSE-data/wannier90_band.dat" using 1:($2 - gw_vbm) w l lc "orange" lw 3 title "GW (wannier90)", \
     "data/band_gw_nn.dat" using 1:($2 - tb_vbm_nn) w l lc "green" lw 2 dt 3 title "TB-NN", \
     "data/band_gw_nn.dat" using 1:($3 - tb_vbm_nn) w l lc "green" lw 2 dt 3 notitle, \
     "data/band_gw_nn.dat" using 1:($4 - tb_vbm_nn) w l lc "green" lw 2 dt 3 notitle, \
     "data/band_gw.dat" using 1:($2 - tb_vbm) w l lc "blue" lw 2 dt 2 title "TB-NN-NNN-TNN", \
     "data/band_gw.dat" using 1:($3 - tb_vbm) w l lc "blue" lw 2 dt 2 notitle, \
     "data/band_gw.dat" using 1:($4 - tb_vbm) w l lc "blue" lw 2 dt 2 notitle
