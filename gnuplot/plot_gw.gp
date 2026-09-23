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
# 两侧的 3 条对比能带都是"拼接后的物理能带"（交叉处不互换身份，见 src/bandtrack.h）
#
#   数据来源：
#     MoS2-GWBSE-data/wannier90_band.dat   GW 原始文件（只作背景）
#     data/gw_reconnect.dat                GW 拼接后的物理能带（由 bin/opt_main 生成）
#     data/band_gw.dat                     TB（NN+NNN+TNN，文献参数）
#     data/band_gw_nn.dat                  TB（只有 NN）
#   后两个由 ./build.sh 的 gw / gw_nn 目标生成，PATH_MODE=1 时输出顺序就是物理能带。
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

# GW 侧的 3 条对比能带用拼接后的物理带：价带是原始第 9 条（跟上下都分得开，不用拼），
# 两条导带取 data/gw_reconnect.dat 的物理带 1 与物理带 3（与拟合用的配对一致）。
# 原始文件只画成淡色细线作背景 —— 它把交叉的两条带按能量排在一起，身份是错的，
# 而 TB 那边（data/band_gw*.dat）也是物理能带顺序，两边口径必须一致。
set dashtype 2 (3,3)
plot "MoS2-GWBSE-data/wannier90_band.dat" using 1:($2 - gw_vbm) w l lc rgb "#f2b070" lw 1 title "GW 其它能带", \
     "MoS2-GWBSE-data/wannier90_band.dat" every ::0:(NV_GW-1):(NK-1):(NV_GW-1) \
          using 1:($2 - gw_vbm) w l lc "orange" lw 3 title "GW 参考（拼接后）", \
     "data/gw_reconnect.dat" using 1:($7 - gw_vbm) w l lc "orange" lw 3 notitle, \
     "data/gw_reconnect.dat" using 1:($9 - gw_vbm) w l lc "orange" lw 3 notitle, \
     "data/band_gw_nn.dat" using 1:($2 - tb_vbm_nn) w l lc "green" lw 2 dt 3 title "TB-NN", \
     "data/band_gw_nn.dat" using 1:($3 - tb_vbm_nn) w l lc "green" lw 2 dt 3 notitle, \
     "data/band_gw_nn.dat" using 1:($4 - tb_vbm_nn) w l lc "green" lw 2 dt 3 notitle, \
     "data/band_gw.dat" using 1:($2 - tb_vbm) w l lc "blue" lw 2 dt 2 title "TB-NN-NNN-TNN", \
     "data/band_gw.dat" using 1:($3 - tb_vbm) w l lc "blue" lw 2 dt 2 notitle, \
     "data/band_gw.dat" using 1:($4 - tb_vbm) w l lc "blue" lw 2 dt 2 notitle
