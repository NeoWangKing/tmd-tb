# ---------------------------------------------------------------------------
# GW (wannier90) vs TB：同一路径 Γ-K-M-Γ、同一横坐标单位 Å⁻¹、两侧各自 VBM 归零
#
#   数据来源：
#     MoS2-GWBSE-data/wannier90_band.dat  —— GW 能带，17 条（无 SOC），按"带"分块存放
#     data/band_gw.dat                    —— TB-NN+NNN+TNN（PATH_MODE=1）
#     data/band_gw_nn.dat                 —— TB-仅NN（PATH_MODE=1）
#   后两个由 ./build.sh 的 gw / gw_nn 目标生成。
# ---------------------------------------------------------------------------
# 输出模式：0 = 弹出 Qt 窗口；1 = 写成 PNG 文件
# 可以用命令行覆盖默认值，不必改文件：gnuplot -e "MODE=0" plot_gw.gp
# 注意：默认值必须是 1 —— build.sh 是裸调 gnuplot（不带 -e）。
# ---------------------------------------------------------------------------
if (!exists("MODE")) MODE = 1
OUTFILE = "img/GW-vs-TB.png"

if (MODE == 1) {
  set terminal pngcairo font "Arial,12" size 900,600 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,12" size 900,600 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}

# ---------------------------------------------------------------------------
# 零点对齐：两侧各自取"价带顶那条带"的极大值
#   GW 侧：MoS2 单层 18 个价电子 -> 9 条价带。注意文件是按【带】分块存的
#          （每块 NK 行、块间以空行分隔），所以"第 NV_GW 条带"要用 every 按块取，
#          块号从 0 开始，且四个字段必须写满：
#              every ::<起始点>:<起始块>:<结束点>:<结束块>
#          （踩过的两个坑：① 不能用 column() 选带的，因为"带"不是列；
#                            ② every ::A::B 这种写法在这里不成立。）
#   TB 侧：无 SOC，1 条价带 -> 第 2 列
# ---------------------------------------------------------------------------
NV_GW = 9
NK    = 285
NV_TB = 1

stats "MoS2-GWBSE-data/wannier90_band.dat" every ::0:(NV_GW-1):(NK-1):(NV_GW-1) using 2 nooutput
gw_vbm = STATS_max
stats "data/band_gw.dat"    using (column(1+NV_TB)) nooutput
tb_vbm    = STATS_max
stats "data/band_gw_nn.dat" using (column(1+NV_TB)) nooutput
tb_vbm_nn = STATS_max

# 高对称点（Å⁻¹），取自 wannier90_band.labelinfo.dat
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
