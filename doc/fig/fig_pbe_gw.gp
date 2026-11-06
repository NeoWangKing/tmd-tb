# 文档 §3 用图：PBE 与 GW 能带在同一段路径（Γ→K）上的叠加
#
# 为什么只画 Γ→K：仓库里的 PBE 数据（MoS2-pbe.txt）在阶段一那条
# M-Γ-K-M'-K' 上，GW 数据（wannier90_band.dat）在 Γ-K-M-Γ 上，
# 两条路径唯一重合的连续段就是 Γ→K。横坐标用这一段的归一化位置 0→1，
# 两份数据点数和单位都不同，归一化以后可以直接叠。
#
#   MoS2-pbe.txt  : 161 行 = 4 段 × 40 步 + 1，第 2 段（行 40..80）是 Γ→K
#   wannier90_band.dat : 每块 285 点，前 121 点（0..120）是 Γ→K
if (!exists("MODE")) MODE = 1
OUTFILE = "doc/fig/fig_pbe_gw.png"
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,17" size 900,620 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,17" size 900,620 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口"
}

NV = 9
stats "MoS2-pbe.txt" using (column(1+NV)) nooutput
pbe_vbm = STATS_max
stats "MoS2-GWBSE-data/wannier90_band.dat" every ::0:(NV-1):284:(NV-1) using 2 nooutput
gw_vbm = STATS_max

PBE   = "#1f77b4"
GW_C  = "#e8890c"
X0    = 0.180988
X1    = 0.389974
KS    = 1.31633

set xlabel "Γ → K（归一化路径）"
set ylabel "Energy (eV, VBM 归零)"
set xrange [-0.07:1.05]
set yrange [-1.4:3.6]
set xtics ("Γ" 0, "K" 1)
set grid
set key top left font "Arial,15"
set arrow from 0.93, 0 to 0.93, 1.675 nohead lc rgb "#1f77b4" dt 2 lw 2.5 front
set arrow from 0.99, 0 to 0.99, 2.642 nohead lc rgb "#e8890c" dt 2 lw 2.5 front
set label 1 "E_g：PBE 1.675 eV，GW 2.642 eV" at 0.12, 0.55 font "Arial,14"

plot for [c=2:17] "MoS2-pbe.txt" every ::40::80 \
       using (($1-X0)/(X1-X0)):(column(c)-pbe_vbm) w l lc rgb PBE lw 2 \
       title (c == 2 ? "PBE" : ""), \
     "MoS2-GWBSE-data/wannier90_band.dat" every ::0:0:120:16 \
       using ($1/KS):($2-gw_vbm) w l lc rgb GW_C lw 2 title "GW"
