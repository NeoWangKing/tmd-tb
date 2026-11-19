# 文档幻灯片用图：DFT 层面的三能带模型 vs GW 准粒子能带，同一坐标系
#   DFT 一侧用文献参数的三能带模型（../../data/band_tnn_gw_path.dat，它拟合的是 PBE 能带）
#   GW 一侧用拼接后的物理带（../../data/gw_bands.dat 的 col2-4）
# 两者都按各自的价带顶归零，画在 Γ-K-M-Γ 同一条 k 轴上。
set terminal pngcairo font "Arial,19" size 760,540 enhanced
set output "dft_vs_gw.png"

set xlabel "k (\AA^{-1})"
set ylabel "E (eV，VBM 归零)"
set xtics ("Γ" 0, "K" 1.316334, "M" 1.974501, "Γ" 3.114479)
set xrange [0:3.114479]
set yrange [-3.4:4.9]
set grid
set arrow from 1.316334, graph 0 to 1.316334, graph 1 nohead lc rgb "black" lw 0.8
set arrow from 1.974501, graph 0 to 1.974501, graph 1 nohead lc rgb "black" lw 0.8
set key bottom right font "Arial,17"
# 在 K 点标出两条导带底的能量差
set arrow 3 from 1.346, 1.658 to 1.346, 2.6416 heads size 0.012,0.16 lw 2.5 lc rgb "#444444" front
# 能量差的数值写在幻灯片图注里，图上只留箭头避免压住曲线

stats "../../data/band_tnn_gw_path.dat" using 2 nooutput
tb_vbm = STATS_max
stats "../../data/gw_bands.dat" using 2 nooutput
gw_vbm = STATS_max

DFT = "#1f77b4"
GW  = "#e8890c"

plot "../../data/band_tnn_gw_path.dat" using 1:($2-tb_vbm) w l lw 2.2 lc rgb DFT t "DFT 层面（三能带模型）", \
     "" using 1:($3-tb_vbm) w l lw 2.2 lc rgb DFT notitle, \
     "" using 1:($4-tb_vbm) w l lw 2.2 lc rgb DFT notitle, \
     "../../data/gw_bands.dat" using 1:($2-gw_vbm) w l lw 2.2 lc rgb GW t "GW 准粒子能带", \
     "" using 1:($3-gw_vbm) w l lw 2.2 lc rgb GW notitle, \
     "" using 1:($4-gw_vbm) w l lw 2.2 lc rgb GW notitle
