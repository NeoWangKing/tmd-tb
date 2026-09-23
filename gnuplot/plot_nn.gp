# ---------------------------------------------------------------------------
# 输出模式：0 = 弹出 Qt 窗口；1 = 写成 PNG 文件
# 可以用命令行覆盖默认值，不必改文件：gnuplot -e "MODE=0" plot_nn.gp
# 注意：默认值必须是 1 —— build.sh 是裸调 gnuplot（不带 -e），改成 0 之后
#       build.sh 就不会再产出 img/ 里的图了，而且本机 Qt 终端开不了窗口时
#       它是静默失败，不会有任何报错。
# ---------------------------------------------------------------------------
if (!exists("MODE")) MODE = 1
OUTFILE = "img/TB-NN.png"

if (MODE == 1) {
  set terminal pngcairo font "Arial,12" size 800,600 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,12" size 800,600 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}

map(x) = (x <= 0.577350) ? 0 + (x-0)*(0.180988-0)/(0.577350-0) : \
         (x <= 1.244017) ? 0.180988 + (x-0.577350)*(0.389974-0.180988)/(1.244017-0.577350) : \
         (x <= 1.577350) ? 0.389974 + (x-1.244017)*(0.494467-0.389974)/(1.577350-1.244017) : \
         (x <= 1.910684) ? 0.494467 + (x-1.577350)*(0.598960-0.494467)/(1.910684-1.577350) : 0.598960

M  = map(0.0)
G  = map(0.577350)
K  = map(1.244017)
Mp = map(1.577350)
Kp = map(1.910684)

# ---------------------------------------------------------------------------
# 零点对齐：价带条数由"电子数 / 模型"决定，两侧用同一条规则取价带顶的极大值。
#   MoS2 单层有 18 个价电子 -> VASP 侧 9 条价带（含 SOC 时每带 2 个自旋分量）
#   本 TB 模型（无 SOC）：1 条价带 + 2 条导带
# 两个坑：
#   1) 必须写 column(1+NV)，写 (1+NV) 会被当成"数值 10"而不是"第 10 列"；
#   2) stats 要放在 set xrange/yrange 之前，否则数据点会被坐标范围过滤掉。
# ---------------------------------------------------------------------------
NV_VASP = 9
NV_TB   = 1

stats "data/band_nn.dat" using (column(1+NV_TB))   nooutput
tb_vbm   = STATS_max
stats "MoS2-pbe.txt"     using (column(1+NV_VASP)) nooutput
vasp_vbm = STATS_max

set xlabel "k-path"
set ylabel "Energy (eV)"
set xrange [M:Kp]
set yrange [-1:5]
set grid

set xtics ("M" 0, "Γ" G, "K" K, "M'" Mp, "K'" Kp)

# 画竖线
set arrow from M,  graph 0 to M,  graph 1 nohead lc black lw 0.5 front
set arrow from G,  graph 0 to G,  graph 1 nohead lc black lw 0.5 front
set arrow from K,  graph 0 to K,  graph 1 nohead lc black lw 0.5 front
set arrow from Mp, graph 0 to Mp, graph 1 nohead lc black lw 0.5 front
set arrow from Kp, graph 0 to Kp, graph 1 nohead lc black lw 0.5 front

# 把两侧的零点标在图上：对齐方式必须可见、可核对
set label 1 sprintf("zero: TB VBM=%.4f / VASP VBM=%.4f eV", tb_vbm, vasp_vbm) \
      at graph 0.012, 0.95 font "Arial,9"

set dashtype 2 (3,3)
plot "MoS2-pbe.txt" using 1:($2 - vasp_vbm) with lines lc "orange" lw 4 title "VASP PBE", \
     for [col=3:17] "MoS2-pbe.txt" using 1:(column(col) - vasp_vbm) with lines lc "orange" lw 4 notitle, \
     "data/band_nn.dat" using (map($1)):($2 - tb_vbm) with lines lc "blue" lw 3 dt 2 title "TB-NN", \
     "data/band_nn.dat" using (map($1)):($3 - tb_vbm) with lines lc "blue" lw 3 dt 2 notitle, \
     "data/band_nn.dat" using (map($1)):($4 - tb_vbm) with lines lc "blue" lw 3 dt 2 notitle
