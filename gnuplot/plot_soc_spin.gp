# 设置终端
# set terminal qt font "Arial,12" size 800,600 enhanced
# 若需要保存为 PNG，取消下面两行的注释，并注释上面的 qt 行
set terminal pngcairo font "Arial,12" size 800,600 enhanced
set output "img/TB-NN-NNN-TNN-SOC-SPIN.png"

map(x) = (x <= 0.577350) ? 0 + (x-0)*(0.180988-0)/(0.577350-0) : \
         (x <= 1.244017) ? 0.180988 + (x-0.577350)*(0.389974-0.180988)/(1.244017-0.577350) : \
         (x <= 1.577350) ? 0.389974 + (x-1.244017)*(0.494467-0.389974)/(1.577350-1.244017) : \
         (x <= 1.910684) ? 0.494467 + (x-1.577350)*(0.598960-0.494467)/(1.910684-1.577350) : 0.598960

# 定义对称点坐标（根据你的程序输出）
M  = map(0.0)
G  = map(0.577350)
K  = map(1.244017)
Mp = map(1.577350)
Kp = map(1.910684)

# ---------------------------------------------------------------------------
# 零点对齐：价带条数由"电子数 / 模型"决定，两侧用同一条规则取价带顶的极大值。
#   MoS2 单层有 18 个价电子 -> VASP 侧 9 条价带；含 SOC 时每带 2 个自旋分量，
#   所以 ncl 文件里价带顶是第 2*9 = 18 条能级（列号 = 1 + 2*NV_VASP = 19）。
#   本 TB 模型（含 SOC）：2 条价带 + 4 条导带。band_soc_spin.dat 的两组自旋是
#   各自升序排列的（第 2~4 列 spin-up，第 5~7 列 spin-down），所以价带顶取
#   两组里第 NV_TB 条的较大者（gnuplot 没有 max() 函数，用三元运算符写）。
# 两个坑：
#   1) 必须写 column(1+NV)，写 (1+NV) 会被当成"数值"而不是"第几列"；
#   2) stats 要放在 set xrange/yrange 之前，否则数据点会被坐标范围过滤掉。
# ---------------------------------------------------------------------------
NV_VASP    = 9
NV_TB      = 2      # 模型里的价带条数（含自旋）：2 条价带 + 4 条导带
NV_TB_SPIN = 1      # 每个自旋块内部的价带条数 = NV_TB/2
                    # （band_soc_spin.dat 里两组自旋各自升序，所以块内第 1 条才是价带）

stats "data/band_soc_spin.dat" \
      using ((column(1+NV_TB_SPIN) > column(4+NV_TB_SPIN)) \
             ? column(1+NV_TB_SPIN) : column(4+NV_TB_SPIN)) nooutput
tb_vbm   = STATS_max
stats "MoS2-pbe-ncl.txt" using (column(1+2*NV_VASP)) nooutput
vasp_vbm = STATS_max

# 设置坐标轴标签和范围
set xlabel "k-path"
set ylabel "Energy (eV)"
set xrange [M:Kp]
set yrange [-1:5]
set grid

# 设置 x 轴刻度标签
set xtics ("M" 0, "Γ" G, "K" K, "M'" Mp, "K'" Kp)

# 画竖线（使用 arrow 样式，无箭头）
set arrow from M,  graph 0 to M,  graph 1 nohead lc black lw 0.5 front
set arrow from G,  graph 0 to G,  graph 1 nohead lc black lw 0.5 front
set arrow from K,  graph 0 to K,  graph 1 nohead lc black lw 0.5 front
set arrow from Mp, graph 0 to Mp, graph 1 nohead lc black lw 0.5 front
set arrow from Kp, graph 0 to Kp, graph 1 nohead lc black lw 0.5 front

# 把两侧的零点标在图上：对齐方式必须可见、可核对
set label 1 sprintf("zero: TB VBM=%.4f / VASP VBM=%.4f eV", tb_vbm, vasp_vbm) \
      at graph 0.012, 0.95 font "Arial,9"

set dashtype 2 (3,3)
plot "MoS2-pbe-ncl.txt" using 1:($2 - vasp_vbm) with lines lc "orange" lw 4 title "VASP PBE", \
     for [col=3:33] "MoS2-pbe-ncl.txt" using 1:(column(col) - vasp_vbm) with lines lc "orange" lw 4 notitle, \
     "data/band_soc_spin.dat" using (map($1)):($2 - tb_vbm) with lines lc "red" lw 3 dt 2 title "spin up", \
     "data/band_soc_spin.dat" using (map($1)):($3 - tb_vbm) with lines lc "red" lw 3 dt 2 notitle, \
     "data/band_soc_spin.dat" using (map($1)):($4 - tb_vbm) with lines lc "red" lw 3 dt 2 notitle, \
     "data/band_soc_spin.dat" using (map($1)):($5 - tb_vbm) with lines lc "blue" lw 3 dt 2 title "spin down", \
     "data/band_soc_spin.dat" using (map($1)):($6 - tb_vbm) with lines lc "blue" lw 3 dt 2 notitle, \
     "data/band_soc_spin.dat" using (map($1)):($7 - tb_vbm) with lines lc "blue" lw 3 dt 2 notitle
