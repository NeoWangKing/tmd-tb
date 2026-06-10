# 设置终端
set terminal qt font "Arial,12" size 800,600
# 若需要保存为 PNG，取消下面两行的注释，并注释上面的 qt 行
# set terminal pngcairo font "Arial,12" size 800,600
# set output "TB-NN-NNN-TNN-SOC-SPIN.png"

# 映射函数（与 C 代码输出的原始横坐标对应）
map(x) = (x <= 0.577350) ? 0 + (x-0)*(0.180988-0)/(0.577350-0) : \
         (x <= 1.244017) ? 0.180988 + (x-0.577350)*(0.389974-0.180988)/(1.244017-0.577350) : \
         (x <= 1.577350) ? 0.389974 + (x-1.244017)*(0.494467-0.389974)/(1.577350-1.244017) : \
         (x <= 1.910684) ? 0.494467 + (x-1.577350)*(0.598960-0.494467)/(1.910684-1.577350) : 0.598960

M  = map(0.0)
G  = map(0.577350)
K  = map(1.244017)
Mp = map(1.577350)
Kp = map(1.910684)

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

# 从注释行读取 TB 的 Gamma 点价带顶
c_vbm = -0.193000

# 读取 VASP 数据文件的价带顶（这里 VASP 有 16 条带，我们取所有带的最大值）
stats "MoS2-pbe.txt" using 10 nooutput
vasp_vbm = STATS_max

plot "MoS2-pbe.txt" using 1:($2 - vasp_vbm) with lines lc "orange" lw 4 title "VASP PBE", \
     for [col=3:17] "MoS2-pbe.txt" using 1:(column(col) - vasp_vbm) with lines lc "orange" lw 4 notitle, \
     "band_nn.dat" using (map($1)):($2 - c_vbm) with lines lc "blue" lw 3 dt 2 title "TB-NN", \
     ""           using (map($1)):($3 - c_vbm) with lines lc "blue" lw 3 dt 2 notitle, \
     ""           using (map($1)):($4 - c_vbm) with lines lc "blue" lw 3 dt 2 notitle
