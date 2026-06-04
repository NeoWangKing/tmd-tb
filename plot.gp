# 设置终端
set terminal qt font "Arial,12" size 800,600
# 若需要保存为 PNG，取消下面两行的注释，并注释上面的 qt 行
# set terminal pngcairo font "Arial,12" size 800,600
# set output "band.png"

# 定义对称点坐标（根据你的程序输出）
M  = 0.0
G  = 0.577350
K  = 1.244017
Mp = 1.577350
Kp = 1.910684

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

# 画横线
set arrow from graph 0,  0 to graph 1,  0 nohead lc black lw 0.5 front

# 绘图：三条能带全部用蓝色（lc "blue"），线宽 3（lw 3）
plot "band.dat" using 1:2 with lines lc "blue" lw 3 title "E1", \
     "band.dat" using 1:3 with lines lc "blue" lw 3 title "E2", \
     "band.dat" using 1:4 with lines lc "blue" lw 3 title "E3"
