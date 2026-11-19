# 幻灯片用图：(a) 物理参考能带　(b) 曲率拼接后的能带
#   数据 data/soc_stitch_compare.dat 由 bin/soc_check 写出
#   col2-7  = 参考的 6 条自旋通道（Up1-3, Dn1-3）
#   col8-13 = 拼接后的 6 条带
set terminal pngcairo font "Arial,20" size 1000,380 enhanced
set output "soc_stitch_cmp.png"

D = "../../data/soc_stitch_compare.dat"
C1 = "#1f77b4"
C2 = "#ff7f0e"
C3 = "#2ca02c"
C4 = "#d62728"
C5 = "#9467bd"
C6 = "#8c564b"

set multiplot layout 1,2 spacing 1.15,1.0

set xrange [0:1.9107]
set xtics ("M" 0, "Γ" 0.577350, "K" 1.244017, "M'" 1.577350, "K'" 1.910684) font "Arial,17"
set yrange [-1.2:4.2]
set ytics 2
set key off
set grid ytics
set arrow 1 from 0.577350, graph 0 to 0.577350, graph 1 nohead lc rgb "#cccccc" lw 1
set arrow 2 from 1.244017, graph 0 to 1.244017, graph 1 nohead lc rgb "#cccccc" lw 1
set arrow 3 from 1.577350, graph 0 to 1.577350, graph 1 nohead lc rgb "#cccccc" lw 1

set title "(a) 物理参考能带" font "Arial,20"
set ylabel "E (eV)" font "Arial,18"
plot D using 1:2 w l lw 2.2 lc rgb C1, \
     D using 1:3 w l lw 2.2 lc rgb C2, \
     D using 1:4 w l lw 2.2 lc rgb C3, \
     D using 1:5 w l lw 2.2 lc rgb C4, \
     D using 1:6 w l lw 2.2 lc rgb C5, \
     D using 1:7 w l lw 2.2 lc rgb C6

set title "(b) 曲率拼接后的能带" font "Arial,20"
unset ylabel
plot D using 1:8  w l lw 2.2 lc rgb C1, \
     D using 1:9  w l lw 2.2 lc rgb C2, \
     D using 1:10 w l lw 2.2 lc rgb C3, \
     D using 1:11 w l lw 2.2 lc rgb C4, \
     D using 1:12 w l lw 2.2 lc rgb C5, \
     D using 1:13 w l lw 2.2 lc rgb C6
unset multiplot
