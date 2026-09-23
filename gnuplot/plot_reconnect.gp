# 参考能带重连检查
#   上：灰点 = 原始"按能量排序"的能级；彩色粗线 = 重连出的每条物理能带
#   下：每条曲线在各 k 取自哪个排序块（10..14）—— 台阶跳变处就是"换身份"的位置
if (!exists("MODE")) MODE = 1
OUTFILE = "img/GW-reconnect.png"
# pngcairo 需要 libcairo，缺了可以用 -e "GTERM='png'" 换成老终端（见 build.sh）
if (!exists("GTERM")) GTERM = "pngcairo"

if (MODE == 1) {
  set terminal @GTERM font "Arial,11" size 950,780 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,11" size 950,780 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}

C1 = "#1f77b4"
C2 = "#d62728"
C3 = "#000000"
C4 = "#9467bd"
C5 = "#2ca02c"

KX = 1.31633
MX = 1.97450
GX = 3.11448
set xrange [0:GX]
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set grid
set key top right
set multiplot layout 2,1 title "参考能带重连：原始(排序) vs 重连(物理)"

# ---- 上图：能量 ----
set ylabel "Energy (eV)"
set yrange [0.8:4.6]
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 1 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 1 front
set label 1 "锚点 (K)" at KX, graph 0.03 center font "Arial,9"
plot "data/gw_reconnect.dat" using 1:2 w p pt 7 ps 0.35 lc rgb "#999999" title "排序块（原始）", \
     "data/gw_reconnect.dat" using 1:3 w p pt 7 ps 0.35 lc rgb "#999999" notitle, \
     "data/gw_reconnect.dat" using 1:4 w p pt 7 ps 0.35 lc rgb "#999999" notitle, \
     "data/gw_reconnect.dat" using 1:5 w p pt 7 ps 0.35 lc rgb "#999999" notitle, \
     "data/gw_reconnect.dat" using 1:6 w p pt 7 ps 0.35 lc rgb "#999999" notitle, \
     "data/gw_reconnect.dat" using 1:7  w l lw 2.5 lc rgb C1 title "物理带 1", \
     "data/gw_reconnect.dat" using 1:8  w l lw 2.5 lc rgb C2 title "物理带 2", \
     "data/gw_reconnect.dat" using 1:9  w l lw 2.5 lc rgb C3 title "物理带 3", \
     "data/gw_reconnect.dat" using 1:10 w l lw 2.5 lc rgb C4 title "物理带 4", \
     "data/gw_reconnect.dat" using 1:11 w l lw 2.5 lc rgb C5 title "物理带 5"

# ---- 下图：每条曲线取自哪个排序块 ----
set ylabel "取自排序块号"
unset label 1
set yrange [9.4:14.6]
set ytics 10,1,14
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 1 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 1 front
set key top right
plot "data/gw_reconnect.dat" using 1:12 w steps lw 2 lc rgb C1 title "物理带 1", \
     "data/gw_reconnect.dat" using 1:13 w steps lw 2 lc rgb C2 title "物理带 2", \
     "data/gw_reconnect.dat" using 1:14 w steps lw 2 lc rgb C3 title "物理带 3", \
     "data/gw_reconnect.dat" using 1:15 w steps lw 2 lc rgb C4 title "物理带 4", \
     "data/gw_reconnect.dat" using 1:16 w steps lw 2 lc rgb C5 title "物理带 5"
unset multiplot
