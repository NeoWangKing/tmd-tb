# 拟合诊断：GW 参考能带（重连成物理能带后）与当前 TB 的残差曲线
# 数据由 bin/fit_main 生成（读 data/band_gw.dat 与 GW 参考能带）
if (!exists("MODE")) MODE = 1
OUTFILE = "img/GW-residual.png"

if (MODE == 1) {
  set terminal pngcairo font "Arial,12" size 900,600 enhanced
  set output OUTFILE
  print "MODE = 1 -> 输出 ", OUTFILE
} else {
  set terminal qt persist font "Arial,12" size 900,600 enhanced
  print "MODE = 0 -> 弹出 Qt 窗口（看不到窗口 = 本机 Qt 终端不可用，请用 -e \"MODE=1\"）"
}

# 原始残差均值（= 最优刚性平移量），从数据里算，不写死
stats "data/fit_residual.dat" using 2 nooutput; mVB  = STATS_mean
stats "data/fit_residual.dat" using 3 nooutput; mCB1 = STATS_mean
stats "data/fit_residual.dat" using 4 nooutput; mCB2 = STATS_mean

KX = 1.31633
MX = 1.97450
GX = 3.11448

set xlabel "k-path (Å^{-1})"
set ylabel "E_GW − E_TB − Δ   (eV)"
set xrange [0:GX]
set yrange [-0.35:0.35]
set grid
set xtics ("Γ" 0, "K" KX, "M" MX, "Γ" GX)
set arrow from KX, graph 0 to KX, graph 1 nohead lc rgb "black" lw 0.5 front
set arrow from MX, graph 0 to MX, graph 1 nohead lc rgb "black" lw 0.5 front
set arrow from 0, 0 to GX, 0 nohead lc rgb "gray60" dt 2
set key top right
set label 1 sprintf("零线=平移后;  原始残差均值 VB %.2f / CB1 %.2f / CB2 %.2f eV", mVB, mCB1, mCB2) \
      at graph 0.015, 0.95 font "Arial,9"

set dashtype 2 (4,3)
plot "data/fit_residual.dat" using 1:5 w l lc "blue"  lw 2.5 title "VB  平移后", \
     "data/fit_residual.dat" using 1:6 w l lc "red"   lw 2.5 title "CB1 平移后", \
     "data/fit_residual.dat" using 1:7 w l lc "green" lw 2.0 title "CB2 平移后", \
     "data/fit_residual.dat" using 1:8 w l lc "blue"  lw 1.0 dt 2 notitle, \
     "data/fit_residual.dat" using 1:9 w l lc "red"   lw 1.0 dt 2 notitle, \
     "data/fit_residual.dat" using 1:10 w l lc "green" lw 1.0 dt 2 notitle
