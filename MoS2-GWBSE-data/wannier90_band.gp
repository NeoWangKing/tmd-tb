set terminal qt font "Arial,12" size 800,600 enhanced

# set terminal pngcairo font "Arial,12" size 800,600 enhanced
# set output "img/TB-NN.png"

set style data dots
set nokey

GAMMA1 = 0
K = 1.31633
M = 1.97450
GAMMA2 = 3.11448

set xlabel "k-path"
set ylabel "Energy (eV)"
# set xrange [GAMMA1: GAMMA2]
# set yrange [-15.81492 : 13.62639]
# set yrange [-3:5]

set arrow from K, graph 0 to K, graph 1 nohead
set arrow from M, graph 0 to M, graph 1 nohead

set xtics ("GAMMA" GAMMA1,"K" K,"M" M,"GAMMA" GAMMA2)

 plot "wannier90_band.dat" using 1:2 title "wannier_band"
