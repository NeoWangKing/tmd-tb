#!/bin/sh

set -xe

clang -Wall -Wextra -o main_nn main.c \
      -I/opt/homebrew/opt/lapack/include \
      -L/opt/homebrew/opt/lapack/lib \
      -llapacke -llapack -lblas \
      -lm -std=c99
./main_nn > band_nn.dat

gnuplot -persist plot_nn.gp

# clang -Wall -Wextra -o main_tnn main.c \
#       -I/opt/homebrew/opt/lapack/include \
#       -L/opt/homebrew/opt/lapack/lib \
#       -DNNN_MODEL=1 \
#       -llapacke -llapack -lblas \
#       -lm -std=c99
# ./main_tnn > band_tnn.dat
#
# gnuplot -persist plot_tnn.gp

# clang -Wall -Wextra -o main_soc main.c \
#       -I/opt/homebrew/opt/lapack/include \
#       -L/opt/homebrew/opt/lapack/lib \
#       -DNNN_MODEL=1 \
#       -DSOC_MODEL=1 \
#       -llapacke -llapack -lblas \
#       -lm -std=c99
# ./main_soc > band_soc.dat
#
# gnuplot -persist plot_soc.gp

# clang -Wall -Wextra -o main_soc_spin main.c \
#       -I/opt/homebrew/opt/lapack/include \
#       -L/opt/homebrew/opt/lapack/lib \
#       -DNNN_MODEL=1 \
#       -DSOC_MODEL=1 \
#       -DSOC_SPIN=1 \
#       -llapacke -llapack -lblas \
#       -lm -std=c99
# ./main_soc_spin > band_soc_spin.dat
#
# gnuplot -persist plot_soc_spin.gp
