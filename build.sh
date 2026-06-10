#!/bin/sh

mkdir -p bin
mkdir -p data
mkdir -p img

set -xe

clang -Wall -Wextra -o bin/main_nn main.c \
      -I/opt/homebrew/opt/lapack/include \
      -L/opt/homebrew/opt/lapack/lib \
      -llapacke -llapack -lblas \
      -lm -std=c99
./bin/main_nn > ./data/band_nn.dat
gnuplot -persist gnuplot/plot_nn.gp

clang -Wall -Wextra -o bin/main_tnn main.c \
      -I/opt/homebrew/opt/lapack/include \
      -L/opt/homebrew/opt/lapack/lib \
      -DNNN_MODEL=1 \
      -llapacke -llapack -lblas \
      -lm -std=c99
./bin/main_tnn > ./data/band_tnn.dat
gnuplot -persist gnuplot/plot_tnn.gp

clang -Wall -Wextra -o bin/main_soc main.c \
      -I/opt/homebrew/opt/lapack/include \
      -L/opt/homebrew/opt/lapack/lib \
      -DNNN_MODEL=1 \
      -DSOC_MODEL=1 \
      -llapacke -llapack -lblas \
      -lm -std=c99
./bin/main_soc > ./data/band_soc.dat
gnuplot -persist gnuplot/plot_soc.gp

clang -Wall -Wextra -o bin/main_soc_spin main.c \
      -I/opt/homebrew/opt/lapack/include \
      -L/opt/homebrew/opt/lapack/lib \
      -DNNN_MODEL=1 \
      -DSOC_MODEL=1 \
      -DSOC_SPIN=1 \
      -llapacke -llapack -lblas \
      -lm -std=c99
./bin/main_soc_spin > ./data/band_soc_spin.dat
gnuplot -persist gnuplot/plot_soc_spin.gp
