#!/bin/sh

set -xe

mkdir -p ./build/
rm -rf ./build/
mkdir -p ./build/

clang -Wall -Wextra -o main main.c \
      -I/opt/homebrew/opt/lapack/include \
      -L/opt/homebrew/opt/lapack/lib \
      -llapacke -llapack -lblas \
      -lm -std=c99
./main > band.dat
./main

gnuplot -persist plot.gp
