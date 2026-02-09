#!/bin/bash
###########################################################################################
####                                  GFORTRAN SET UP                                  ####
###########################################################################################
fcompile="gfortran-12 -Ofast \
                      -march=native \
                      -mno-vzeroupper \
                      -mprefer-vector-width=512 \
                      -fno-bounds-check \
                      -fomit-frame-pointer \
                      -fno-stack-protector \
                      -fdefault-real-8 \
                      -flto=auto \
                      -fwhole-program \
                      -fopenmp \
                      -D$memory \
                      -cpp"