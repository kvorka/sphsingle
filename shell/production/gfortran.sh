#!/bin/bash
###########################################################################################
####                                  GFORTRAN SET UP                                  ####
###########################################################################################
fcompile="gfortran-12 -Ofast \
                      -march=native \
                      -mfma \
                      -mno-vzeroupper \
                      -mprefer-vector-width=512 \
                      -fno-bounds-check \
                      -fargument-noalias-global \
                      -fstrict-aliasing \
                      -fomit-frame-pointer \
                      -fno-stack-protector \
                      -fdefault-real-8 \
                      -funroll-loops \
                      -flto=auto \
                      -fwhole-program \
                      -fopenmp \
                      -D$memory \
                      -cpp"