#!/bin/bash
###########################################################################################
####                                    IFORT SET UP                                   ####
###########################################################################################
if [[ "$omp" == "enabled" ]];
  then
    fcompile="ifort -fast \
                    -ffast-math \
                    -fomit-frame-pointer \
                    -fno-stack-protector \
                    -xHost \
                    -qopt-zmm-usage=high \
                    -qopenmp \
                    -diag-disable=10448 \
                    -diag-disable=11021 \
                    -D$memory \
                    -cpp"
  else
    fcompile="ifort -fast \
                    -ffast-math \
                    -fomit-frame-pointer \
                    -fno-stack-protector \
                    -xHost \
                    -qopt-zmm-usage=high \
                    -qopenmp-simd \
                    -diag-disable=10448 \
                    -diag-disable=11021 \
                    -D$memory \
                    -cpp"
fi
