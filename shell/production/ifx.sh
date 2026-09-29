#!/bin/bash
###########################################################################################
####                                     IFX SET UP                                    ####
###########################################################################################
if [[ "$omp" == "enabled" ]];
  then
    fcompile="ifx -fast \
                  -ffast-math \
                  -fomit-frame-pointer \
                  -fno-stack-protector \
                  -xHost \
                  -qopt-zmm-usage=high \
                  -qno-openmp-simd \
                  -diag-disable=8711 \
                  -D$memory \
                  -cpp"
  else
    fcompile="ifx -fast \
                  -ffast-math \
                  -fomit-frame-pointer \
                  -fno-stack-protector \
                  -xHost \
                  -qopt-zmm-usage=high \
                  -diag-disable=8711 \
                  -D$memory \
                  -cpp"
fi