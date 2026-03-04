#!/bin/bash
###########################################################################################
####                                     IFX SET UP                                    ####
###########################################################################################
fcompile="ifx -fast \
              -ffast-math \
              -fomit-frame-pointer \
              -fno-stack-protector \
              -xHost \
              -qopt-zmm-usage=high \
              -qopenmp \
              -diag-disable=8711 \
              -D$memory \
              -cpp"