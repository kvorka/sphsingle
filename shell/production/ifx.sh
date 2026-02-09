#!/bin/bash
###########################################################################################
####                                     IFX SET UP                                    ####
###########################################################################################
fcompile="ifx -fast \
              -Ofast \
              -ipo \
              -ffast-math \
              -funroll-loops \
              -fomit-frame-pointer \
              -fno-stack-protector \
              -xHost \
              -qopt-zmm-usage=high \
              -qopenmp \
              -D$memory \
              -cpp"