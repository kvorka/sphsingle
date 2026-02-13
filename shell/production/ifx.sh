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
              -fvec-remainder-loops \
              -flto \
              -xHost \
              -qopt-zmm-usage=high \
              -assume contiguous_assumed_shape \
              -assume contiguous_pointer \
              -assume nodummy_aliases \
              -qopenmp \
              -D$memory \
              -cpp"