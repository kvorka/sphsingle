#!/bin/bash
###########################################################################################
####                                    IFORT SET UP                                   ####
###########################################################################################
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
