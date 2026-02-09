#!/bin/bash
###########################################################################################
####                                   DEBUG SET UP                                    ####
###########################################################################################
fcompile="gfortran-12 -Og \
                      -fimplicit-none \
                      -Wall \
                      -Wline-truncation \
                      -Wcharacter-truncation \
                      -Wsurprising \
                      -Waliasing \
                      -Wimplicit-interface \
                      -Wunused-parameter \
                      -fwhole-file \
                      -fcheck=all \
                      -std=f2008 \
                      -pedantic \
                      -fbacktrace \
                      -D$memory \
                      -cpp"