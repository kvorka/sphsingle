#!/bin/bash
###########################################################################################
####                                 COMPILER SET UP                                   ####
###########################################################################################
case $compiler in
    ifort)
        source ./shell/production/ifort.sh
    ;;
    
    ifx)
        source ./shell/production/ifx.sh
    ;;
    
    gfortran)
        source ./shell/production/gfortran.sh
    ;;
    
    gdebug)
        source ./shell/debug/gdebug.sh
    ;;
esac

###########################################################################################
####                                  COMPILER FUNCTION                                ####
###########################################################################################
function libfcompile() {
    $fcompile -c $1/$2
    $fcompile -c $(find $1/smod/. -type f)
}

###########################################################################################
####                                  CLEANING FUNCTION                                ####
###########################################################################################
function libfclean() {
    rm *.smod *.mod *.o || true
}