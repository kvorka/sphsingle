#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxzm3a_c( const int m,
               const int k,
               const int l,
                     double *restrict x,
               const double *restrict t )

{
    
    // FFT adjustement
    const int l3 = l / 3;
    
    // Walking pointer difference between real and imag part
    const ptrdiff_t step  = vlen4 * m;
    const ptrdiff_t step2 = vlen4 * m * l3 * 4;
    
    // Memory addresses
    double *restrict px0re = x + step * ( 0 + 2 * l3 * 0 );
    double *restrict px0im = x + step * ( 1 + 2 * l3 * 0 );
    double *restrict px1re = x + step * ( 0 + 2 * l3 * 1 );
    double *restrict px1im = x + step * ( 1 + 2 * l3 * 1 );
    double *restrict px2re = x + step * ( 0 + 2 * l3 * 2 );
    double *restrict px2im = x + step * ( 1 + 2 * l3 * 2 );
    
    // FFT constants
    const __td rC31 = _t_set1_pd( -0.50000000000000000000 );
    const __td rC32 = _t_set1_pd( +0.86602540378443864676 );
    
    // Registers to be used
    __td rt1re, rt1im, rt2re, rt2im, r0re, r0im, r1re, r1im, r2re, r2im, r01, r02;
    
    for ( int i4 = 0; i4 < k; i4++ ) {
        
        rt1re = _t_set1_pd( *( t + 0 + 4 * i4 ) );
        rt1im = _t_set1_pd( *( t + 1 + 4 * i4 ) );
        rt2re = _t_set1_pd( *( t + 2 + 4 * i4 ) );
        rt2im = _t_set1_pd( *( t + 3 + 4 * i4 ) );
        
        for ( int i3 = 0; i3 < l3; i3++ ) {
            
            for ( int i2 = 0; i2 < m; i2++ ) {
                
                for ( int i1 = 0; i1 < 4; i1++ ) {
                    
                    r0re = _t_load_pd( px1re );
                    r0im = _t_load_pd( px1im );
                    
                    r01  = _t_mul_pd( rt1im, r0im );
                    r02  = _t_mul_pd( rt1im, r0re );
                    
                    r2re = _t_load_pd( px2re );
                    r2im = _t_load_pd( px2im );
                    
                    r1re = _t_mul_pd( rt2re, r2re );
                    r1im = _t_mul_pd( rt2re, r2im );
                    
                    #if defined (__FMA__)
                    r01  = _t_fmsub_pd(  rt1re, r0re, r01 );
                    r02  = _t_fmadd_pd(  rt1re, r0im, r02 );
                    r1re = _t_fnmadd_pd( rt2im, r2im, r1re );
                    r1im = _t_fmadd_pd(  rt2im, r2re, r1im );
                    #else
                    r0re = _t_mul_pd( rt1re, r0re );
                    r0im = _t_mul_pd( rt1re, r0im );
                    r2im = _t_mul_pd( rt2im, r2im );
                    r2re = _t_mul_pd( rt2im, r2re );
                    
                    r01  = _t_sub_pd( r0re, r01 );
                    r02  = _t_add_pd( r0im, r02 );
                    r1re = _t_sub_pd( r1re, r2im );
                    r1im = _t_add_pd( r1im, r2re );
                    #endif
                    
                    r1re = _t_sub_pd( r01, r1re );
                    r1im = _t_sub_pd( r02, r1im );
                    r2re = _t_add_pd( r01, r01 );
                    r2im = _t_add_pd( r02, r02 );
                    
                    r01  = _t_sub_pd( r2re, r1re );
                    r02  = _t_sub_pd( r2im, r1im );
                    r0re = _t_load_pd( px0re );
                    r0im = _t_load_pd( px0im );
                    
                    #if defined (__FMA__)
                    r2re = _t_fmadd_pd( rC31, r01, r0re );
                    r2im = _t_fmadd_pd( rC31, r02, r0im );
                    #else
                    r2re = _t_mul_pd( rC31, r01 );
                    r2im = _t_mul_pd( rC31, r02 );
                    
                    r2re = _t_add_pd( r2re, r0re );
                    r2im = _t_add_pd( r2im, r0im );
                    #endif
                    
                    r0re = _t_add_pd( r0re, r01 );
                    r0im = _t_add_pd( r0im, r02 );
                    
                    _t_store_pd( px0re, r0re );
                    _t_store_pd( px0im, r0im );
                    
                    #if defined (__FMA__)
                    r1im = _t_fmadd_pd(  rC32, r1im, r2re );
                    r1re = _t_fnmadd_pd( rC32, r1re, r2im );
                    #else
                    r1im = _t_mul_pd( rC32, r1im );
                    r1re = _t_mul_pd( rC32, r1re );
                    
                    r1im = _t_add_pd( r2re, r1im );
                    r1re = _t_sub_pd( r2im, r1re );
                    #endif
                    
                    r2re = _t_add_pd( r2re, r2re );
                    r2im = _t_add_pd( r2im, r2im );
                    
                    _t_store_pd( px2re, r1im );
                    _t_store_pd( px2im, r1re );
                    
                    r2re = _t_sub_pd( r2re, r1im );
                    r2im = _t_sub_pd( r2im, r1re );
                    
                    _t_store_pd( px1re, r2re );
                    _t_store_pd( px1im, r2im );
                    
                    // Walking to next SIMD line before next
                    // i1 cycle iteration.
                    px0re += vlen;
                    px0im += vlen;
                    px1re += vlen;
                    px1im += vlen;
                    px2re += vlen;
                    px2im += vlen;
                    
                }
                
                // No walking needed in here, because basic simd line,
                // i1 and even i2 are contiguous in memory.
                
            }
            
            // After i2 cycle, the address offset is already step, meaning 
            // px0re is where px0im initially started. Another move in addresses 
            // is required in order to move to next real/imag pair.
            px0re += step;
            px0im += step;
            px1re += step;
            px1im += step;
            px2re += step;
            px2im += step;
            
        }
        
        // After i3 cycle, the address offset is 2*l3*step. Overall step before
        // next stage needed is 6*l3*step, therefore more walking.
        px0re += step2;
        px0im += step2;
        px1re += step2;
        px1im += step2;
        px2re += step2;
        px2im += step2;
        
    }
    
}