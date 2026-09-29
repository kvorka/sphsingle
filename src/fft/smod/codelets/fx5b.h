#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxzm5b_c( const int m,
               const int l,
                     double *restrict x )

{
    
    // FFT adjustement
    const int l5 = l / 5;
    
    // Walking pointer difference between real and imag part
    const ptrdiff_t step  = vlen4 * m;
    
    // Memory addresses
    double *restrict px0re = x + step * ( 0 + 2 * l5 * 0 );
    double *restrict px0im = x + step * ( 1 + 2 * l5 * 0 );
    double *restrict px1re = x + step * ( 0 + 2 * l5 * 1 );
    double *restrict px1im = x + step * ( 1 + 2 * l5 * 1 );
    double *restrict px2re = x + step * ( 0 + 2 * l5 * 2 );
    double *restrict px2im = x + step * ( 1 + 2 * l5 * 2 );
    double *restrict px3re = x + step * ( 0 + 2 * l5 * 3 );
    double *restrict px3im = x + step * ( 1 + 2 * l5 * 3 );
    double *restrict px4re = x + step * ( 0 + 2 * l5 * 4 );
    double *restrict px4im = x + step * ( 1 + 2 * l5 * 4 );
    
    // FFT constants
    const __td rC51 = _t_set1_pd( +0.2500000000000000000 );
    const __td rC52 = _t_set1_pd( +0.5590169943749474241 );
    const __td rC53 = _t_set1_pd( +0.6180339887498948482 );
    const __td rC54 = _t_set1_pd( -0.9510565162951535721 );
    
    // Registers to be used
    __td r0re, r0im, r1re, r1im, r2re, r2im, r3re, r3im, r4re, r4im, r01, r02;
    
    for ( int i3 = 0; i3 < l5; i3++ ) {
            
        for ( int i2 = 0; i2 < m; i2++ ) {
            
            for ( int i1 = 0; i1 < 4; i1++ ) {
                
                r0im = _t_load_pd( px1re );
                r0re = _t_load_pd( px4re );
                
                r1re = _t_add_pd( r0im, r0re );
                r4re = _t_sub_pd( r0im, r0re );
                
                r2im = _t_load_pd( px1im );
                r2re = _t_load_pd( px4im );
                
                r1im = _t_add_pd( r2im, r2re );
                r4im = _t_sub_pd( r2im, r2re );
                
                r01 = _t_load_pd( px2re );
                r02 = _t_load_pd( px3re );
                
                r3re = _t_sub_pd( r01, r02 );
                r0re = _t_add_pd( r01, r02 );
                
                r2im = _t_load_pd( px2im );
                r2re = _t_load_pd( px3im );
                
                r3im = _t_sub_pd( r2im, r2re );
                r0im = _t_add_pd( r2im, r2re );
                
                #if defined (__FMA__)
                r2re = _t_fmadd_pd( rC53, r3re, r4re );
                r2im = _t_fmadd_pd( rC53, r3im, r4im );
                
                r3re = _t_fmsub_pd( rC53, r4re, r3re );
                r3im = _t_fmsub_pd( rC53, r4im, r3im );
                #else
                r2re = _t_mul_pd( rC53, r3re );
                r2im = _t_mul_pd( rC53, r3im );
                r01  = _t_mul_pd( rC53, r4re );
                r02  = _t_mul_pd( rC53, r4im );
                
                r2re = _t_add_pd( r2re, r4re );
                r2im = _t_add_pd( r2im, r4im );
                r3re = _t_sub_pd( r01, r3re );
                r3im = _t_sub_pd( r02, r3im );
                #endif
                
                r4re = _t_add_pd( r1re, r0re );
                r4im = _t_add_pd( r1im, r0im );
                
                r1re = _t_sub_pd( r1re, r0re );
                r1im = _t_sub_pd( r1im, r0im );
                
                r0re = _t_load_pd( px0re );
                r0im = _t_load_pd( px0im );
                
                r01 = _t_add_pd( r0re, r4re );
                r02 = _t_add_pd( r0im, r4im );
                
                _t_store_pd( px0re, r01 );
                _t_store_pd( px0im, r02 );
                
                #if defined (__FMA__)
                r0re = _t_fnmadd_pd( rC51, r4re, r0re );
                r0im = _t_fnmadd_pd( rC51, r4im, r0im );
                
                r1re = _t_fnmadd_pd( rC52, r1re, r0re );
                r1im = _t_fnmadd_pd( rC52, r1im, r0im );
                #else
                r01  = _t_mul_pd( rC51, r4re );
                r02  = _t_mul_pd( rC51, r4im );
                r1re = _t_mul_pd( rC52, r1re );
                r1im = _t_mul_pd( rC52, r1im );
                
                r0re = _t_sub_pd( r0re, r01 );
                r0im = _t_sub_pd( r0im, r02 );
                
                r1re = _t_sub_pd( r0re, r1re );
                r1im = _t_sub_pd( r0im, r1im );
                #endif
                
                r0re = _t_add_pd( r0re, r0re );
                r0im = _t_add_pd( r0im, r0im );
                
                r0re = _t_sub_pd( r0re, r1re );
                r0im = _t_sub_pd( r0im, r1im );
                
                #if defined (__FMA__)
                r3re = _t_fmadd_pd(  rC54, r3re, r1im );
                r3im = _t_fnmadd_pd( rC54, r3im, r1re );
                r2im = _t_fnmadd_pd( rC54, r2im, r0re );
                r2re = _t_fmadd_pd(  rC54, r2re, r0im );
                #else
                r3re = _t_mul_pd( rC54, r3re );
                r3im = _t_mul_pd( rC54, r3im );
                r2im = _t_mul_pd( rC54, r2im );
                r2re = _t_mul_pd( rC54, r2re );
                
                r3re = _t_add_pd( r1im, r3re );
                r3im = _t_sub_pd( r1re, r3im );
                r2im = _t_sub_pd( r0re, r2im );
                r2re = _t_add_pd( r0im, r2re );
                #endif
                
                _t_store_pd( px4re, r2im );
                _t_store_pd( px4im, r2re );
                
                r1re = _t_add_pd( r1re, r1re );
                r1im = _t_add_pd( r1im, r1im );
                
                _t_store_pd( px3re, r3im );
                _t_store_pd( px3im, r3re );
                
                r0re = _t_add_pd( r0re, r0re );
                r0im = _t_add_pd( r0im, r0im );
                
                r1re = _t_sub_pd( r1re, r3im );
                r1im = _t_sub_pd( r1im, r3re );
                
                _t_store_pd( px2re, r1re );
                _t_store_pd( px2im, r1im );
                
                r0re = _t_sub_pd( r0re, r2im );
                r0im = _t_sub_pd( r0im, r2re );
                
                _t_store_pd( px1re, r0re );
                _t_store_pd( px1im, r0im );
                
                // Walking to next SIMD line before next
                // i1 cycle iteration.
                px0re += vlen;
                px0im += vlen;
                px1re += vlen;
                px1im += vlen;
                px2re += vlen;
                px2im += vlen;
                px3re += vlen;
                px3im += vlen;
                px4re += vlen;
                px4im += vlen;
                
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
        px3re += step;
        px3im += step;
        px4re += step;
        px4im += step;
        
    }
    
}