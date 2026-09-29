#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxzm3b_c( const int m,
               const int l,
                     double *restrict x )

{
    
    // FFT adjustement
    const int l3 = l / 3;
    
    // Walking pointer difference between real and imag part
    const ptrdiff_t step = vlen4 * m;
    
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
    __td r0re, r0im, r1re, r1im, r2re, r2im, r01, r02, r03, r04;
    
    for ( int i3 = 0; i3 < l3; i3++ ) {
        
        for ( int i2 = 0; i2 < m; i2++ ) {
            
            for ( int i1 = 0; i1 < 4; i1++ ) {
                
                r1re = _t_load_pd( px1re );
                r2re = _t_load_pd( px2re );
                
                r01 = _t_sub_pd( r1re, r2re );
                r03 = _t_add_pd( r1re, r2re );
                
                r1im = _t_load_pd( px1im );
                r2im = _t_load_pd( px2im );
                
                r02 = _t_sub_pd( r1im, r2im );
                r04 = _t_add_pd( r1im, r2im );
                
                r0re = _t_load_pd( px0re );
                r0im = _t_load_pd( px0im );
                
                #if defined (__FMA__)
                r1re = _t_fmadd_pd( rC31, r03, r0re );
                r1im = _t_fmadd_pd( rC31, r04, r0im );
                
                r2re = _t_fmadd_pd(  rC32, r02, r1re );
                r2im = _t_fnmadd_pd( rC32, r01, r1im );
                #else
                r1re = _t_mul_pd( rC31, r03 );
                r1im = _t_mul_pd( rC31, r04 );
                r2re = _t_mul_pd( rC32, r02 );
                r2im = _t_mul_pd( rC32, r01 );
                
                r1re = _t_add_pd( r1re, r0re );
                r1im = _t_add_pd( r1im, r0im );
                
                r2re = _t_add_pd( r1re, r2re );
                r2im = _t_sub_pd( r1im, r2im );
                #endif
                
                r0re = _t_add_pd( r0re, r03 );
                r0im = _t_add_pd( r0im, r04 );
                
                _t_store_pd( px0re, r0re );
                _t_store_pd( px0im, r0im );
                
                r1re = _t_add_pd( r1re, r1re );
                r1im = _t_add_pd( r1im, r1im );
                
                _t_store_pd( px2re, r2re );
                _t_store_pd( px2im, r2im );
                
                r1re = _t_sub_pd( r1re, r2re );
                r1im = _t_sub_pd( r1im, r2im );
                
                _t_store_pd( px1re, r1re );
                _t_store_pd( px1im, r1im );
                
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
    
}