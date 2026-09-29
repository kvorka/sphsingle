#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxzm2b_c( const int m,
               const int l,
               double *restrict x )

{
    
    // FFT adjustement
    const int l2 = l / 2;
    
    // Walking pointer difference between real and imag part
    const ptrdiff_t step = vlen4 * m;
    
    // Memory addresses
    double *restrict px0re = x + step * ( 0 + 2 * l2 * 0 );
    double *restrict px0im = x + step * ( 1 + 2 * l2 * 0 );
    double *restrict px1re = x + step * ( 0 + 2 * l2 * 1 );
    double *restrict px1im = x + step * ( 1 + 2 * l2 * 1 );
    
    // Registers to be used
    __td r0re, r0im, r1re, r1im, r01, r02;
    
    for ( int i3 = 0; i3 < l2; i3++ ) {
        
        for ( int i2 = 0; i2 < m; i2++ ) {
            
            for ( int i1 = 0; i1 < 4; i1++ ) {
                
                r0re = _t_load_pd( px0re );
                r0im = _t_load_pd( px0im );
                r1re = _t_load_pd( px1re );
                r1im = _t_load_pd( px1im );
                
                r1re = _t_sub_pd( r0re, r1re );
                r1im = _t_sub_pd( r0im, r1im );
                r01  = _t_add_pd( r0re, r0re );
                r02  = _t_add_pd( r0im, r0im );
                
                r0re = _t_sub_pd( r01, r1re );
                r0im = _t_sub_pd( r02, r1im );
                
                _t_store_pd( px0re, r0re );
                _t_store_pd( px0im, r0im );
                _t_store_pd( px1re, r1re );
                _t_store_pd( px1im, r1im );
                
                // Walking to next SIMD line before next
                // i1 cycle iteration.
                px0re += vlen;
                px0im += vlen;
                px1re += vlen;
                px1im += vlen;
                
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
        
    }
    
}