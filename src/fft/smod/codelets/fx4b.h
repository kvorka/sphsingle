#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxzm4b_c( const int m,
               const int l,
                     double *restrict x )

{
    
    // FFT adjustement
    const int l4 = l / 4;
    
    // Walking pointer difference between real and imag part
    const ptrdiff_t step  = vlen4 * m;
    
    // Memory addresses
    double *restrict px0re = x + step * ( 0 + 2 * l4 * 0 );
    double *restrict px0im = x + step * ( 1 + 2 * l4 * 0 );
    double *restrict px1re = x + step * ( 0 + 2 * l4 * 1 );
    double *restrict px1im = x + step * ( 1 + 2 * l4 * 1 );
    double *restrict px2re = x + step * ( 0 + 2 * l4 * 2 );
    double *restrict px2im = x + step * ( 1 + 2 * l4 * 2 );
    double *restrict px3re = x + step * ( 0 + 2 * l4 * 3 );
    double *restrict px3im = x + step * ( 1 + 2 * l4 * 3 );
    
    // Registers to be used
    __td r0re, r0im, r1re, r1im, r2re, r2im, r3re, r3im, r01, r02, r03, r04;
    
    for ( int i3 = 0; i3 < l4; i3++ ) {
            
        for ( int i2 = 0; i2 < m; i2++ ) {
            
            for ( int i1 = 0; i1 < 4; i1++ ) {
                
                r01 = _t_load_pd( px0re );
                r03 = _t_load_pd( px2re );
                
                r0re = _t_add_pd( r01, r03 );
                r2re = _t_sub_pd( r01, r03 );
                
                r02 = _t_load_pd( px0im );
                r04 = _t_load_pd( px2im );
                
                r0im = _t_add_pd( r02, r04 );
                r2im = _t_sub_pd( r02, r04 );
                
                r01 = _t_load_pd( px1re );
                r03 = _t_load_pd( px3re );
                
                r1re = _t_add_pd( r01, r03 );
                r3re = _t_sub_pd( r01, r03 );
                
                r02 = _t_load_pd( px1im );
                r04 = _t_load_pd( px3im );
                
                r1im = _t_add_pd( r02, r04 );
                r3im = _t_sub_pd( r02, r04 );
                
                r1re = _t_sub_pd( r0re, r1re );
                r1im = _t_sub_pd( r0im, r1im );
                r01  = _t_add_pd( r0re, r0re );
                r02  = _t_add_pd( r0im, r0im );
                
                _t_store_pd( px2re, r1re );
                _t_store_pd( px2im, r1im );
                
                r0re = _t_sub_pd( r01,  r1re );
                r0im = _t_sub_pd( r02,  r1im );
                r3im = _t_sub_pd( r2re, r3im );
                r3re = _t_add_pd( r2im, r3re );
                
                r01  = _t_add_pd( r2re, r2re );
                r02  = _t_add_pd( r2im, r2im );
                
                _t_store_pd( px0re, r0re );
                _t_store_pd( px0im, r0im );
                _t_store_pd( px1re, r3im );
                _t_store_pd( px1im, r3re );
                
                r2re = _t_sub_pd( r01, r3im );
                r2im = _t_sub_pd( r02, r3re );
                
                _t_store_pd( px3re, r2re );
                _t_store_pd( px3im, r2im );
                
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
        
    }
    
}