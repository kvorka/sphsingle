#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxzm2a_c( const int m,
               const int k,
               const int l,
                     double *restrict x,
               const double *restrict t )

{
    
    // FFT adjustement
    const int l2 = l / 2;
    
    // Walking pointer difference between real and imag part
    const ptrdiff_t step  = vlen4 * m;
    const ptrdiff_t step2 = vlen4 * m * l2 * 2;
    
    // Memory addresses
    double *restrict px0re = x + step * ( 0 + 2 * l2 * 0 );
    double *restrict px0im = x + step * ( 1 + 2 * l2 * 0 );
    double *restrict px1re = x + step * ( 0 + 2 * l2 * 1 );
    double *restrict px1im = x + step * ( 1 + 2 * l2 * 1 );
    
    // Registers to be used
    __td rtre, rtim, r0re, r0im, r1re, r1im, r01, r02, r03, r04;
    
    for ( int i4 = 0; i4 < k; i4++ ) {
        
        rtre = _t_set1_pd( *( t + 0 + 2 * i4 ) );
        rtim = _t_set1_pd( *( t + 1 + 2 * i4 ) );
        
        for ( int i3 = 0; i3 < l2; i3++ ) {
            
            for ( int i2 = 0; i2 < m; i2++ ) {
                
                for ( int i1 = 0; i1 < 4; i1++ ) {
                    
                    r0re = _t_load_pd( px0re );
                    r0im = _t_load_pd( px0im );
                    r1re = _t_load_pd( px1re );
                    r1im = _t_load_pd( px1im );
                    
                    #if defined (__FMA__)
                    r03 = _t_fnmadd_pd( rtre, r1re, r0re );
                    r04 = _t_fnmadd_pd( rtim, r1re, r0im );
                    
                    r03  = _t_fmadd_pd(  rtim, r1im, r03 );
                    r04  = _t_fnmadd_pd( rtre, r1im, r04 );
                    #else
                    r03 = _t_mul_pd( rtre, r1re );
                    r04 = _t_mul_pd( rtim, r1re );
                    r01 = _t_mul_pd( rtim, r1im );
                    r02 = _t_mul_pd( rtre, r1im );
                    
                    r03 = _t_sub_pd( r0re, r03  );
                    r04 = _t_sub_pd( r0im, r04  );
                    
                    r03  = _t_add_pd( r03, r01 );
                    r04  = _t_sub_pd( r04, r02 );
                    #endif
                    
                    r0re = _t_add_pd( r0re, r0re );
                    r0im = _t_add_pd( r0im, r0im );
                    
                    _t_store_pd( px1re, r03 );
                    _t_store_pd( px1im, r04 );
                    
                    r0re = _t_sub_pd( r0re, r03 );
                    r0im = _t_sub_pd( r0im, r04 );
                    
                    _t_store_pd( px0re, r0re );
                    _t_store_pd( px0im, r0im );
                    
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
        
        // After i3 cycle, the address offset is 2*l2*step. Overall step before
        // next stage needed is 4*l2*step, therefore more walking.
        px0re += step2;
        px0im += step2;
        px1re += step2;
        px1im += step2;
        
    }
    
}