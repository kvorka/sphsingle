#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxcpy_c( const int m,
              const double *restrict arr_from,
                    double *restrict arr_to )

{
    
    // Memory references
    const double *restrict pf = arr_from;
          double *restrict pt = arr_to;
    
    // Registers to be used
    __td r0, r1, r2, r3,
         r4, r5, r6, r7;
    
    // Main cycle
    for ( int i = 0; i < m; i++ ) { 
        
        r0 = _t_load_pd( pf + vlen0 );
        r1 = _t_load_pd( pf + vlen1 );
        
        _t_store_pd( pt + vlen0, r0 );
        _t_store_pd( pt + vlen1, r1 );
        
        r2 = _t_load_pd( pf + vlen2 );
        r3 = _t_load_pd( pf + vlen3 );
        
        _t_store_pd( pt + vlen2, r2 );
        _t_store_pd( pt + vlen3, r3 );
        
        r4 = _t_load_pd( pf + vlen4 );
        r5 = _t_load_pd( pf + vlen5 );
        
        _t_store_pd( pt + vlen4, r4 );
        _t_store_pd( pt + vlen5, r5 );
        
        r6 = _t_load_pd( pf + vlen6 );
        r7 = _t_load_pd( pf + vlen7 );
        
        _t_store_pd( pt + vlen6, r6 );
        _t_store_pd( pt + vlen7, r7 );
        
        pt += vlen8;
        pf += vlen8;
        
    }
    
}