#include "../../cvec.h"

extern inline __attribute__((always_inline))
void zero_rarray_c( const int length,
                          double *restrict arr )

{
    
    // Memory address
    double *restrict parr = arr;
    
    // Registers to be used
    const __td r00 = _t_setzero_pd();
    
    // Iterator
    int i = 0;
    
    // Main loop
    for ( ; i <= length - vlen4; i += vlen4 ) {
        
        _t_store_pd( parr + vlen0, r00 );
        _t_store_pd( parr + vlen1, r00 );
        _t_store_pd( parr + vlen2, r00 );
        _t_store_pd( parr + vlen3, r00 );
        
        parr += vlen4;
        
    }
    
    // Loop remainder
    for ( ; i <= length - vlen; i += vlen ) { 
        
        _t_store_pd( parr, r00 );
        
        parr += vlen;
        
    }
    
}