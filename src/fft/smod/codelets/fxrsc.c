#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void fxrsc_c( const int m,
              const double fac,
                    double *restrict arr )

{
    
    // Memory references
    double *restrict parr = arr;
    
    // Registers to be used
    const __td rfac = _t_set1_pd( fac );
          __td r0, r1, r2, r3;
    
    // Main cycle
    for ( int i = 0; i < m; i++ ) {
        
        r0 = _t_load_pd( parr + vlen0 );
        r1 = _t_load_pd( parr + vlen1 );
        
        r0 = _t_mul_pd( rfac, r0 );
        r1 = _t_mul_pd( rfac, r1 );
        
        r2 = _t_load_pd( parr + vlen2 );
        r3 = _t_load_pd( parr + vlen3 );
        
        _t_store_pd( parr + vlen0, r0 );
        _t_store_pd( parr + vlen1, r1 );
        
        r2 = _t_mul_pd( rfac, r2 );
        r3 = _t_mul_pd( rfac, r3 );
        
        _t_store_pd( parr + vlen2, r2 );
        _t_store_pd( parr + vlen3, r3 );
        
        parr += vlen4;
        
    }
    
}