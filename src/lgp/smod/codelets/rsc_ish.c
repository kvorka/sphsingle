#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void is_rescale_c( const int n,
                   const double *restrict amj,
                         double *restrict rcab )

{
    
    // Registers to be used
    __td rc, ra;
    
    // Main cycle
    #pragma omp unroll partial (6)
    for ( int i = 0; i < n; i++ ) {
        
        ra = _t_set1_pd( *( amj + i ) );
        rc = _t_loadu_pd( rcab + 4*i );
        
        rc = _t_mul_pd( ra, rc );
        
        _t_storeu_pd( rcab + 4*i, rc );
        
    }
    
}