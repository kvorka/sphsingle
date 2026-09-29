#include "../cvec.h"

extern inline __attribute__((always_inline))
void zero_rarray_c( const int length,
                          double *restrict arr )

{
    
    // Main loop
    #pragma omp unroll partial (vlen4) simd
    for ( int i = 0; i < length; i++ ) { arr[i] = 0.; }
    
}