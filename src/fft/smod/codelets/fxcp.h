#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fxcpy_c( const int m,
              const double *restrict arr_from,
                    double *restrict arr_to )

{
    
    #pragma omp unroll (vlen8) simd uniform (fac) aligned (arr_to,arr_from:alignement)
    for ( int i = 0; i < vlen8 * m; i++ ) { arr_to[i] = arr_from[i]; }
    
}