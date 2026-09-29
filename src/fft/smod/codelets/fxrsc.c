#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void fxrsc_c( const int m,
              const double fac,
                    double *restrict arr )

{
    
    #pragma omp unroll (vlen4) simd uniform (fac) aligned (arr:alignement)
    for ( int i = 0; i < vlen4 * m; i++ ) { arr[i] = fac * arr[i]; }
    
}