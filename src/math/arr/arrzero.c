#include <stdlib.h>
#include <immintrin.h>

extern inline void zero_rarray_c( int n, double *restrict arr )

#if defined ( avx ) || defined( fma )
{
    
    // constant needed for loop unrolling
    const int n32 = (n/32)*32;
    
    // zero avx variable
    const __m256d rzero = _mm256_setzero_pd();
    
    // zero: cycle is unrolled by 32, factor of 4 is handled by an explicit
    // vectorization, factor of 8 is added in order to unroll the cycle a bit more for efficiency
    for ( int i2 = 0; i2 < n32; i2+=32 ) {
      for ( int i1 = 0; i1 < 32; i1+=4 ) {
        
        _mm256_store_pd( arr+i2+i1, rzero );
        
      }
    }
    
    // zero: remainder cases
    for ( int i1 = 0; i1 < n-n32; i1+=4 ) {
        
      _mm256_store_pd( arr+n32+i1, rzero );
      
    }
    
}
#endif