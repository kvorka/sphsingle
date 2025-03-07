#include <stdlib.h>
#include <immintrin.h>

extern inline void copy_rarray_c( int n, 
                                  const double *restrict arrfrom, 
                                  double *restrict arrto )

#if defined( avx ) || defined( fma )
{
    
    // constant needed for loop unrolling
    const int n32 = (n/32)*32;
    
    // copy: cycle over arrays, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
    // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
    for ( int i2 = 0; i2 < n32; i2+=32 ) {
      for ( int i1 = 0; i1 < 32; i1+=4 ) {
        
        _mm256_store_pd( arrto+i1+i2, _mm256_load_pd( arrfrom+i1+i2 ) );
        
      }
    }
    
    // copy: remainder cases
    for ( int i1 = 0; i1 < n-n32; i1+=4) {
      
      _mm256_store_pd( arrto+i1+n32, _mm256_load_pd( arrfrom+i1+n32 ) );
      
    }
    
}
#endif