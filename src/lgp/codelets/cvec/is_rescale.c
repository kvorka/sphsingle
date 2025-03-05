#include <stdlib.h>
#include <immintrin.h>

extern inline void is_rescale_c( const int n,
                                 const double *restrict amj,
                                 double *restrict cc )

#if defined( avx ) || defined( fma )
{
    
    // rescale loop: unroll by 4 for efficiency
    for ( int i1 = 0; i1 < (n/4)*4; i1+=4 ) {
      
      _mm256_store_pd( cc+4*i1   , _mm256_mul_pd( _mm256_broadcast_sd( amj+i1  ), _mm256_load_pd( cc+4*i1    ) ) );
      _mm256_store_pd( cc+4*i1+ 4, _mm256_mul_pd( _mm256_broadcast_sd( amj+i1+1), _mm256_load_pd( cc+4*i1+ 4 ) ) );
      _mm256_store_pd( cc+4*i1+ 8, _mm256_mul_pd( _mm256_broadcast_sd( amj+i1+2), _mm256_load_pd( cc+4*i1+ 8 ) ) );
      _mm256_store_pd( cc+4*i1+12, _mm256_mul_pd( _mm256_broadcast_sd( amj+i1+3), _mm256_load_pd( cc+4*i1+12 ) ) );
      
    }
    
    // remainder cases
    if ( n%4 != 0 ) {
      
      for ( int i1 = (n/4)*4; i1 < n; i1++) {
        _mm256_store_pd( cc+4*i1, _mm256_mul_pd( _mm256_broadcast_sd( amj+i1 ), _mm256_load_pd( cc+4*i1 ) ) );
      }
      
    }
    
}
#endif