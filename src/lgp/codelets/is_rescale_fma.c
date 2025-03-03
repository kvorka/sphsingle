#include <stdlib.h>
#include <immintrin.h>
#include <stdio.h>

extern inline void is_rescale_c( const int n,
                                 const double *restrict amj,
                                 double *restrict cc ) {
  
  // rescale loop
  for ( int i1 = 0; i1 < n; i1++ ) {
    
    // rescale
    _mm256_store_pd( cc+4*i1, _mm256_mul_pd( _mm256_broadcast_sd( amj+i1 ), _mm256_load_pd( cc+4*i1 ) ) );
    
  }
  
}