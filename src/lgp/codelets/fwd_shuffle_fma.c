#include <stdlib.h>
#include <immintrin.h>

extern inline void fwd_shuffle_c( const int n,                  // howmany roots (step)
                                  const double *restrict cosx,  // roots
                                  const double *restrict wght,  // weights
                                  double *restrict sumN,        // North sum
                                  double *restrict sumS,        // South sum
                                  double *restrict swork ) {    // partial sum
  
  // avx vars for cosine and weight values, North and South sums
  __m256d rwcsx, rwght, rsumN[2], rsumS[2];
  
  // ??
  for ( int i1 = 0; i1 < n; i1+=8 ) {
    for ( int j = 0; j < 2; j++) {
      
      rwght =                _mm256_load_pd( wght+i1+4*j );
      rwcsx = _mm256_mul_pd( _mm256_load_pd( cosx+i1+4*j ), rwght );
      
      rsumN[0] = _mm256_load_pd( sumN+i1+4*j   );
      rsumN[1] = _mm256_load_pd( sumN+i1+4*j+n );
      rsumS[0] = _mm256_load_pd( sumS+i1+4*j   );
      rsumS[1] = _mm256_load_pd( sumS+i1+4*j+n );
      
      _mm256_store_pd( swork+i1+4*j    , _mm256_mul_pd( _mm256_sub_pd( rsumN[0], rsumS[0] ), rwght ) );
      _mm256_store_pd( swork+i1+4*j+  n, _mm256_mul_pd( _mm256_sub_pd( rsumN[1], rsumS[1] ), rwght ) );
      _mm256_store_pd( swork+i1+4*j+2*n, _mm256_mul_pd( _mm256_add_pd( rsumN[0], rsumS[0] ), rwcsx ) );
      _mm256_store_pd( swork+i1+4*j+3*n, _mm256_mul_pd( _mm256_add_pd( rsumN[1], rsumS[1] ), rwcsx ) );
      
    }
  }
  
}
