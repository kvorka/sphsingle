#include <stdlib.h>
#include <immintrin.h>

extern inline void bwd_shuffle_c( const int n,                   // howmany roots (step)
                                  const double *restrict cosx,   // roots
                                  const double *restrict swork,  // partial sums to be shuffled
                                  double *restrict sumN,         // North Legendre polys sum
                                  double *restrict sumS ) {      // South Legendre polys sum
  
  // avx vars for cosine values and partial sums
  __m256d rcosx, rssym[2], rasym[2];
  
  // ??
  for ( int i1 = 0; i1 < n; i1+=8 ) {
    for ( int j = 0; j < 2; j++ ) {
      
      rcosx = _mm256_load_pd( cosx+i1+4*j );
      
      rasym[0] = _mm256_load_pd( swork+i1+4*j     );
      rasym[1] = _mm256_load_pd( swork+i1+4*j+  n );
      rssym[0] = _mm256_load_pd( swork+i1+4*j+2*n );
      rssym[1] = _mm256_load_pd( swork+i1+4*j+3*n );
      
      _mm256_store_pd( sumN+i1+4*j  , _mm256_fmadd_pd( rcosx, rssym[0], rasym[0] ) );
      _mm256_store_pd( sumS+i1+4*j  , _mm256_fmsub_pd( rcosx, rssym[0], rasym[0] ) );
      _mm256_store_pd( sumN+i1+4*j+n, _mm256_fmadd_pd( rcosx, rssym[1], rasym[1] ) );
      _mm256_store_pd( sumS+i1+4*j+n, _mm256_fmsub_pd( rcosx, rssym[1], rasym[1] ) );
      
    }
  }
  
}
