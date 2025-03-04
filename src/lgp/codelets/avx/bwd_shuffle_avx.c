#include <stdlib.h>
#include <immintrin.h>

extern inline void bwd_shuffle_c( const int n,                   // howmany roots (step)
                                  const double *restrict cosx,   // roots
                                  const double *restrict swork,  // partial sums to be shuffled
                                  double *restrict grid) {       // Legendre polys sum
  
  // constant needed for loop unrolling
  const int n32 = (n/32)*32;
  
  // avx vars for cosine values and partial sums
  __m256d rcosx, rssym[2], rasym[2];
  
  // bwd shuffle: cycle over the roots, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
  // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
  for ( int i2 = 0; i2 < n32; i2+=32 ) {
    for ( int i1 = 0; i1 < 32; i1+=4 ) {
      
      rcosx = _mm256_load_pd( cosx+i1+i2 );
      
      rasym[0] =                       _mm256_load_pd( swork+i2+i1     );
      rasym[1] =                       _mm256_load_pd( swork+i2+i1+  n );
      rssym[0] = _mm256_mul_pd( rcosx, _mm256_load_pd( swork+i2+i1+2*n ) );
      rssym[1] = _mm256_mul_pd( rcosx, _mm256_load_pd( swork+i2+i1+3*n ) );
      
      _mm256_store_pd( grid+i2+i1    , _mm256_add_pd( rssym[0], rasym[0] ) );
      _mm256_store_pd( grid+i2+i1+  n, _mm256_sub_pd( rssym[0], rasym[0] ) );
      _mm256_store_pd( grid+i2+i1+2*n, _mm256_add_pd( rssym[1], rasym[1] ) );
      _mm256_store_pd( grid+i2+i1+3*n, _mm256_sub_pd( rssym[1], rasym[1] ) );
      
    }
  }
  
  // bwd shuffle: remainder cases
  for ( int i1 = 0; i1 < n-n32; i1+=4) {
    
    rcosx = _mm256_load_pd( cosx+i1+n32 );
    
    rasym[0] =                       _mm256_load_pd( swork+i1+n32     );
    rasym[1] =                       _mm256_load_pd( swork+i1+n32+  n );
    rssym[0] = _mm256_mul_pd( rcosx, _mm256_load_pd( swork+i1+n32+2*n ) );
    rssym[1] = _mm256_mul_pd( rcosx, _mm256_load_pd( swork+i1+n32+3*n ) );
    
    _mm256_store_pd( grid+i1+n32    , _mm256_add_pd( rssym[0], rasym[0] ) );
    _mm256_store_pd( grid+i1+n32+  n, _mm256_sub_pd( rssym[0], rasym[0] ) );
    _mm256_store_pd( grid+i1+n32+2*n, _mm256_add_pd( rssym[1], rasym[1] ) );
    _mm256_store_pd( grid+i1+n32+3*n, _mm256_sub_pd( rssym[1], rasym[1] ) );
    
  }
  
}
