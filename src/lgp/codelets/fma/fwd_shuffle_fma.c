#include <stdlib.h>
#include <immintrin.h>

extern inline void fwd_shuffle_c( const int n,                  // howmany roots (step)
                                  const double *restrict cosx,  // roots
                                  const double *restrict wght,  // weights
                                  double *restrict grid,        // Legendre sums
                                  double *restrict swork ) {    // partial sum
  
  // constant needed for loop unrolling
  const int n32 = (n/32)*32;
  
  // avx vars for cosine and weight values, North and South sums
  __m256d rwcsx, rwght, rsumN[2], rsumS[2];
  
  // fwd shuffle: cycle over the roots, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
  // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
  for ( int i1 = 0; i1 < n32; i1+=32 ) {
    for ( int i2 = 0; i2 < 32; i2+=4 ) {
      
      rwght =                _mm256_load_pd( wght+i1+i2 );
      rwcsx = _mm256_mul_pd( _mm256_load_pd( cosx+i1+i2 ), rwght );
      
      rsumN[0] = _mm256_load_pd( grid+i1+i2     );
      rsumS[0] = _mm256_load_pd( grid+i1+i2+  n );
      rsumN[1] = _mm256_load_pd( grid+i1+i2+2*n );
      rsumS[1] = _mm256_load_pd( grid+i1+i2+3*n );
      
      _mm256_store_pd( swork+i1+i2    , _mm256_mul_pd( _mm256_sub_pd( rsumN[0], rsumS[0] ), rwght ) );
      _mm256_store_pd( swork+i1+i2+  n, _mm256_mul_pd( _mm256_sub_pd( rsumN[1], rsumS[1] ), rwght ) );
      _mm256_store_pd( swork+i1+i2+2*n, _mm256_mul_pd( _mm256_add_pd( rsumN[0], rsumS[0] ), rwcsx ) );
      _mm256_store_pd( swork+i1+i2+3*n, _mm256_mul_pd( _mm256_add_pd( rsumN[1], rsumS[1] ), rwcsx ) );
      
    }
  }
  
  // fwd shuffle: remainder cases
  for ( int i2 = 0; i2 < n-n32; i2+=4 ) {
      
    rwght =                _mm256_load_pd( wght+n32+i2 );
    rwcsx = _mm256_mul_pd( _mm256_load_pd( cosx+n32+i2 ), rwght );
    
    rsumN[0] = _mm256_load_pd( grid+n32+i2     );
    rsumS[0] = _mm256_load_pd( grid+n32+i2+  n );
    rsumN[1] = _mm256_load_pd( grid+n32+i2+2*n );
    rsumS[1] = _mm256_load_pd( grid+n32+i2+3*n );
    
    _mm256_store_pd( swork+n32+i2    , _mm256_mul_pd( _mm256_sub_pd( rsumN[0], rsumS[0] ), rwght ) );
    _mm256_store_pd( swork+n32+i2+  n, _mm256_mul_pd( _mm256_sub_pd( rsumN[1], rsumS[1] ), rwght ) );
    _mm256_store_pd( swork+n32+i2+2*n, _mm256_mul_pd( _mm256_add_pd( rsumN[0], rsumS[0] ), rwcsx ) );
    _mm256_store_pd( swork+n32+i2+3*n, _mm256_mul_pd( _mm256_add_pd( rsumN[1], rsumS[1] ), rwcsx ) );
    
  }
  
}
