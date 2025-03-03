#include <stdlib.h>
#include <immintrin.h>

extern inline void fwd_sum_c( const int n,                     // howmany roots (n)
                              const double *restrict pmj,      // Legendre polynomials
                              const double *restrict swork,    // partial sums
                              double *restrict cc) {           // sph coeffs
  
  __m256d rpmj;    // Legendre polynomials 
  __m256d rcc[4];  // sph coeffs accumulators
  __m128d rsum;    // temporary sum of sph coeffs
  
   // set accumulators to zero (suboptimal)
   for ( int j = 0; j < 4; j++ ) { rcc[j] = _mm256_setzero_pd(); }
  
  // sums: the inner loop is unrolled by 4 with the same reasoning and strategy as above,
  // the outer loop is unrolled by 4, because n2 is guaranteed to be a multiple of 4
  for ( int i2 = 0; i2 < n; i2+=8 ) {
    for ( int j = 0; j < 2; j++ ) {
      
      // load Legendre polynomials
      rpmj = _mm256_load_pd( pmj+i2+4*j );
      
      // sum over roots and sph coeffs
      rcc[0] = _mm256_add_pd( _mm256_mul_pd( rpmj, _mm256_load_pd( swork+    i2+4*j ) ), rcc[0] );
      rcc[1] = _mm256_add_pd( _mm256_mul_pd( rpmj, _mm256_load_pd( swork+  n+i2+4*j ) ), rcc[1] );
      rcc[2] = _mm256_add_pd( _mm256_mul_pd( rpmj, _mm256_load_pd( swork+2*n+i2+4*j ) ), rcc[2] );
      rcc[3] = _mm256_add_pd( _mm256_mul_pd( rpmj, _mm256_load_pd( swork+3*n+i2+4*j ) ), rcc[3] );
      
    }
  }
  
  // horizontal sum of the accumulators, store results in the sph coeffs
  for ( int j = 0; j < 4; j++ ) {
    
    rsum   = _mm_add_pd( _mm256_extractf128_pd( rcc[j], 0 ), _mm256_extractf128_pd( rcc[j], 1 ) );
    cc[j] += _mm_cvtsd_f64( _mm_add_pd( rsum, _mm_unpackhi_pd( rsum, rsum ) ) );
    
  }
  
}