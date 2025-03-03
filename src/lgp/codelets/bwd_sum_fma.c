#include <stdlib.h>
#include <immintrin.h>

extern inline void bwd_sum_c( const int n,                 // howmany roots (n)
                              const double *restrict pmj,  // Legendre polynomials
                              const double *restrict cc,   // sph coeffs
                              double *restrict swork) {    // partial sums
  
  __m256d rpmj;    // Legendre polynomials
  __m256d rcc[4];  // sph coeffs
  
  // load and broadcast sph coeffs
  for ( int j = 0; j < 4; j++ ) { rcc[j] = _mm256_broadcast_sd( cc+j ); }
  
  // bwd sum: cycle over the roots, the outer cycle is unrolled by 8, factor of 4 is handled by an explicit
  // vectorization, factor of 2 is added in order to unroll the cycle a bit more for efficiency
  for ( int i2 = 0; i2 < n; i2+=8 ) {
    for ( int j = 0; j < 2; j++ ) {
      
      // load Legendre polynomials
      rpmj = _mm256_load_pd( pmj+i2+4*j );
      
      // sum over roots and sph coeffs
      _mm256_store_pd( swork+    i2+4*j, _mm256_fmadd_pd( rpmj, rcc[0], _mm256_load_pd( swork+    i2+4*j ) ) );
      _mm256_store_pd( swork+  n+i2+4*j, _mm256_fmadd_pd( rpmj, rcc[1], _mm256_load_pd( swork+  n+i2+4*j ) ) );
      _mm256_store_pd( swork+2*n+i2+4*j, _mm256_fmadd_pd( rpmj, rcc[2], _mm256_load_pd( swork+2*n+i2+4*j ) ) );
      _mm256_store_pd( swork+3*n+i2+4*j, _mm256_fmadd_pd( rpmj, rcc[3], _mm256_load_pd( swork+3*n+i2+4*j ) ) );
      
    }
  }
  
}