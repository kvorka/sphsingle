#include <stdlib.h>
#include <immintrin.h>

extern inline void bwd_sum_c( const int n,                 // howmany roots (n)
                              const double *restrict pmj,  // Legendre polynomials
                              const double *restrict cc,   // sph coeffs
                              double *restrict swork)      // partial sums
  
#if defined( avx )
{
        
    // constant needed for loop unrolling
    const int n32 = (n/32)*32;
    
    // avx vars for Legendre polynomials and sph coeffs
    __m256d rpmj, rcc[4];
    
    // load and broadcast sph coeffs
    for ( int j = 0; j < 4; j++ ) { rcc[j] = _mm256_broadcast_sd( cc+j ); }
    
    // bwd sum: cycle over the roots, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
    // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
    for ( int i2 = 0; i2 < n32; i2+=32 ) {
      for ( int i1 = 0; i1 < 32; i1+=4 ) {
        
        // load Legendre polynomials
        rpmj = _mm256_load_pd( pmj+i2+i1 );
        
        // sum over roots and sph coeffs
        _mm256_store_pd( swork+    i2+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[0] ), _mm256_load_pd( swork+    i2+i1 ) ) );
        _mm256_store_pd( swork+  n+i2+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[1] ), _mm256_load_pd( swork+  n+i2+i1 ) ) );
        _mm256_store_pd( swork+2*n+i2+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[2] ), _mm256_load_pd( swork+2*n+i2+i1 ) ) );
        _mm256_store_pd( swork+3*n+i2+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[3] ), _mm256_load_pd( swork+3*n+i2+i1 ) ) );
        
      }
    }
    
    // bwd sum: remainer cases
    for ( int i1 = 0; i1 < n-n32; i1+=4) {
      
      // load Legendre polynomials
      rpmj = _mm256_load_pd( pmj+n32+i1 );
        
      // sum over roots and sph coeffs
      _mm256_store_pd( swork+    n32+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[0] ), _mm256_load_pd( swork+    n32+i1 ) ) );
      _mm256_store_pd( swork+  n+n32+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[1] ), _mm256_load_pd( swork+  n+n32+i1 ) ) );
      _mm256_store_pd( swork+2*n+n32+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[2] ), _mm256_load_pd( swork+2*n+n32+i1 ) ) );
      _mm256_store_pd( swork+3*n+n32+i1, _mm256_add_pd( _mm256_mul_pd( rpmj, rcc[3] ), _mm256_load_pd( swork+3*n+n32+i1 ) ) );
      
    }
    
}
#elif defined ( fma )
{    
    
    // constant needed for loop unrolling
    const int n32 = (n/32)*32;
    
    // avx vars for Legendre polynomials and sph coeffs
    __m256d rpmj, rcc[4];
    
    // load and broadcast sph coeffs
    for ( int j = 0; j < 4; j++ ) { rcc[j] = _mm256_broadcast_sd( cc+j ); }
    
    // bwd sum: cycle over the roots, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
    // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
    for ( int i2 = 0; i2 < n32; i2+=32 ) {
      for ( int i1 = 0; i1 < 32; i1+=4 ) {
        
        // load Legendre polynomials
        rpmj = _mm256_load_pd( pmj+i2+i1 );
        
        // sum over roots and sph coeffs
        _mm256_store_pd( swork+    i2+i1, _mm256_fmadd_pd( rpmj, rcc[0], _mm256_load_pd( swork+    i2+i1 ) ) );
        _mm256_store_pd( swork+  n+i2+i1, _mm256_fmadd_pd( rpmj, rcc[1], _mm256_load_pd( swork+  n+i2+i1 ) ) );
        _mm256_store_pd( swork+2*n+i2+i1, _mm256_fmadd_pd( rpmj, rcc[2], _mm256_load_pd( swork+2*n+i2+i1 ) ) );
        _mm256_store_pd( swork+3*n+i2+i1, _mm256_fmadd_pd( rpmj, rcc[3], _mm256_load_pd( swork+3*n+i2+i1 ) ) );
        
      }
    }
    
    // bwd sum: remainder cases
    for ( int i1 = 0; i1 < n-n32; i1+=4 ) {
        
      // load Legendre polynomials
      rpmj = _mm256_load_pd( pmj+n32+i1 );
      
      // sum over roots and sph coeffs
      _mm256_store_pd( swork+    n32+i1, _mm256_fmadd_pd( rpmj, rcc[0], _mm256_load_pd( swork+    n32+i1 ) ) );
      _mm256_store_pd( swork+  n+n32+i1, _mm256_fmadd_pd( rpmj, rcc[1], _mm256_load_pd( swork+  n+n32+i1 ) ) );
      _mm256_store_pd( swork+2*n+n32+i1, _mm256_fmadd_pd( rpmj, rcc[2], _mm256_load_pd( swork+2*n+n32+i1 ) ) );
      _mm256_store_pd( swork+3*n+n32+i1, _mm256_fmadd_pd( rpmj, rcc[3], _mm256_load_pd( swork+3*n+n32+i1 ) ) );
      
    }
    
}
#endif