#include <stdlib.h>
#include <immintrin.h>

extern inline void mj_rec_c( const int n,                  // howmany roots (step)
                            const double *restrict cff,    // recursion coeffs
                            const double *restrict cosx2,  // roots**2
                            double *restrict pmj1,         // Lege polys from previous step
                            double *restrict pmj ) {       // Lege polys
  
  // constant needed for loop unrolling
  const int n32 = (n/32)*32;
  
  // avx vars for recursion coeffs and Legendre polynomials
  const __m256d rcff1 = _mm256_broadcast_sd( cff   );
  const __m256d rcff2 = _mm256_broadcast_sd( cff+1 );
        __m256d rpmj;
  
  // recursion: cycle over the roots, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
  // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
  for ( int i2 = 0; i2 < n32; i2+=32 ) {
    for ( int i1 = 0; i1 < 32; i1+=4 ) {
      
      rpmj = _mm256_sub_pd( _mm256_mul_pd( rcff1, _mm256_load_pd( cosx2+i2+i1 ) ), rcff2 );
      rpmj = _mm256_sub_pd( _mm256_mul_pd( rpmj,  _mm256_load_pd( pmj1 +i2+i1 ) ), _mm256_load_pd( pmj+i2+i1 ) );
      
      _mm256_store_pd( pmj+i2+i1, rpmj );
      
    }
  }
  
  // recursion: remainer cases
  for ( int i1 = 0; i1 < n-n32; i1+=4 ) {
      
    rpmj = _mm256_sub_pd( _mm256_mul_pd( rcff1, _mm256_load_pd( cosx2+n32+i1 ) ), rcff2 );
    rpmj = _mm256_sub_pd( _mm256_mul_pd( rpmj,  _mm256_load_pd( pmj1 +n32+i1 ) ), _mm256_load_pd( pmj+n32+i1 ) );
    
    _mm256_store_pd( pmj+n32+i1, rpmj );
    
  }
  
}