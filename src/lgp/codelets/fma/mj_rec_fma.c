#include <stdlib.h>
#include <immintrin.h>

extern inline void mj_rec_c( const int n,                    // howmany roots (step)
                            const double *restrict cff,     // recursion coeffs
                            const double *restrict cosx2,   // roots**2
                            double *restrict pmj1,          // Lege polys from previous step
                            double *restrict pmj ) {        // Lege polys
  
  // constant needed for loop unrolling
  const int n16 = (n/16)*16;
  
  // avx vars for recursion coeffs and Legendre polynomials
  const __m256d rcff1 = _mm256_broadcast_sd( cff   );
  const __m256d rcff2 = _mm256_broadcast_sd( cff+1 );
        __m256d rpmj;
  
  // recursion: cycle over the roots, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
  // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
  for ( int i2 = 0; i2 < n16; i2+=16 ) {
    for ( int i1 = 0; i1 < 16; i1+=4 ) {
      
      rpmj = _mm256_fmsub_pd( rcff1, _mm256_load_pd( cosx2+i1+i2 ), rcff2 );
      rpmj = _mm256_fmsub_pd( rpmj,  _mm256_load_pd( pmj1 +i1+i2 ), _mm256_load_pd( pmj+i1+i2 ) );
      
      _mm256_store_pd( pmj+i1+i2, rpmj );
      
    }
  }
  
  // recursion: remainder cases
  for ( int i1 = 0; i1 < 8; i1+=4 ) {
      
    rpmj = _mm256_fmsub_pd( rcff1, _mm256_load_pd( cosx2+i1+n16 ), rcff2 );
    rpmj = _mm256_fmsub_pd( rpmj,  _mm256_load_pd( pmj1 +i1+n16 ), _mm256_load_pd( pmj+i1+n16 ) );
    
    _mm256_store_pd( pmj+i1+n16, rpmj );
    
  }
  
}