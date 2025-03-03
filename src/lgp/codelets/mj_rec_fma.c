#include <stdlib.h>
#include <immintrin.h>

extern inline void mj_rec_c( const int n,                    // howmany roots (step)
                            const double *restrict cff,     // recursion coeffs
                            const double *restrict cosx2,   // roots**2
                            double *restrict pmj1,          // Lege polys from previous step
                            double *restrict pmj ) {        // Lege polys
  
  // avx variables for coefficients and polynomials
  const __m256d rcff1 = _mm256_broadcast_sd( cff   );
  const __m256d rcff2 = _mm256_broadcast_sd( cff+1 );
  __m256d rpmj;
  
  // recursion
  for ( int i1 = 0; i1 < n; i1+=8 ) {
    for ( int j = 0; j < 2; j++ ) {
      
      rpmj = _mm256_fmsub_pd( rcff1, _mm256_load_pd( cosx2+i1+4*j ), rcff2 );
      rpmj = _mm256_fmsub_pd( rpmj,  _mm256_load_pd( pmj1 +i1+4*j ), _mm256_load_pd( pmj+i1+4*j ) );
      
      _mm256_store_pd( pmj+i1+4*j, rpmj );
      
    }
  }
  
}