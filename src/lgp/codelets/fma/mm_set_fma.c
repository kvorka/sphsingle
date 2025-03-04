#include <stdlib.h>
#include <immintrin.h>

extern inline void mm_set_c( const int ma,                  // identifier for m=0 case
                             const int n,                   // howmany roots (step)
                             const double cff,              // recursion coeffs
                             const double *restrict cosx,   // roots
                             const double *restrict sinx,   // sqrt(1-roots**2)
                             double *restrict pmm,          // Lege polys m=0
                             double *restrict pmj1,         // Lege polys for previous step holder
                             double *restrict pmj ) {       // Lege polys
  
  // avx variables for coefficients and polynomials
  __m256d rcff, rpmm;
  
  // set pmj1
  rcff = _mm256_setzero_pd();
  
  for ( int i2 = 0; i2 < n; i2+=8 ) {
    for ( int i1 = 0; i1 < 8; i1+=4 ) {
    
      _mm256_store_pd( pmj1+i2+i1, rcff );
      
    }
  }
  
  // set pmm, pmj
  rcff = _mm256_set1_pd( cff );
  
  switch ( ma ) {
    
    case 1:
    {
      
      for ( int i2 = 0; i2 < n; i2+=8 ) {
        for ( int i1 = 0; i1 < 8; i1+=4 ) {
          
          _mm256_store_pd( pmm+i2+i1, rcff );
          _mm256_store_pd( pmj+i2+i1, _mm256_div_pd( rcff, _mm256_load_pd( cosx+i2+i1 ) ) );
          
        }
      }
      
    }
    break;
    
    default:
    {
      
      for ( int i2 = 0; i2 < n; i2+=8 ) {
        for ( int i1 = 0; i1 < 8; i1+=4 ) {
          
          rpmm = _mm256_mul_pd( rcff, _mm256_load_pd( sinx+i2+i1 ) );
          rpmm = _mm256_mul_pd( rpmm, _mm256_load_pd( pmm +i2+i1 ) );
          
          _mm256_store_pd( pmm+i2+i1, rpmm );
          _mm256_store_pd( pmj+i2+i1, _mm256_div_pd( rpmm, _mm256_load_pd( cosx+i2+i1 ) ) );
          
        }
      }
      
    }
    break;
    
  }
  
}
