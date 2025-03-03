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
  
  for ( int i1 = 0; i1 < n; i1+=8 ) {
    for ( int j = 0; j < 2; j++) {
    
      _mm256_store_pd( pmj1+i1+4*j, rcff );
      
    }
  }
  
  // set pmm, pmj: n1 >= 16 by design, the loop is unrolled by 4
  rcff = _mm256_set1_pd( cff );
  
  switch ( ma ) {
    
    case 1:
    {
      
      for ( int i1 = 0; i1 < n; i1+=8 ) {
        for ( int j = 0; j < 2; j++ ) {
          
          _mm256_store_pd( pmm+i1+4*j, rcff );
          _mm256_store_pd( pmj+i1+4*j, _mm256_div_pd( rcff, _mm256_load_pd( cosx+i1+4*j ) ) );
          
        }
      }
      
    }
    break;
    
    default:
    {
      
      for ( int i1 = 0; i1 < n; i1+=8 ) {
        for ( int j = 0; j < 2; j++ ) {
          
          rpmm = _mm256_mul_pd( rcff, _mm256_load_pd( sinx+i1+4*j ) );
          rpmm = _mm256_mul_pd( rpmm, _mm256_load_pd( pmm +i1+4*j ) );
          
          _mm256_store_pd( pmm+i1+4*j, rpmm );
          _mm256_store_pd( pmj+i1+4*j, _mm256_div_pd( rpmm, _mm256_load_pd( cosx+i1+4*j ) ) );
          
        }
      }
      
    }
    break;
    
  }
  
}
