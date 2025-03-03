#include <stdlib.h>
#include <immintrin.h>

extern inline void zero_rarray_c( const int n, double *restrict arr ) {
  
  const __m256d rzero = _mm256_setzero_pd();
  
  for ( int i1 = 0; i1 < n; i1+=4 ) { _mm256_store_pd( arr+i1, rzero ); }
  
}