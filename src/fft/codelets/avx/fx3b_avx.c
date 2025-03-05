#include <stdlib.h>
#include <immintrin.h>

extern inline void fxzm3b_c( const int n,
                             const int l,
                             double *restrict x ) {
  
  const __m256d rc31 = _mm256_set1_pd( 0.50000000000000000000 );
  const __m256d rc32 = _mm256_set1_pd( 0.86602540378443864676 );
  
  __m256d rxre[6], rxim[6];
  
  for ( int i = 0; i < l; i++) {
    for ( int iv = 0; iv < n; iv+=16 ) {
      for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
        
        rxre[4] = _mm256_load_pd( x+iv1+iv+  2*i*n+2*l*n );
        rxim[4] = _mm256_load_pd( x+iv1+iv+n+2*i*n+2*l*n );
        rxre[5] = _mm256_load_pd( x+iv1+iv+  2*i*n+4*l*n );
        rxim[5] = _mm256_load_pd( x+iv1+iv+n+2*i*n+4*l*n );
        
        rxre[0] = _mm256_add_pd( rxre[4], rxre[5] );
        rxim[0] = _mm256_add_pd( rxim[4], rxim[5] );
        rxre[1] = _mm256_sub_pd( rxre[4], rxre[5] );
        rxim[1] = _mm256_sub_pd( rxim[4], rxim[5] );
        rxre[2] = _mm256_mul_pd( rc31, rxre[0] );
        rxim[2] = _mm256_mul_pd( rc31, rxim[0] );
        
        rxre[3] = _mm256_load_pd( x+iv1+iv+  2*i*n );
        rxim[3] = _mm256_load_pd( x+iv1+iv+n+2*i*n );
        
        rxre[0] = _mm256_add_pd( rxre[0], rxre[3] );
        rxim[0] = _mm256_add_pd( rxim[0], rxim[3] );
        rxre[2] = _mm256_add_pd( rxre[2], rxre[3] );
        rxim[2] = _mm256_add_pd( rxim[2], rxim[3] );
        
        _mm256_store_pd( x+iv1+iv+  2*i*n, rxre[0] );
        _mm256_store_pd( x+iv1+iv+n+2*i*n, rxim[0] );
        
        rxre[1] = _mm256_mul_pd( rc32, rxre[1] );
        rxim[1] = _mm256_mul_pd( rc32, rxim[1] );
        
        rxre[0] = _mm256_add_pd( rxre[2], rxim[1] );
        rxim[0] = _mm256_sub_pd( rxim[2], rxre[1] );
        
        _mm256_store_pd( x+iv1+iv+  2*i*n+4*l*n, rxre[0] );
        _mm256_store_pd( x+iv1+iv+n+2*i*n+4*l*n, rxim[0] );
        
        rxre[2] = _mm256_add_pd( rxre[2], rxre[2] );
        rxim[2] = _mm256_add_pd( rxim[2], rxim[2] );
        
        rxre[2] = _mm256_sub_pd( rxre[2], rxre[5] );
        rxim[2] = _mm256_sub_pd( rxim[2], rxim[5] );
        
        _mm256_store_pd( x+iv1+iv+  2*i*n+2*l*n, rxre[1] );
        _mm256_store_pd( x+iv1+iv+n+2*i*n+2*l*n, rxim[1] );
        
      }
    }
  }
  
}