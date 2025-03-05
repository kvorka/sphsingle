#include <stdlib.h>
#include <immintrin.h>

extern inline void fxzm3a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t)

#if defined( avx )
{
    
    const __m256d rc31 = _mm256_set1_pd( -0.50000000000000000000 );
    const __m256d rc32 = _mm256_set1_pd( +0.86602540378443864676 );
    
    __m256d rtre[2], rtim[2], rxre[5], rxim[5];
    
    for ( int j = 0; j < k; j++) {
      
      rtre[0] = _mm256_broadcast_sd( t+4*j   );
      rtim[0] = _mm256_broadcast_sd( t+4*j+1 );
      rtre[1] = _mm256_broadcast_sd( t+4*j+2 );
      rtim[1] = _mm256_broadcast_sd( t+4*j+3 );
      
      for ( int i = 0; i < l; i++) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
            
            rxre[3] = _mm256_load_pd( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n );
            rxim[3] = _mm256_load_pd( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n );
            rxre[4] = _mm256_load_pd( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n );
            rxim[4] = _mm256_load_pd( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n );
            
            rxre[0] = _mm256_mul_pd( rtre[0], rxre[3] );
            rxim[0] = _mm256_mul_pd( rtre[0], rxim[3] );
            rxre[2] = _mm256_mul_pd( rtim[0], rxim[3] );
            rxim[2] = _mm256_mul_pd( rtim[0], rxre[3] );
            
            rxre[1] = _mm256_mul_pd( rtre[1], rxre[4] );
            rxim[1] = _mm256_mul_pd( rtre[1], rxim[4] );
            rxre[3] = _mm256_mul_pd( rtim[1], rxim[4] );
            rxim[3] = _mm256_mul_pd( rtim[1], rxre[4] );
            
            rxre[0] = _mm256_sub_pd( rxre[0], rxre[2] );
            rxim[0] = _mm256_add_pd( rxim[0], rxim[2] );
            
            rxre[1] = _mm256_sub_pd( rxre[3], rxre[1] );
            rxim[1] = _mm256_add_pd( rxim[3], rxim[1] );
            
            rxre[1] = _mm256_add_pd( rxre[0], rxre[1] );
            rxim[1] = _mm256_sub_pd( rxim[0], rxim[1] );
            
            rxre[0] = _mm256_add_pd( rxre[0], rxre[0] );
            rxim[0] = _mm256_add_pd( rxim[0], rxim[0] );
            
            rxre[0] = _mm256_sub_pd( rxre[0], rxre[1] );
            rxim[0] = _mm256_sub_pd( rxim[0], rxim[1] );
            
            rxre[2] = _mm256_mul_pd( rc31, rxre[0] );
            rxim[2] = _mm256_mul_pd( rc31, rxim[0] );
            rxre[3] = _mm256_load_pd( x+iv1+iv  +2*i*n+6*j*l*n );
            rxim[3] = _mm256_load_pd( x+iv1+iv+n+2*i*n+6*j*l*n );
            
            rxre[0] = _mm256_add_pd( rxre[0], rxre[3] );
            rxim[0] = _mm256_add_pd( rxim[0], rxim[3] );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+6*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+6*j*l*n, rxim[0] );
            
            rxre[1] = _mm256_mul_pd( rc32, rxre[1] );
            rxim[1] = _mm256_mul_pd( rc32, rxim[1] );
            rxre[2] = _mm256_add_pd( rxre[2], rxre[3] );
            rxim[2] = _mm256_add_pd( rxim[2], rxim[3] );
            
            rxre[0] = _mm256_add_pd( rxre[2], rxim[1] );
            rxim[0] = _mm256_sub_pd( rxim[2], rxre[1] );
            rxre[4] = _mm256_add_pd( rxre[2], rxre[2] );
            rxim[4] = _mm256_add_pd( rxim[2], rxim[2] );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n, rxim[0] );
            
            rxre[4] = _mm256_sub_pd( rxre[4], rxre[0] );
            rxim[4] = _mm256_sub_pd( rxim[4], rxim[0] );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n, rxre[4] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n, rxim[4] );
            
          }
        }
      }
      
    }
    
}
#elif defined( fma )
{
    
    const __m256d rtwo = _mm256_set1_pd( +2.0 );
    const __m256d rm00 = _mm256_set1_pd( -0.0 );
    const __m256d rc31 = _mm256_set1_pd( -0.50000000000000000000 );
    const __m256d rc32 = _mm256_set1_pd( +0.86602540378443864676 );
    
    __m256d rtre[2], rtim[2], rxre[5], rxim[5];
    
    for ( int j = 0; j < k; j++) {
      
      rtre[0] = _mm256_broadcast_sd( t+4*j   );
      rtim[0] = _mm256_broadcast_sd( t+4*j+1 );
      rtre[1] = _mm256_broadcast_sd( t+4*j+2 );
      rtim[1] = _mm256_broadcast_sd( t+4*j+3 );
      
      for ( int i = 0; i < l; i++) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
            
            rxre[3] = _mm256_load_pd( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n );
            rxim[3] = _mm256_load_pd( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n );
            rxre[4] = _mm256_load_pd( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n );
            rxim[4] = _mm256_load_pd( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n );
            
            rxre[0] = _mm256_mul_pd( rtim[0], rxim[3] );
            rxim[0] = _mm256_mul_pd( rtim[0], rxre[3] );
            rxre[1] = _mm256_mul_pd( rtim[1], rxim[4] );
            rxim[1] = _mm256_mul_pd( rtim[1], rxre[4] );
            
            rxre[0] = _mm256_fmsub_pd( rtre[0], rxre[3], rxre[0] );
            rxim[0] = _mm256_fmadd_pd( rtre[0], rxim[3], rxim[0] );
            rxre[1] = _mm256_fmsub_pd( rtre[1], rxre[4], rxre[1] );
            rxim[1] = _mm256_fmadd_pd( rtre[1], rxim[4], rxim[1] );
            
            rxre[1] = _mm256_sub_pd( rxre[0], rxre[1] );
            rxim[1] = _mm256_sub_pd( rxim[0], rxim[1] );
            
            rxre[0] = _mm256_fmsub_pd( rtwo, rxre[0], rxre[1] );
            rxim[0] = _mm256_fmsub_pd( rtwo, rxim[0], rxim[1] );
            
            rxre[3] = _mm256_load_pd( x+iv1+iv  +2*i*n+6*j*l*n );
            rxim[3] = _mm256_load_pd( x+iv1+iv+n+2*i*n+6*j*l*n );
            
            rxre[2] = _mm256_fmadd_pd( rc31, rxre[0], rxre[3] );
            rxim[2] = _mm256_fmadd_pd( rc31, rxim[0], rxim[3] );
            
            rxre[0] = _mm256_add_pd( rxre[0], rxre[3] );
            rxim[0] = _mm256_add_pd( rxim[0], rxim[3] );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+6*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+6*j*l*n, rxim[0] );
            
            rxre[0] = _mm256_fmadd_pd( rc32, rxim[1], rxre[2] );
            rxim[0] = _mm256_fmsub_pd( rc32, rxre[1], rxim[2] );
            rxim[0] = _mm256_xor_pd( rxim[0], rm00 );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n, rxim[0] );
            
            rxre[0] = _mm256_fmsub_pd( rtwo, rxre[2], rxre[0] );
            rxim[0] = _mm256_fmsub_pd( rtwo, rxim[2], rxim[0] );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n, rxim[0] );
            
          }
        }
      }
      
    }
  
}
#endif