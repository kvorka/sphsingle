#include <stdlib.h>
#include <immintrin.h>

extern inline void fxzm2a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t)

#if defined( avx )
{
  
  __m256d rtre, rtim, rxre[3], rxim[3];
  
  for ( int j = 0; j < k; j++) {
    
    rtre = _mm256_broadcast_sd( t+2*j   );
    rtim = _mm256_broadcast_sd( t+2*j+1 );
    
    for ( int i = 0; i < l; i++) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
          
          rxre[0] = _mm256_load_pd( x+iv1+iv  +2*i*n+      4*j*l*n );
          rxim[0] = _mm256_load_pd( x+iv1+iv+n+2*i*n+      4*j*l*n );
          rxre[2] = _mm256_load_pd( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n );
          rxim[2] = _mm256_load_pd( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n );
          
          rxre[1] = _mm256_mul_pd( rtre, rxre[2] );
          rxim[1] = _mm256_mul_pd( rtim, rxre[2] );
          
          rxre[1] = _mm256_sub_pd( rxre[0], rxre[1] );
          rxim[1] = _mm256_sub_pd( rxim[0], rxim[1] );
          rxre[2] = _mm256_mul_pd( rtim, rxim[2] );
          rxim[2] = _mm256_mul_pd( rtre, rxim[2] );
          
          rxre[0] = _mm256_add_pd( rxre[0], rxre[0] );
          rxim[0] = _mm256_add_pd( rxim[0], rxim[0] );
          rxre[1] = _mm256_add_pd( rxre[1], rxre[2] );
          rxim[1] = _mm256_sub_pd( rxim[1], rxim[2] );
          
          rxre[0] = _mm256_sub_pd( rxre[0], rxre[1] );
          rxim[0] = _mm256_sub_pd( rxim[0], rxim[1] );
          
          _mm256_store_pd( x+iv1+iv  +2*i*n+      4*j*l*n, rxre[0] );
          _mm256_store_pd( x+iv1+iv+n+2*i*n+      4*j*l*n, rxim[0] );
          _mm256_store_pd( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n, rxre[1] );
          _mm256_store_pd( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n, rxim[1] );
          
        }
      }
    }
    
  }
    
}
#elif defined( fma )
{
    
    const __m256d rtwo = _mm256_set1_pd( +2.0 );
    const __m256d rm00 = _mm256_set1_pd( -0.0 );
    
    __m256d rtre, rtim, rxre[3], rxim[3];
    
    for ( int j = 0; j < k; j++) {
      
      rtre = _mm256_broadcast_sd( t+2*j   );
      rtim = _mm256_broadcast_sd( t+2*j+1 );
      
      for ( int i = 0; i < l; i++) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
            
            rxre[0] = _mm256_load_pd( x+iv1+iv  +2*i*n+      4*j*l*n );
            rxim[0] = _mm256_load_pd( x+iv1+iv+n+2*i*n+      4*j*l*n );
            rxre[2] = _mm256_load_pd( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n );
            rxim[2] = _mm256_load_pd( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n );
            
            rxre[1] = _mm256_fmsub_pd( rtre, rxre[2], rxre[0] );
            rxim[1] = _mm256_fmsub_pd( rtim, rxre[2], rxim[0] );
            
            rxre[1] = _mm256_fmsub_pd( rtim, rxim[2], rxre[1] );
            rxim[1] = _mm256_fmadd_pd( rtre, rxim[2], rxim[1] );
            rxim[1] = _mm256_xor_pd( rxim[1], rm00 );
            
            rxre[0] = _mm256_fmsub_pd( rtwo, rxre[0], rxre[1] );
            rxim[0] = _mm256_fmsub_pd( rtwo, rxim[0], rxim[1] );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+      4*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+      4*j*l*n, rxim[0] );
            _mm256_store_pd( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n, rxre[1] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n, rxim[1] );
            
          }
        }
      }
      
    }
    
}
#endif