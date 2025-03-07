#include <stdlib.h>
#include <immintrin.h>

extern inline void fxzm4a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t )

#if defined( avx )
{
    
    __m256d rtre[2], rtim[2], rxre[6], rxim[6];
    
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = _mm256_broadcast_sd( t+6*j   );
      rtim[0] = _mm256_broadcast_sd( t+6*j+1 );
      rtre[1] = _mm256_broadcast_sd( t+6*j+2 );
      rtim[1] = _mm256_broadcast_sd( t+6*j+3 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
            
            rxre[4] = _mm256_load_pd(x+iv1+iv  +2*i*n+4*l*n+8*j*l*n);
            rxim[4] = _mm256_load_pd(x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n);
            rxre[5] = _mm256_load_pd(x+iv1+iv  +2*i*n+6*l*n+8*j*l*n);
            rxim[5] = _mm256_load_pd(x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n);
            
            rxre[2] = _mm256_mul_pd( rtim[1], rxim[4] );
            rxim[2] = _mm256_mul_pd( rtre[1], rxim[4] );
            rxre[3] = _mm256_mul_pd( rtim[1], rxim[5] );
            rxim[3] = _mm256_mul_pd( rtre[1], rxim[5] );
            
            rxre[0] = _mm256_mul_pd( rtre[1], rxre[4] );
            rxim[0] = _mm256_mul_pd( rtim[1], rxre[4] );
            rxre[1] = _mm256_mul_pd( rtre[1], rxre[5] );
            rxim[1] = _mm256_mul_pd( rtim[1], rxre[5] );
            
            rxre[2] = _mm256_sub_pd( rxre[2], rxre[0] );
            rxim[2] = _mm256_add_pd( rxim[2], rxim[0] );
            rxre[3] = _mm256_sub_pd( rxre[3], rxre[1] );
            rxim[3] = _mm256_add_pd( rxim[3], rxim[1] );
            
            rxre[0] = _mm256_load_pd(x+iv1+iv  +2*i*n+      8*j*l*n);
            rxim[0] = _mm256_load_pd(x+iv1+iv+n+2*i*n+      8*j*l*n);
            rxre[1] = _mm256_load_pd(x+iv1+iv  +2*i*n+2*l*n+8*j*l*n);
            rxim[1] = _mm256_load_pd(x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n);
            
            rxre[2] = _mm256_add_pd( rxre[0], rxre[2] );
            rxim[2] = _mm256_sub_pd( rxim[0], rxim[2] );
            rxre[3] = _mm256_add_pd( rxre[1], rxre[3] );
            rxim[3] = _mm256_sub_pd( rxim[1], rxim[3] );
            
            rxre[4] = _mm256_add_pd( rxre[0], rxre[0] );
            rxim[4] = _mm256_add_pd( rxim[0], rxim[0] );
            rxre[5] = _mm256_add_pd( rxre[1], rxre[1] );
            rxim[5] = _mm256_add_pd( rxim[1], rxim[1] );
            
            rxre[0] = _mm256_sub_pd( rxre[4], rxre[2] );
            rxim[0] = _mm256_sub_pd( rxim[4], rxim[2] );
            rxre[1] = _mm256_sub_pd( rxre[5], rxre[3] );
            rxim[1] = _mm256_sub_pd( rxim[5], rxim[3] );
            
            rxre[4] = _mm256_mul_pd( rtim[0], rxim[1] );
            rxim[4] = _mm256_mul_pd( rtre[0], rxim[1] );
            rxre[5] = _mm256_mul_pd( rtre[0], rxim[3] );
            rxim[5] = _mm256_mul_pd( rtre[0], rxre[3] );
            
            rxim[1] = _mm256_mul_pd( rtre[0], rxre[1] );
            rxre[1] = _mm256_mul_pd( rtim[0], rxre[1] );
            rxre[3] = _mm256_mul_pd( rtim[0], rxre[3] );
            rxim[3] = _mm256_mul_pd( rtim[0], rxim[3] );
            
            rxre[4] = _mm256_sub_pd( rxre[4], rxim[1] );
            rxim[4] = _mm256_add_pd( rxim[4], rxre[1] );
            rxre[5] = _mm256_add_pd( rxre[5], rxre[3] );
            rxim[5] = _mm256_sub_pd( rxim[5], rxim[3] );
            
            rxre[4] = _mm256_add_pd( rxre[0], rxre[4] );
            rxim[4] = _mm256_sub_pd( rxim[0], rxim[4] );
            rxre[5] = _mm256_sub_pd( rxre[2], rxre[5] );
            rxim[5] = _mm256_add_pd( rxim[2], rxim[5] );
            
            rxre[0] = _mm256_add_pd( rxre[0], rxre[0] );
            rxim[0] = _mm256_add_pd( rxim[0], rxim[0] );
            rxre[2] = _mm256_add_pd( rxre[2], rxre[2] );
            rxim[2] = _mm256_add_pd( rxim[2], rxim[2] );
            
            rxre[0] = _mm256_sub_pd( rxre[0], rxre[4] );
            rxim[0] = _mm256_sub_pd( rxim[0], rxim[4] );
            rxre[2] = _mm256_sub_pd( rxre[2], rxre[5] );
            rxim[2] = _mm256_sub_pd( rxim[2], rxim[5] );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+4*l*n+8*j*l*n, rxre[4] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n, rxim[4] );
            _mm256_store_pd( x+iv1+iv  +2*i*n+      8*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+      8*j*l*n, rxim[0] );
            _mm256_store_pd( x+iv1+iv  +2*i*n+2*l*n+8*j*l*n, rxre[5] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n, rxim[5] );
            _mm256_store_pd( x+iv1+iv  +2*i*n+6*l*n+8*j*l*n, rxre[2] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n, rxim[2] );
            
          }
        }
      }
      
    }
    
}
#elif defined( fma )
{
    
    const __m256d rtwo = _mm256_set1_pd( +2.0 );
    const __m256d rm00 = _mm256_set1_pd( -0.0 );
    
    __m256d rtre[2], rtim[2], rxre[6], rxim[6];
    
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = _mm256_broadcast_sd( t+6*j   );
      rtim[0] = _mm256_broadcast_sd( t+6*j+1 );
      rtre[1] = _mm256_broadcast_sd( t+6*j+2 );
      rtim[1] = _mm256_broadcast_sd( t+6*j+3 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
            
            rxre[4] = _mm256_load_pd(x+iv1+iv  +2*i*n+4*l*n+8*j*l*n);
            rxim[4] = _mm256_load_pd(x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n);
            rxre[5] = _mm256_load_pd(x+iv1+iv  +2*i*n+6*l*n+8*j*l*n);
            rxim[5] = _mm256_load_pd(x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n);
            
            rxre[0] = _mm256_load_pd(x+iv1+iv  +2*i*n+      8*j*l*n);
            rxim[0] = _mm256_load_pd(x+iv1+iv+n+2*i*n+      8*j*l*n);
            rxre[1] = _mm256_load_pd(x+iv1+iv  +2*i*n+2*l*n+8*j*l*n);
            rxim[1] = _mm256_load_pd(x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n);
            
            rxre[2] = _mm256_fmsub_pd( rtre[1], rxre[4], rxre[0] );
            rxim[2] = _mm256_fmsub_pd( rtim[1], rxre[4], rxim[0] );
            rxre[3] = _mm256_fmsub_pd( rtre[1], rxre[5], rxre[1] );
            rxim[3] = _mm256_fmsub_pd( rtim[1], rxre[5], rxim[1] );
            
            rxre[2] = _mm256_fmsub_pd( rtim[1], rxim[4], rxre[2] );
            rxim[2] = _mm256_fmadd_pd( rtre[1], rxim[4], rxim[2] );
            rxre[3] = _mm256_fmsub_pd( rtim[1], rxim[5], rxre[3] );
            rxim[3] = _mm256_fmadd_pd( rtre[1], rxim[5], rxim[3] );
            
            rxre[0] = _mm256_fmsub_pd( rtwo, rxre[0], rxre[2] );
            rxim[0] = _mm256_fmadd_pd( rtwo, rxim[0], rxim[2] );
            rxre[1] = _mm256_fmsub_pd( rtwo, rxre[1], rxre[3] );
            rxim[1] = _mm256_fmadd_pd( rtwo, rxim[1], rxim[3] );
            
            rxre[4] = _mm256_fmsub_pd( rtre[0], rxre[1], rxre[0] );
            rxim[4] = _mm256_fmsub_pd( rtim[0], rxre[1], rxim[0] );
            rxre[5] = _mm256_fmadd_pd( rtre[0], rxim[3], rxre[2] );
            rxim[5] = _mm256_fmsub_pd( rtre[0], rxre[3], rxim[2] );
            
            rxre[4] = _mm256_fmsub_pd( rtim[0], rxim[1], rxre[4] );
            rxim[4] = _mm256_fmadd_pd( rtre[0], rxim[1], rxim[4] );
            rxre[5] = _mm256_fmsub_pd( rtim[0], rxre[3], rxre[5] );
            rxim[5] = _mm256_fmadd_pd( rtim[0], rxim[3], rxim[5] );
            
            rxre[0] = _mm256_fmsub_pd( rtwo, rxre[0], rxre[4] );
            rxim[0] = _mm256_fmadd_pd( rtwo, rxim[0], rxim[4] );
            rxre[2] = _mm256_fmadd_pd( rtwo, rxre[2], rxre[5] );
            rxim[2] = _mm256_fmadd_pd( rtwo, rxim[2], rxim[5] );
            
            rxim[2] = _mm256_xor_pd( rxim[2], rm00 );
            rxim[4] = _mm256_xor_pd( rxim[4], rm00 );
            rxre[5] = _mm256_xor_pd( rxre[5], rm00 );
            
            _mm256_store_pd( x+iv1+iv  +2*i*n+4*l*n+8*j*l*n, rxre[4] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n, rxim[4] );
            _mm256_store_pd( x+iv1+iv  +2*i*n+      8*j*l*n, rxre[0] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+      8*j*l*n, rxim[0] );
            _mm256_store_pd( x+iv1+iv  +2*i*n+2*l*n+8*j*l*n, rxre[5] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n, rxim[5] );
            _mm256_store_pd( x+iv1+iv  +2*i*n+6*l*n+8*j*l*n, rxre[2] );
            _mm256_store_pd( x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n, rxim[2] );
            
          }
        }
      }
      
    }
    
}
#endif