#include <stdlib.h>
#include <immintrin.h>

extern inline void fxzm4a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t ) {
  
  __m256d rt1re, rt1im, rt2re, rt2im;
  __m256d rxre[6], rxim[6];
  
  for ( int j = 0; j < k; j++) {
    
    rt1re = _mm256_broadcast_sd( t+3*j   );
    rt1im = _mm256_broadcast_sd( t+3*j+1 );
    rt2re = _mm256_broadcast_sd( t+3*j+2 );
    rt2im = _mm256_broadcast_sd( t+3*j+3 );
    
    for ( int i = 0; i < l; i++) {
      for ( int iv = 0; iv < n; iv+=4) {
        
        rxre[0] = _mm256_load_pd(x+iv  +2*i*n+      8*j*l*n);
        rxim[0] = _mm256_load_pd(x+iv+n+2*i*n+      8*j*l*n);
        rxre[1] = _mm256_load_pd(x+iv  +2*i*n+2*l*n+8*j*l*n);
        rxim[1] = _mm256_load_pd(x+iv+n+2*i*n+2*l*n+8*j*l*n);
        
        rxre[2] = _mm256_mul_pd( rt2im, _mm256_load_pd(x+iv+n+2*i*n+4*l*n+8*j*l*n) );
        rxim[2] = _mm256_mul_pd( rt2re, _mm256_load_pd(x+iv+n+2*i*n+4*l*n+8*j*l*n) );
        rxre[3] = _mm256_mul_pd( rt2im, _mm256_load_pd(x+iv+n+2*i*n+6*l*n+8*j*l*n) );
        rxim[3] = _mm256_mul_pd( rt2re, _mm256_load_pd(x+iv+n+2*i*n+6*l*n+8*j*l*n) );
        
        rxre[2] = _mm256_sub_pd( rxre[2], _mm256_mul_pd( rt2re, _mm256_load_pd(x+iv+2*i*n+4*l*n+8*j*l*n) ) );
        rxim[2] = _mm256_add_pd( rxim[2], _mm256_mul_pd( rt2im, _mm256_load_pd(x+iv+2*i*n+4*l*n+8*j*l*n) ) );
        rxre[3] = _mm256_sub_pd( rxre[3], _mm256_mul_pd( rt2re, _mm256_load_pd(x+iv+2*i*n+6*l*n+8*j*l*n) ) );
        rxim[3] = _mm256_add_pd( rxim[3], _mm256_mul_pd( rt2im, _mm256_load_pd(x+iv+2*i*n+6*l*n+8*j*l*n) ) );
        
        rxre[2] = _mm256_add_pd( rxre[0], rxre[2] );
        rxim[2] = _mm256_sub_pd( rxim[0], rxim[2] );
        rxre[3] = _mm256_add_pd( rxre[1], rxre[3] );
        rxim[3] = _mm256_sub_pd( rxim[1], rxim[3] );
        
        rxre[0] = _mm256_sub_pd( _mm256_add_pd( rxre[0], rxre[0] ), rxre[2] );
        rxim[0] = _mm256_sub_pd( _mm256_add_pd( rxim[0], rxim[0] ), rxim[2] );
        rxre[1] = _mm256_sub_pd( _mm256_add_pd( rxre[1], rxre[1] ), rxre[3] );
        rxim[1] = _mm256_sub_pd( _mm256_add_pd( rxim[1], rxim[1] ), rxim[3] );
        
        rxre[5] = _mm256_sub_pd( _mm256_mul_pd( rt1im, rxim[1] ), _mm256_mul_pd( rt1re, rxre[1] ) );
        rxim[5] = _mm256_add_pd( _mm256_mul_pd( rt1re, rxim[1] ), _mm256_mul_pd( rt1im, rxre[1] ) );
        rxre[6] = _mm256_add_pd( _mm256_mul_pd( rt1re, rxim[3] ), _mm256_mul_pd( rt1im, rxre[3] ) );
        rxim[6] = _mm256_sub_pd( _mm256_mul_pd( rt1re, rxre[3] ), _mm256_mul_pd( rt1im, rxim[3] ) );
        
        rxre[5] = _mm256_add_pd( rxre[0], rxre[5] );
        rxim[5] = _mm256_sub_pd( rxim[0], rxim[5] );
        rxre[6] = _mm256_sub_pd( rxre[2], rxre[6] );
        rxim[6] = _mm256_add_pd( rxim[2], rxim[6] );
        
        rxre[0] = _mm256_sub_pd( _mm256_add_pd( rxre[0], rxre[0] ), rxre[5] );
        rxim[0] = _mm256_sub_pd( _mm256_add_pd( rxim[0], rxim[0] ), rxim[5] );
        rxre[2] = _mm256_sub_pd( _mm256_add_pd( rxre[2], rxre[2] ), rxre[6] );
        rxim[2] = _mm256_sub_pd( _mm256_add_pd( rxim[2], rxim[2] ), rxim[6] );
        
        _mm256_store_pd( x+iv+  2*i*n+4*l*n+8*j*l*n, rxre[5] );
        _mm256_store_pd( x+iv+n+2*i*n+4*l*n+8*j*l*n, rxim[5] );
        _mm256_store_pd( x+iv  +2*i*n+      8*j*l*n, rxre[0] );
        _mm256_store_pd( x+iv+n+2*i*n+      8*j*l*n, rxim[0] );
        _mm256_store_pd( x+iv  +2*i*n+2*l*n+8*j*l*n, rxre[6] );
        _mm256_store_pd( x+iv+n+2*i*n+2*l*n+8*j*l*n, rxim[6] );
        _mm256_store_pd( x+iv+  2*i*n+6*l*n+8*j*l*n, rxre[2] );
        _mm256_store_pd( x+iv+n+2*i*n+6*l*n+8*j*l*n, rxim[2] );
        
      }
    }
    
  }
  
}