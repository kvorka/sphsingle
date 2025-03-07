#include <stdlib.h>
#include <immintrin.h>

extern inline void fxzm2b_c( const int n,
                             const int l,
                             double *restrict x )

#if defined( avx )
{
    
    __m256d rxre[2], rxim[2];
    
    for ( int i = 0; i < l; i++) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
          
          rxre[0] = _mm256_load_pd( x+iv+iv1  +2*i*n       );
          rxim[0] = _mm256_load_pd( x+iv+iv1+n+2*i*n       );
          rxre[1] = _mm256_load_pd( x+iv+iv1  +2*i*n+2*l*n );
          rxim[1] = _mm256_load_pd( x+iv+iv1+n+2*i*n+2*l*n );
          
          rxre[1] = _mm256_sub_pd( rxre[0], rxre[1] );
          rxim[1] = _mm256_sub_pd( rxre[0], rxim[1] );
          
          rxre[0] = _mm256_add_pd( rxre[0], rxre[0] );
          rxim[0] = _mm256_add_pd( rxre[0], rxim[0] );
          
          rxre[0] = _mm256_sub_pd( rxre[0], rxre[1] );
          rxim[0] = _mm256_sub_pd( rxre[0], rxim[1] );
          
          _mm256_store_pd( x+iv+iv1  +2*i*n,       rxre[0] );
          _mm256_store_pd( x+iv+iv1+n+2*i*n,       rxim[0] );
          _mm256_store_pd( x+iv+iv1  +2*i*n+2*l*n, rxre[1] );
          _mm256_store_pd( x+iv+iv1+n+2*i*n+2*l*n, rxim[1] );
          
        }
      }
    }
    
}
#elif defined( fma )
{
    
    const __m256d rtwo = _mm256_set1_pd( +2.0 );
    
    __m256d rxre[2], rxim[2];
    
    for ( int i = 0; i < l; i++) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=4 ) {
          
          rxre[0] = _mm256_load_pd( x+iv+iv1  +2*i*n       );
          rxim[0] = _mm256_load_pd( x+iv+iv1+n+2*i*n       );
          rxre[1] = _mm256_load_pd( x+iv+iv1  +2*i*n+2*l*n );
          rxim[1] = _mm256_load_pd( x+iv+iv1+n+2*i*n+2*l*n );
          
          rxre[1] = _mm256_sub_pd( rxre[0], rxre[1] );
          rxim[1] = _mm256_sub_pd( rxre[0], rxim[1] );
          
          rxre[0] = _mm256_fmsub_pd( rtwo, rxre[0], rxre[1] );
          rxim[0] = _mm256_fmsub_pd( rtwo, rxre[0], rxim[1] );
          
          _mm256_store_pd( x+iv+iv1  +2*i*n,       rxre[0] );
          _mm256_store_pd( x+iv+iv1+n+2*i*n,       rxim[0] );
          _mm256_store_pd( x+iv+iv1  +2*i*n+2*l*n, rxre[1] );
          _mm256_store_pd( x+iv+iv1+n+2*i*n+2*l*n, rxim[1] );
          
        }
      }
    }
    
}
#endif