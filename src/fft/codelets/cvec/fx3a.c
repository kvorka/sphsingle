#include "../../../math/cvec.h"

extern inline void fxzm3a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t)

#if defined( avx ) || defined( avx512 )
{
    
    const mmreg rc31 = setdbl( -0.50000000000000000000 );
    const mmreg rc32 = setdbl( +0.86602540378443864676 );
    
    mmreg rtre[2], rtim[2], rxre[5], rxim[5];
    
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = broadcast( t+4*j   );
      rtim[0] = broadcast( t+4*j+1 );
      rtre[1] = broadcast( t+4*j+2 );
      rtim[1] = broadcast( t+4*j+3 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
            
            rxre[3] = load( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n );
            rxim[3] = load( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n );
            rxre[4] = load( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n );
            rxim[4] = load( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n );
            
            rxre[0] = mul( rtre[0], rxre[3] );
            rxim[0] = mul( rtre[0], rxim[3] );
            rxre[2] = mul( rtim[0], rxim[3] );
            rxim[2] = mul( rtim[0], rxre[3] );
            
            rxre[1] = mul( rtre[1], rxre[4] );
            rxim[1] = mul( rtre[1], rxim[4] );
            rxre[3] = mul( rtim[1], rxim[4] );
            rxim[3] = mul( rtim[1], rxre[4] );
            
            rxre[0] = sub( rxre[0], rxre[2] );
            rxim[0] = add( rxim[0], rxim[2] );
            
            rxre[1] = sub( rxre[3], rxre[1] );
            rxim[1] = add( rxim[3], rxim[1] );
            
            rxre[1] = add( rxre[0], rxre[1] );
            rxim[1] = sub( rxim[0], rxim[1] );
            
            rxre[0] = add( rxre[0], rxre[0] );
            rxim[0] = add( rxim[0], rxim[0] );
            
            rxre[0] = sub( rxre[0], rxre[1] );
            rxim[0] = sub( rxim[0], rxim[1] );
            
            rxre[2] = mul( rc31, rxre[0] );
            rxim[2] = mul( rc31, rxim[0] );
            rxre[3] = load( x+iv1+iv  +2*i*n+6*j*l*n );
            rxim[3] = load( x+iv1+iv+n+2*i*n+6*j*l*n );
            
            rxre[0] = add( rxre[0], rxre[3] );
            rxim[0] = add( rxim[0], rxim[3] );
            
            store( x+iv1+iv  +2*i*n+6*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+6*j*l*n, rxim[0] );
            
            rxre[1] = mul( rc32, rxre[1] );
            rxim[1] = mul( rc32, rxim[1] );
            rxre[2] = add( rxre[2], rxre[3] );
            rxim[2] = add( rxim[2], rxim[3] );
            
            rxre[0] = add( rxre[2], rxim[1] );
            rxim[0] = sub( rxim[2], rxre[1] );
            rxre[4] = add( rxre[2], rxre[2] );
            rxim[4] = add( rxim[2], rxim[2] );
            
            store( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n, rxim[0] );
            
            rxre[4] = sub( rxre[4], rxre[0] );
            rxim[4] = sub( rxim[4], rxim[0] );
            
            store( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n, rxre[4] );
            store( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n, rxim[4] );
            
          }
        }
      }
      
    }
    
}
#elif defined( fma )
{
    
    const mmreg rtwo = setdbl( +2.0 );
    const mmreg rm00 = setdbl( -0.0 );
    const mmreg rc31 = setdbl( -0.50000000000000000000 );
    const mmreg rc32 = setdbl( +0.86602540378443864676 );
    
    mmreg rtre[2], rtim[2], rxre[5], rxim[5];
    
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = broadcast( t+4*j   );
      rtim[0] = broadcast( t+4*j+1 );
      rtre[1] = broadcast( t+4*j+2 );
      rtim[1] = broadcast( t+4*j+3 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
            
            rxre[3] = load( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n );
            rxim[3] = load( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n );
            rxre[4] = load( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n );
            rxim[4] = load( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n );
            
            rxre[0] = mul( rtim[0], rxim[3] );
            rxim[0] = mul( rtim[0], rxre[3] );
            rxre[1] = mul( rtim[1], rxim[4] );
            rxim[1] = mul( rtim[1], rxre[4] );
            
            rxre[0] = _mm256_fmsub_pd( rtre[0], rxre[3], rxre[0] );
            rxim[0] = _mm256_fmadd_pd( rtre[0], rxim[3], rxim[0] );
            rxre[1] = _mm256_fmsub_pd( rtre[1], rxre[4], rxre[1] );
            rxim[1] = _mm256_fmadd_pd( rtre[1], rxim[4], rxim[1] );
            
            rxre[1] = sub( rxre[0], rxre[1] );
            rxim[1] = sub( rxim[0], rxim[1] );
            
            rxre[0] = _mm256_fmsub_pd( rtwo, rxre[0], rxre[1] );
            rxim[0] = _mm256_fmsub_pd( rtwo, rxim[0], rxim[1] );
            
            rxre[3] = load( x+iv1+iv  +2*i*n+6*j*l*n );
            rxim[3] = load( x+iv1+iv+n+2*i*n+6*j*l*n );
            
            rxre[2] = _mm256_fmadd_pd( rc31, rxre[0], rxre[3] );
            rxim[2] = _mm256_fmadd_pd( rc31, rxim[0], rxim[3] );
            
            rxre[0] = add( rxre[0], rxre[3] );
            rxim[0] = add( rxim[0], rxim[3] );
            
            store( x+iv1+iv  +2*i*n+6*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+6*j*l*n, rxim[0] );
            
            rxre[0] = _mm256_fmadd_pd( rc32, rxim[1], rxre[2] );
            rxim[0] = _mm256_fmsub_pd( rc32, rxre[1], rxim[2] );
            rxim[0] = _mm256_xor_pd( rxim[0], rm00 );
            
            store( x+iv1+iv  +2*i*n+4*l*n+6*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+4*l*n+6*j*l*n, rxim[0] );
            
            rxre[0] = _mm256_fmsub_pd( rtwo, rxre[2], rxre[0] );
            rxim[0] = _mm256_fmsub_pd( rtwo, rxim[2], rxim[0] );
            
            store( x+iv1+iv  +2*i*n+2*l*n+6*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+2*l*n+6*j*l*n, rxim[0] );
            
          }
        }
      }
      
    }
  
}
#endif