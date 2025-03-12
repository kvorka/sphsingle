#include "fft.h"

extern inline void fxzm2a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t)

#if defined( avx ) || defined( avx512 )
{
  
  mmreg rtre, rtim, rxre[3], rxim[3];
  
  for ( int j = 0; j < k; j++ ) {
    
    rtre = broadcast( t+2*j   );
    rtim = broadcast( t+2*j+1 );
    
    for ( int i = 0; i < l; i++ ) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
          
          rxre[0] = load( x+iv1+iv  +2*i*n+      4*j*l*n );
          rxim[0] = load( x+iv1+iv+n+2*i*n+      4*j*l*n );
          rxre[2] = load( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n );
          rxim[2] = load( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n );
          
          rxre[1] = mul( rtre, rxre[2] );
          rxim[1] = mul( rtim, rxre[2] );
          
          rxre[1] = sub( rxre[0], rxre[1] );
          rxim[1] = sub( rxim[0], rxim[1] );
          rxre[2] = mul( rtim, rxim[2] );
          rxim[2] = mul( rtre, rxim[2] );
          
          rxre[0] = add( rxre[0], rxre[0] );
          rxim[0] = add( rxim[0], rxim[0] );
          rxre[1] = add( rxre[1], rxre[2] );
          rxim[1] = sub( rxim[1], rxim[2] );
          
          rxre[0] = sub( rxre[0], rxre[1] );
          rxim[0] = sub( rxim[0], rxim[1] );
          
          store( x+iv1+iv  +2*i*n+      4*j*l*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n+      4*j*l*n, rxim[0] );
          store( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n, rxre[1] );
          store( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n, rxim[1] );
          
        }
      }
    }
    
  }
    
}
#elif defined( fma ) || defined( avx512fma )
{
    
    mmreg rtre, rtim, rxre[3], rxim[3];
    
    for ( int j = 0; j < k; j++ ) {
      
      rtre = broadcast( t+2*j   );
      rtim = broadcast( t+2*j+1 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
            
            rxre[0] = load( x+iv1+iv  +2*i*n+      4*j*l*n );
            rxim[0] = load( x+iv1+iv+n+2*i*n+      4*j*l*n );
            rxre[2] = load( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n );
            rxim[2] = load( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n );
            
            rxre[1] = fmsub( rtre, rxre[2], rxre[0] );
            rxim[1] = fmsub( rtim, rxre[2], rxim[0] );
            
            rxre[1] = fmsub( rtim, rxim[2], rxre[1] );
            rxim[1] = fmadd( rtre, rxim[2], rxim[1] );
            rxim[1] = xor( rxim[1], rm00 );
            
            rxre[0] = fmsub( rtwo, rxre[0], rxre[1] );
            rxim[0] = fmsub( rtwo, rxim[0], rxim[1] );
            
            store( x+iv1+iv  +2*i*n+      4*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+      4*j*l*n, rxim[0] );
            store( x+iv1+iv  +2*i*n+2*l*n+4*j*l*n, rxre[1] );
            store( x+iv1+iv+n+2*i*n+2*l*n+4*j*l*n, rxim[1] );
            
          }
        }
      }
      
    }
    
}
#endif