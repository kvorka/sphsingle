#include "fft.h"

extern inline void fxzm4a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t )

#if defined( avx ) || defined( avx512 )
{
    
    mmreg rtre[2], rtim[2], rxre[6], rxim[6];
    
    #pragma omp parallel for private (rtre,rtim,rxre,rxim)
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = broadcast( t+6*j   );
      rtim[0] = broadcast( t+6*j+1 );
      rtre[1] = broadcast( t+6*j+2 );
      rtim[1] = broadcast( t+6*j+3 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
            
            rxre[4] = load(x+iv1+iv  +2*i*n+4*l*n+8*j*l*n);
            rxim[4] = load(x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n);
            rxre[5] = load(x+iv1+iv  +2*i*n+6*l*n+8*j*l*n);
            rxim[5] = load(x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n);
            
            rxre[2] = mul( rtim[1], rxim[4] );
            rxim[2] = mul( rtre[1], rxim[4] );
            rxre[3] = mul( rtim[1], rxim[5] );
            rxim[3] = mul( rtre[1], rxim[5] );
            
            rxre[0] = mul( rtre[1], rxre[4] );
            rxim[0] = mul( rtim[1], rxre[4] );
            rxre[1] = mul( rtre[1], rxre[5] );
            rxim[1] = mul( rtim[1], rxre[5] );
            
            rxre[2] = sub( rxre[2], rxre[0] );
            rxim[2] = add( rxim[2], rxim[0] );
            rxre[3] = sub( rxre[3], rxre[1] );
            rxim[3] = add( rxim[3], rxim[1] );
            
            rxre[0] = load(x+iv1+iv  +2*i*n+      8*j*l*n);
            rxim[0] = load(x+iv1+iv+n+2*i*n+      8*j*l*n);
            rxre[1] = load(x+iv1+iv  +2*i*n+2*l*n+8*j*l*n);
            rxim[1] = load(x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n);
            
            rxre[2] = add( rxre[0], rxre[2] );
            rxim[2] = sub( rxim[0], rxim[2] );
            rxre[3] = add( rxre[1], rxre[3] );
            rxim[3] = sub( rxim[1], rxim[3] );
            
            rxre[4] = add( rxre[0], rxre[0] );
            rxim[4] = add( rxim[0], rxim[0] );
            rxre[5] = add( rxre[1], rxre[1] );
            rxim[5] = add( rxim[1], rxim[1] );
            
            rxre[0] = sub( rxre[4], rxre[2] );
            rxim[0] = sub( rxim[4], rxim[2] );
            rxre[1] = sub( rxre[5], rxre[3] );
            rxim[1] = sub( rxim[5], rxim[3] );
            
            rxre[4] = mul( rtim[0], rxim[1] );
            rxim[4] = mul( rtre[0], rxim[1] );
            rxre[5] = mul( rtre[0], rxim[3] );
            rxim[5] = mul( rtre[0], rxre[3] );
            
            rxim[1] = mul( rtre[0], rxre[1] );
            rxre[1] = mul( rtim[0], rxre[1] );
            rxre[3] = mul( rtim[0], rxre[3] );
            rxim[3] = mul( rtim[0], rxim[3] );
            
            rxre[4] = sub( rxre[4], rxim[1] );
            rxim[4] = add( rxim[4], rxre[1] );
            rxre[5] = add( rxre[5], rxre[3] );
            rxim[5] = sub( rxim[5], rxim[3] );
            
            rxre[4] = add( rxre[0], rxre[4] );
            rxim[4] = sub( rxim[0], rxim[4] );
            rxre[5] = sub( rxre[2], rxre[5] );
            rxim[5] = add( rxim[2], rxim[5] );
            
            rxre[0] = add( rxre[0], rxre[0] );
            rxim[0] = add( rxim[0], rxim[0] );
            rxre[2] = add( rxre[2], rxre[2] );
            rxim[2] = add( rxim[2], rxim[2] );
            
            rxre[0] = sub( rxre[0], rxre[4] );
            rxim[0] = sub( rxim[0], rxim[4] );
            rxre[2] = sub( rxre[2], rxre[5] );
            rxim[2] = sub( rxim[2], rxim[5] );
            
            store( x+iv1+iv  +2*i*n+4*l*n+8*j*l*n, rxre[4] );
            store( x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n, rxim[4] );
            store( x+iv1+iv  +2*i*n+      8*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+      8*j*l*n, rxim[0] );
            store( x+iv1+iv  +2*i*n+2*l*n+8*j*l*n, rxre[5] );
            store( x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n, rxim[5] );
            store( x+iv1+iv  +2*i*n+6*l*n+8*j*l*n, rxre[2] );
            store( x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n, rxim[2] );
            
          }
        }
      }
      
    }
    
}
#elif defined( fma ) || defined( avx512fma )
{
    
    mmreg rtre[2], rtim[2], rxre[6], rxim[6];
    
    #pragma omp parallel for private ( rtre, rtim, rxre, rxim )
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = broadcast( t+6*j   );
      rtim[0] = broadcast( t+6*j+1 );
      rtre[1] = broadcast( t+6*j+2 );
      rtim[1] = broadcast( t+6*j+3 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
            
            rxre[4] = load(x+iv1+iv  +2*i*n+4*l*n+8*j*l*n);
            rxim[4] = load(x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n);
            rxre[5] = load(x+iv1+iv  +2*i*n+6*l*n+8*j*l*n);
            rxim[5] = load(x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n);
            
            rxre[0] = load(x+iv1+iv  +2*i*n+      8*j*l*n);
            rxim[0] = load(x+iv1+iv+n+2*i*n+      8*j*l*n);
            rxre[1] = load(x+iv1+iv  +2*i*n+2*l*n+8*j*l*n);
            rxim[1] = load(x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n);
            
            rxre[2] = fmsub( rtre[1], rxre[4], rxre[0] );
            rxim[2] = fmsub( rtim[1], rxre[4], rxim[0] );
            rxre[3] = fmsub( rtre[1], rxre[5], rxre[1] );
            rxim[3] = fmsub( rtim[1], rxre[5], rxim[1] );
            
            rxre[2] = fmsub( rtim[1], rxim[4], rxre[2] );
            rxim[2] = fmadd( rtre[1], rxim[4], rxim[2] );
            rxre[3] = fmsub( rtim[1], rxim[5], rxre[3] );
            rxim[3] = fmadd( rtre[1], rxim[5], rxim[3] );
            
            rxre[0] = fmsub( rtwo, rxre[0], rxre[2] );
            rxim[0] = fmadd( rtwo, rxim[0], rxim[2] );
            rxre[1] = fmsub( rtwo, rxre[1], rxre[3] );
            rxim[1] = fmadd( rtwo, rxim[1], rxim[3] );
            
            rxre[4] = fmsub( rtre[0], rxre[1], rxre[0] );
            rxim[4] = fmsub( rtim[0], rxre[1], rxim[0] );
            rxre[5] = fmadd( rtre[0], rxim[3], rxre[2] );
            rxim[5] = fmsub( rtre[0], rxre[3], rxim[2] );
            
            rxre[4] = fmsub( rtim[0], rxim[1], rxre[4] );
            rxim[4] = fmadd( rtre[0], rxim[1], rxim[4] );
            rxre[5] = fmsub( rtim[0], rxre[3], rxre[5] );
            rxim[5] = fmadd( rtim[0], rxim[3], rxim[5] );
            
            rxre[0] = fmsub( rtwo, rxre[0], rxre[4] );
            rxim[0] = fmadd( rtwo, rxim[0], rxim[4] );
            rxre[2] = fmadd( rtwo, rxre[2], rxre[5] );
            rxim[2] = fmadd( rtwo, rxim[2], rxim[5] );
            
            rxim[2] = xor( rxim[2], rm00 );
            rxim[4] = xor( rxim[4], rm00 );
            rxre[5] = xor( rxre[5], rm00 );
            
            store( x+iv1+iv  +2*i*n+4*l*n+8*j*l*n, rxre[4] );
            store( x+iv1+iv+n+2*i*n+4*l*n+8*j*l*n, rxim[4] );
            store( x+iv1+iv  +2*i*n+      8*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+      8*j*l*n, rxim[0] );
            store( x+iv1+iv  +2*i*n+2*l*n+8*j*l*n, rxre[5] );
            store( x+iv1+iv+n+2*i*n+2*l*n+8*j*l*n, rxim[5] );
            store( x+iv1+iv  +2*i*n+6*l*n+8*j*l*n, rxre[2] );
            store( x+iv1+iv+n+2*i*n+6*l*n+8*j*l*n, rxim[2] );
            
          }
        }
      }
      
    }
    
}
#endif