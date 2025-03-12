#include "fft.h"

extern inline void fxzm5a_c( const int n,
                             const int k,
                             const int l,
                             double *restrict x,
                             const double *restrict t )

#if defined( avx ) || defined( avx512 )
{
    
    mmreg rtre[4], rtim[4], rxre[6], rxim[6];
    
    #pragma omp parallel for private ( rtre, rtim, rxre, rxim )
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = broadcast( t+8*j   );
      rtim[0] = broadcast( t+8*j+1 );
      rtre[1] = broadcast( t+8*j+2 );
      rtim[1] = broadcast( t+8*j+3 );
      rtre[2] = broadcast( t+8*j+4 );
      rtim[2] = broadcast( t+8*j+5 );
      rtre[3] = broadcast( t+8*j+6 );
      rtim[3] = broadcast( t+8*j+7 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
            
            rxre[0] = load( x+iv1+iv  +2*i*n+2*l*n+10*j*l*n );
            rxim[0] = load( x+iv1+iv+n+2*i*n+2*l*n+10*j*l*n );
            rxre[3] = load( x+iv1+iv  +2*i*n+4*l*n+10*j*l*n );
            rxim[3] = load( x+iv1+iv+n+2*i*n+4*l*n+10*j*l*n );
            rxre[4] = load( x+iv1+iv  +2*i*n+6*l*n+10*j*l*n );
            rxim[4] = load( x+iv1+iv+n+2*i*n+6*l*n+10*j*l*n );
            rxre[5] = load( x+iv1+iv  +2*i*n+8*l*n+10*j*l*n );
            rxim[5] = load( x+iv1+iv+n+2*i*n+8*l*n+10*j*l*n );
            
            rxre[1] = mul( rtre[0], rxre[0] );
            rxim[1] = mul( rtre[0], rxim[0] );
            rxre[2] = mul( rtre[1], rxre[3] );
            rxim[2] = mul( rtre[1], rxim[3] );
            
            rxre[0] = mul( rtim[0], rxre[0] );
            rxim[0] = mul( rtim[0], rxim[0] );
            rxre[3] = mul( rtim[1], rxre[3] );
            rxim[3] = mul( rtim[1], rxim[3] );
            
            rxre[1] = sub( rxre[1], rxim[0] );
            rxim[1] = add( rxim[1], rxre[0] );
            rxre[2] = sub( rxre[2], rxim[3] );
            rxim[2] = add( rxim[2], rxre[3] );
            
            rxre[3] = mul( rtre[2], rxre[4] );
            rxim[3] = mul( rtre[2], rxim[4] );
            rxre[0] = mul( rtre[3], rxre[5] );
            rxim[0] = mul( rtre[3], rxim[5] );
            
            rxre[4] = mul( rtim[2], rxre[4] );
            rxim[4] = mul( rtim[2], rxim[4] );
            rxre[5] = mul( rtim[3], rxre[5] );
            rxim[5] = mul( rtim[3], rxim[5] );
            
            rxre[3] = sub( rxre[3], rxim[4] );
            rxim[3] = add( rxim[3], rxre[4] );
            rxre[0] = sub( rxre[0], rxim[5] );
            rxim[0] = add( rxim[0], rxre[5] );
            
            rxre[3] = sub( rxre[2], rxre[3] );
            rxim[3] = sub( rxim[2], rxim[3] );
            rxre[0] = sub( rxre[1], rxre[0] );
            rxim[0] = sub( rxim[1], rxim[0] );
            
            rxre[1] = add( rxre[1], rxre[1] );
            rxim[1] = add( rxim[1], rxim[1] );
            rxre[4] = add( rxre[2], rxre[2] );
            rxim[4] = add( rxim[2], rxim[2] );
            
            rxre[2] = mul( rc53,    rxre[3] );
            rxim[2] = mul( rc53,    rxim[3] );
            rxre[5] = mul( rc53,    rxre[0] );
            rxim[5] = mul( rc53,    rxim[0] );
            
            rxre[1] = sub( rxre[1], rxre[0] );
            rxim[1] = sub( rxim[1], rxim[0] );
            rxre[4] = sub( rxre[4], rxre[3] );
            rxim[4] = sub( rxim[4], rxim[3] );
            rxre[2] = add( rxre[2], rxre[0] );
            rxim[2] = add( rxim[2], rxim[0] );
            
            rxre[3] = sub( rxre[5], rxre[3] );
            rxim[3] = sub( rxim[5], rxim[3] );
            rxre[0] = add( rxre[1], rxre[4] );
            rxim[0] = add( rxim[1], rxim[4] );
            
            rxre[1] = sub( rxre[1], rxre[4] );
            rxim[1] = sub( rxim[1], rxim[4] );
            rxre[4] = mul( rc51,    rxre[0] );
            rxim[4] = mul( rc51,    rxim[0] );
            
            rxre[5] = load( x+iv1+iv  +2*i*n+10*j*l*n );
            rxim[5] = load( x+iv1+iv+n+2*i*n+10*j*l*n );
            
            rxre[4] = sub( rxre[5], rxre[4] );
            rxim[4] = sub( rxim[5], rxim[4] );
            rxre[1] = mul( rc52,    rxre[1] );
            rxim[1] = mul( rc52,    rxim[1] );
            
            rxre[1] = sub( rxre[4], rxre[1] );
            rxim[1] = sub( rxim[4], rxim[1] );
            rxre[4] = add( rxre[4], rxre[4] );
            rxim[4] = add( rxim[4], rxim[4] );
            
            rxre[4] = sub( rxre[4], rxre[1] );
            rxim[4] = sub( rxim[4], rxim[1] );
            rxre[0] = add( rxre[0], rxre[5] );
            rxim[0] = add( rxim[0], rxim[5] );
            
            store( x+iv1+iv  +2*i*n+10*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+10*j*l*n, rxim[0] );
            
            rxre[2] = mul( rc54, rxre[2] );
            rxim[2] = mul( rc54, rxim[2] );
            rxre[3] = mul( rc54, rxre[3] );
            rxim[3] = mul( rc54, rxim[3] );
            
            rxre[0] = sub( rxre[1], rxim[3] );
            rxim[0] = add( rxim[1], rxre[3] );
            rxre[5] = sub( rxre[4], rxim[2] );
            rxim[5] = add( rxim[4], rxre[2] );
            
            store( x+iv1+iv  +2*i*n+6*l*n+10*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+6*l*n+10*j*l*n, rxim[0] );
            store( x+iv1+iv  +2*i*n+8*l*n+10*j*l*n, rxre[5] );
            store( x+iv1+iv+n+2*i*n+8*l*n+10*j*l*n, rxim[5] );
            
            rxre[1] = add( rxre[1], rxre[1] );
            rxim[1] = add( rxim[1], rxim[1] );
            rxre[4] = add( rxre[4], rxre[4] );
            rxim[4] = add( rxim[4], rxim[4] );
            
            rxre[1] = sub( rxre[1], rxre[0] );
            rxim[1] = sub( rxim[1], rxim[0] );
            rxre[4] = sub( rxre[4], rxre[5] );
            rxim[4] = sub( rxim[4], rxim[5] );
            
            store( x+iv1+iv  +2*i*n+4*l*n+10*j*l*n, rxre[1] );
            store( x+iv1+iv+n+2*i*n+4*l*n+10*j*l*n, rxim[1] );
            store( x+iv1+iv  +2*i*n+2*l*n+10*j*l*n, rxre[4] );
            store( x+iv1+iv+n+2*i*n+2*l*n+10*j*l*n, rxim[4] );
            
          }
        }
      }
      
    }
    
}
#elif defined( fma ) || defined( avx512fma )
{
    
    mmreg rtre[4], rtim[4], rxre[6], rxim[6];
    
    for ( int j = 0; j < k; j++ ) {
      
      rtre[0] = broadcast( t+8*j   );
      rtim[0] = broadcast( t+8*j+1 );
      rtre[1] = broadcast( t+8*j+2 );
      rtim[1] = broadcast( t+8*j+3 );
      rtre[2] = broadcast( t+8*j+4 );
      rtim[2] = broadcast( t+8*j+5 );
      rtre[3] = broadcast( t+8*j+6 );
      rtim[3] = broadcast( t+8*j+7 );
      
      for ( int i = 0; i < l; i++ ) {
        for ( int iv = 0; iv < n; iv+=16 ) {
          for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
            
            rxre[0] = load( x+iv1+iv  +2*i*n+2*l*n+10*j*l*n );
            rxim[0] = load( x+iv1+iv+n+2*i*n+2*l*n+10*j*l*n );
            rxre[3] = load( x+iv1+iv  +2*i*n+4*l*n+10*j*l*n );
            rxim[3] = load( x+iv1+iv+n+2*i*n+4*l*n+10*j*l*n );
            rxre[4] = load( x+iv1+iv  +2*i*n+6*l*n+10*j*l*n );
            rxim[4] = load( x+iv1+iv+n+2*i*n+6*l*n+10*j*l*n );
            rxre[5] = load( x+iv1+iv  +2*i*n+8*l*n+10*j*l*n );
            rxim[5] = load( x+iv1+iv+n+2*i*n+8*l*n+10*j*l*n );
            
            rxre[1] = mul( rtim[0], rxim[0] );
            rxim[1] = mul( rtim[0], rxre[0] );
            rxre[2] = mul( rtim[1], rxim[3] );
            rxim[2] = mul( rtim[1], rxre[3] );
            
            rxre[1] = fmsub( rtre[0], rxre[0], rxre[1] );
            rxim[1] = fmadd( rtre[0], rxim[0], rxim[1] );
            rxre[2] = fmsub( rtre[1], rxre[3], rxre[2] );
            rxim[2] = fmadd( rtre[1], rxim[3], rxim[2] );
            
            rxre[3] = mul( rtim[2], rxim[4] );
            rxim[3] = mul( rtim[2], rxre[4] );
            rxre[0] = mul( rtim[3], rxim[5] );
            rxim[0] = mul( rtim[3], rxre[5] );
            
            rxre[3] = fmsub( rtre[2], rxre[4], rxre[3] );
            rxim[3] = fmadd( rtre[2], rxim[4], rxim[3] );
            rxre[0] = fmsub( rtre[3], rxre[5], rxre[0] );
            rxim[0] = fmadd( rtre[3], rxim[5], rxim[0] );
            
            rxre[3] = sub( rxre[2], rxre[3] );
            rxim[3] = sub( rxim[2], rxim[3] );
            rxre[0] = sub( rxre[1], rxre[0] );
            rxim[0] = sub( rxim[1], rxim[0] );
            
            rxre[1] = fmsub( rtwo, rxre[1], rxre[0] );
            rxim[1] = fmsub( rtwo, rxim[1], rxim[0] );
            rxre[4] = fmsub( rtwo, rxre[2], rxre[3] );
            rxim[4] = fmsub( rtwo, rxim[2], rxim[3] );
            rxre[2] = fmadd( rc53, rxre[3], rxre[0] );
            rxim[2] = fmadd( rc53, rxim[3], rxim[0] );
            
            rxre[3] = fmsub( rc53, rxre[0], rxre[3] );
            rxim[3] = fmsub( rc53, rxim[0], rxim[3] );
            
            rxre[0] = add( rxre[1], rxre[4] );
            rxim[0] = add( rxim[1], rxim[4] );
            rxre[1] = sub( rxre[1], rxre[4] );
            rxim[1] = sub( rxim[1], rxim[4] );
            
            rxre[5] = load( x+iv1+iv  +2*i*n+10*j*l*n );
            rxim[5] = load( x+iv1+iv+n+2*i*n+10*j*l*n );
            
            rxre[4] = fmsub( rc51, rxre[0], rxre[5] );
            rxim[4] = fmsub( rc51, rxim[0], rxim[5] );
            
            rxre[1] = fmadd( rc52, rxre[1], rxre[4] );
            rxim[1] = fmadd( rc52, rxim[1], rxim[4] );
            
            rxre[4] = fmsub( rtwo, rxre[4], rxre[1] );
            rxim[4] = fmsub( rtwo, rxim[4], rxim[1] );
            
            rxre[1] = xor( rxre[1], rm00 );
            rxim[1] = xor( rxim[1], rm00 );
            rxre[4] = xor( rxre[4], rm00 );
            rxim[4] = xor( rxim[4], rm00 );
            
            rxre[0] = add( rxre[0], rxre[5] );
            rxim[0] = add( rxim[0], rxim[5] );
            
            store( x+iv1+iv  +2*i*n+10*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+10*j*l*n, rxim[0] );
            
            rxre[0] = fmsub( rc54, rxim[3], rxre[1] );
            rxim[0] = fmadd( rc54, rxre[3], rxim[1] );
            rxre[5] = fmsub( rc54, rxim[2], rxre[4] );
            rxim[5] = fmadd( rc54, rxre[2], rxim[4] );
            
            rxre[0] = xor( rxre[0], rm00 );
            rxre[5] = xor( rxre[5], rm00 );
            
            store( x+iv1+iv  +2*i*n+6*l*n+10*j*l*n, rxre[0] );
            store( x+iv1+iv+n+2*i*n+6*l*n+10*j*l*n, rxim[0] );
            store( x+iv1+iv  +2*i*n+8*l*n+10*j*l*n, rxre[5] );
            store( x+iv1+iv+n+2*i*n+8*l*n+10*j*l*n, rxim[5] );
            
            rxre[1] = fmsub( rtwo, rxre[1], rxre[0] );
            rxim[1] = fmsub( rtwo, rxim[1], rxim[0] );
            rxre[4] = fmsub( rtwo, rxre[4], rxre[5] );
            rxim[4] = fmsub( rtwo, rxim[4], rxim[5] );
            
            store( x+iv1+iv  +2*i*n+4*l*n+10*j*l*n, rxre[1] );
            store( x+iv1+iv+n+2*i*n+4*l*n+10*j*l*n, rxim[1] );
            store( x+iv1+iv  +2*i*n+2*l*n+10*j*l*n, rxre[4] );
            store( x+iv1+iv+n+2*i*n+2*l*n+10*j*l*n, rxim[4] );
            
          }
        }
      }
    
    }
    
}
#endif
