#include "fft.h"

extern inline void fxzm5b_c( const int n,
                             const int l,
                             double *restrict x )

#if defined( avx ) || defined( avx512 )
{
    
    mmreg rxre[6], rxim[6];
    
    for ( int i = 0; i < l; i++ ) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
          
          rxre[0] = load( x+iv1+iv  +2*i*n+2*l*n );
          rxim[0] = load( x+iv1+iv+n+2*i*n+2*l*n );
          rxre[3] = load( x+iv1+iv  +2*i*n+4*l*n );
          rxim[3] = load( x+iv1+iv+n+2*i*n+4*l*n );
          rxre[5] = load( x+iv1+iv  +2*i*n+6*l*n );
          rxim[5] = load( x+iv1+iv+n+2*i*n+6*l*n );
          rxre[2] = load( x+iv1+iv  +2*i*n+8*l*n );
          rxim[2] = load( x+iv1+iv+n+2*i*n+8*l*n );
          
          rxre[1] = add( rxre[0], rxre[2] );
          rxim[1] = add( rxim[0], rxim[2] );
          rxre[4] = add( rxre[3], rxre[5] );
          rxim[4] = add( rxim[3], rxim[5] );
          
          rxre[0] = sub( rxre[0], rxre[2] );
          rxim[0] = sub( rxim[0], rxim[2] );
          rxre[3] = sub( rxre[3], rxre[5] );
          rxim[3] = sub( rxim[3], rxim[5] );
          
          rxre[2] = mul( rc53, rxre[3] );
          rxim[2] = mul( rc53, rxim[3] );
          rxre[5] = mul( rc53, rxre[0] );
          rxim[5] = mul( rc53, rxim[0] );
          
          rxre[2] = add( rxre[2], rxre[0] );
          rxim[2] = add( rxim[2], rxim[0] );
          rxre[3] = sub( rxre[5], rxre[3] );
          rxim[3] = sub( rxim[5], rxim[3] );
          
          rxre[0] = add( rxre[1], rxre[4] );
          rxim[0] = add( rxim[1], rxim[4] );
          rxre[5] = sub( rxre[1], rxre[4] );
          rxim[5] = sub( rxim[1], rxim[4] );
          
          rxre[5] = mul( rc52, rxre[5] );
          rxim[5] = mul( rc52, rxim[5] );
          rxre[4] = mul( rc51, rxre[0] );
          rxim[4] = mul( rc51, rxim[0] );
          
          rxre[1] = load( x+iv1+iv  +2*i*n );
          rxim[1] = load( x+iv1+iv+n+2*i*n );
          
          rxre[0] = add( rxre[1], rxre[0] );
          rxim[0] = add( rxim[1], rxim[0] );
          rxre[4] = sub( rxre[1], rxre[4] );
          rxim[4] = sub( rxim[1], rxim[4] );
          
          store( x+iv1+iv  +2*i*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n, rxim[0] );
          
          rxre[5] = sub( rxre[4], rxre[5] );
          rxim[5] = sub( rxim[4], rxim[5] );
          rxre[0] = add( rxre[4], rxre[4] );
          rxim[0] = add( rxim[4], rxim[4] );
          
          rxre[2] = mul( rc54, rxre[2] );
          rxim[2] = mul( rc54, rxim[2] );
          rxre[3] = mul( rc54, rxre[3] );
          rxim[3] = mul( rc54, rxim[3] );
          rxre[4] = sub( rxre[0], rxre[5] );
          rxim[4] = sub( rxim[0], rxim[5] );
          
          rxre[0] = sub( rxre[5], rxim[3] );
          rxim[0] = add( rxim[5], rxre[3] );
          rxre[1] = sub( rxre[4], rxim[2] );
          rxim[1] = add( rxim[4], rxre[2] );
          
          store( x+iv1+iv  +2*i*n+6*l*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n+6*l*n, rxim[0] );
          store( x+iv1+iv  +2*i*n+8*l*n, rxre[1] );
          store( x+iv1+iv+n+2*i*n+8*l*n, rxim[1] );
          
          rxre[5] = add( rxre[5], rxre[5] );
          rxim[5] = add( rxim[5], rxim[5] );
          rxre[4] = add( rxre[4], rxre[4] );
          rxim[4] = add( rxim[4], rxim[4] );
          
          rxre[5] = sub( rxre[5], rxre[0] );
          rxim[5] = sub( rxim[5], rxim[0] );
          rxre[4] = sub( rxre[4], rxre[1] );
          rxim[4] = sub( rxim[4], rxim[1] );
          
          store( x+iv1+iv  +2*i*n+2*l*n, rxre[4] );
          store( x+iv1+iv+n+2*i*n+2*l*n, rxim[4] );
          store( x+iv1+iv  +2*i*n+4*l*n, rxre[5] );
          store( x+iv1+iv+n+2*i*n+4*l*n, rxim[5] );
          
        }
      }
    }
    
}
#elif defined( fma ) || defined( avx512fma )
{
    
    mmreg rxre[6], rxim[6];
    
    for ( int i = 0; i < l; i++ ) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
          
          rxre[0] = load( x+iv1+iv  +2*i*n+2*l*n );
          rxim[0] = load( x+iv1+iv+n+2*i*n+2*l*n );
          rxre[3] = load( x+iv1+iv  +2*i*n+4*l*n );
          rxim[3] = load( x+iv1+iv+n+2*i*n+4*l*n );
          rxre[5] = load( x+iv1+iv  +2*i*n+6*l*n );
          rxim[5] = load( x+iv1+iv+n+2*i*n+6*l*n );
          rxre[2] = load( x+iv1+iv  +2*i*n+8*l*n );
          rxim[2] = load( x+iv1+iv+n+2*i*n+8*l*n );
          
          rxre[1] = add( rxre[0], rxre[2] );
          rxim[1] = add( rxim[0], rxim[2] );
          rxre[4] = add( rxre[3], rxre[5] );
          rxim[4] = add( rxim[3], rxim[5] );
          
          rxre[0] = sub( rxre[0], rxre[2] );
          rxim[0] = sub( rxim[0], rxim[2] );
          rxre[5] = sub( rxre[3], rxre[5] );
          rxim[5] = sub( rxim[3], rxim[5] );
          
          rxre[2] = fmadd( rc53, rxre[5], rxre[0] );
          rxim[2] = fmadd( rc53, rxim[5], rxim[0] );
          rxre[3] = fmsub( rc53, rxre[0], rxre[5] );
          rxim[3] = fmsub( rc53, rxim[0], rxim[5] );
          
          rxre[0] = add( rxre[1], rxre[4] );
          rxim[0] = add( rxim[1], rxim[4] );
          rxre[5] = sub( rxre[1], rxre[4] );
          rxim[5] = sub( rxim[1], rxim[4] );
          
          rxre[1] = load( x+iv1+iv  +2*i*n );
          rxim[1] = load( x+iv1+iv+n+2*i*n );
          
          rxre[4] = fmsub( rc51, rxre[0], rxre[1] );
          rxim[4] = fmsub( rc51, rxim[0], rxim[1] );
          
          rxre[5] = fmadd( rc52, rxre[5], rxre[4] );
          rxim[5] = fmadd( rc52, rxim[5], rxim[4] );
          
          rxre[4] = xor( rxre[4], rm00 );
          rxim[4] = xor( rxim[4], rm00 );
          rxre[5] = xor( rxre[5], rm00 );
          rxim[5] = xor( rxim[5], rm00 );
          
          rxre[4] = fmsub( rtwo, rxre[4], rxre[5] );
          rxim[4] = fmsub( rtwo, rxim[4], rxim[5] );
          
          rxre[0] = add( rxre[0], rxre[1] );
          rxim[0] = add( rxim[0], rxim[1] );
          
          store( x+iv1+iv  +2*i*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n, rxim[0] );
          
          rxre[0] = fmsub( rc54, rxim[3], rxre[5] );
          rxim[0] = fmadd( rc54, rxre[3], rxim[5] );
          rxre[1] = fmsub( rc54, rxim[2], rxre[4] );
          rxim[1] = fmadd( rc54, rxre[2], rxim[4] );
          
          rxre[0] = xor( rxre[0], rm00 );
          rxre[1] = xor( rxre[1], rm00 );
          
          store( x+iv1+iv  +2*i*n+6*l*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n+6*l*n, rxim[0] );
          store( x+iv1+iv  +2*i*n+8*l*n, rxre[1] );
          store( x+iv1+iv+n+2*i*n+8*l*n, rxim[1] );
          
          rxre[4] = fmsub( rtwo, rxre[4], rxre[1] );
          rxim[4] = fmsub( rtwo, rxim[4], rxim[1] );
          rxre[5] = fmsub( rtwo, rxre[5], rxre[0] );
          rxim[5] = fmsub( rtwo, rxim[5], rxim[0] );
          
          store( x+iv1+iv  +2*i*n+2*l*n, rxre[4] );
          store( x+iv1+iv+n+2*i*n+2*l*n, rxim[4] );
          store( x+iv1+iv  +2*i*n+4*l*n, rxre[5] );
          store( x+iv1+iv+n+2*i*n+4*l*n, rxim[5] );
          
        }
      }
    }
    
}
#endif