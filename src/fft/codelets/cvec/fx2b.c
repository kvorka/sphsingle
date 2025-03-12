#include "fft.h"

extern inline void fxzm2b_c( const int n,
                             const int l,
                             double *restrict x )

#if defined( avx ) || defined( avx512 )
{
    
    mmreg rxre[2], rxim[2];
    
    for ( int i = 0; i < l; i++ ) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
          
          rxre[0] = load( x+iv+iv1  +2*i*n       );
          rxim[0] = load( x+iv+iv1+n+2*i*n       );
          rxre[1] = load( x+iv+iv1  +2*i*n+2*l*n );
          rxim[1] = load( x+iv+iv1+n+2*i*n+2*l*n );
          
          rxre[1] = sub( rxre[0], rxre[1] );
          rxim[1] = sub( rxre[0], rxim[1] );
          
          rxre[0] = add( rxre[0], rxre[0] );
          rxim[0] = add( rxre[0], rxim[0] );
          
          rxre[0] = sub( rxre[0], rxre[1] );
          rxim[0] = sub( rxre[0], rxim[1] );
          
          store( x+iv+iv1  +2*i*n,       rxre[0] );
          store( x+iv+iv1+n+2*i*n,       rxim[0] );
          store( x+iv+iv1  +2*i*n+2*l*n, rxre[1] );
          store( x+iv+iv1+n+2*i*n+2*l*n, rxim[1] );
          
        }
      }
    }
    
}
#elif defined( fma ) || defined( avx512fma )
{
    
    mmreg rxre[2], rxim[2];
    
    for ( int i = 0; i < l; i++ ) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
          
          rxre[0] = load( x+iv+iv1  +2*i*n       );
          rxim[0] = load( x+iv+iv1+n+2*i*n       );
          rxre[1] = load( x+iv+iv1  +2*i*n+2*l*n );
          rxim[1] = load( x+iv+iv1+n+2*i*n+2*l*n );
          
          rxre[1] = sub( rxre[0], rxre[1] );
          rxim[1] = sub( rxre[0], rxim[1] );
          
          rxre[0] = fmsub( rtwo, rxre[0], rxre[1] );
          rxim[0] = fmsub( rtwo, rxre[0], rxim[1] );
          
          store( x+iv+iv1  +2*i*n,       rxre[0] );
          store( x+iv+iv1+n+2*i*n,       rxim[0] );
          store( x+iv+iv1  +2*i*n+2*l*n, rxre[1] );
          store( x+iv+iv1+n+2*i*n+2*l*n, rxim[1] );
          
        }
      }
    }
    
}
#endif