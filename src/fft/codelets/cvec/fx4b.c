#include "../../../math/cvec.h"

extern inline void fxzm4b_c( const int n,
                             const int l,
                             double *restrict x ) {
  
  mmreg rxre[6], rxim[6];
  
  for ( int i = 0; i < l; i++ ) {
    for ( int iv = 0; iv < n; iv+=16 ) {
      for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
        
        rxre[4] = load( x+iv1+iv+  2*i*n       );
        rxim[4] = load( x+iv1+iv+n+2*i*n       );
        rxre[5] = load( x+iv1+iv+  2*i*n+4*l*n );
        rxim[5] = load( x+iv1+iv+n+2*i*n+4*l*n );
        
        rxre[0] = add( rxre[4], rxre[5] );
        rxim[0] = add( rxim[4], rxim[5] );
        rxre[2] = sub( rxre[4], rxre[5] );
        rxim[2] = sub( rxim[4], rxim[5] );
        
        rxre[4] = load( x+iv1+iv+  2*i*n+2*l*n );
        rxim[4] = load( x+iv1+iv+n+2*i*n+2*l*n );
        rxre[5] = load( x+iv1+iv+  2*i*n+6*l*n );
        rxim[5] = load( x+iv1+iv+n+2*i*n+6*l*n );
        
        rxre[1] = add( rxre[4], rxre[5] );
        rxim[1] = add( rxim[4], rxim[5] );
        rxre[3] = sub( rxre[4], rxre[5] );
        rxim[3] = sub( rxim[4], rxim[5] );
        
        store( x+iv1+iv+  2*i*n      , add( rxre[0], rxre[1] ) );
        store( x+iv1+iv+n+2*i*n      , add( rxim[0], rxim[1] ) );
        store( x+iv1+iv+  2*i*n+2*l*n, sub( rxre[2], rxim[3] ) );
        store( x+iv1+iv+n+2*i*n+2*l*n, add( rxim[2], rxre[3] ) );
        store( x+iv1+iv+  2*i*n+4*l*n, sub( rxre[0], rxre[1] ) );
        store( x+iv1+iv+n+2*i*n+4*l*n, sub( rxim[0], rxim[1] ) );
        store( x+iv1+iv+  2*i*n+6*l*n, add( rxre[2], rxim[3] ) );
        store( x+iv1+iv+n+2*i*n+6*l*n, sub( rxim[2], rxre[3] ) );
        
      }
    }
  }
  
}