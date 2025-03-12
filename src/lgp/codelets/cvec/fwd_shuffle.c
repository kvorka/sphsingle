#include "../../../math/cvec.h"

extern inline void fwd_shuffle_c( const int n,                  // howmany roots (step)
                                  const double *restrict cosx,  // roots
                                  const double *restrict wght,  // weights
                                  double *restrict grid,        // Legendre sums
                                  double *restrict swork ) {    // partial sum

    
    // constant needed for loop unrolling
    const int n32 = (n/32)*32;
    
    // avx vars for cosine and weight values, North and South sums
    mmreg rwcsx, rwght, rsumN[2], rsumS[2];
    
    // fwd shuffle: cycle over the roots, the outer cycle is unrolled by 16, factor of 4 is handled by an explicit
    // vectorization, factor of 4 is added in order to unroll the cycle a bit more for efficiency
    for ( int i1 = 0; i1 < n32; i1+=32 ) {
      for ( int i2 = 0; i2 < 32; i2+=incr ) {
        
        rwght =      load( wght+i1+i2 );
        rwcsx = mul( load( cosx+i1+i2 ), rwght );
        
        rsumN[0] = load( grid+i1+i2     );
        rsumS[0] = load( grid+i1+i2+  n );
        rsumN[1] = load( grid+i1+i2+2*n );
        rsumS[1] = load( grid+i1+i2+3*n );
        
        store( swork+i1+i2    , mul( sub( rsumN[0], rsumS[0] ), rwght ) );
        store( swork+i1+i2+  n, mul( sub( rsumN[1], rsumS[1] ), rwght ) );
        store( swork+i1+i2+2*n, mul( add( rsumN[0], rsumS[0] ), rwcsx ) );
        store( swork+i1+i2+3*n, mul( add( rsumN[1], rsumS[1] ), rwcsx ) );
        
      }
    }
    
    // fwd shuffle: remainder cases
    for ( int i2 = 0; i2 < n-n32; i2+=incr ) {
        
      rwght =      load( wght+n32+i2 );
      rwcsx = mul( load( cosx+n32+i2 ), rwght );
      
      rsumN[0] = load( grid+n32+i2     );
      rsumS[0] = load( grid+n32+i2+  n );
      rsumN[1] = load( grid+n32+i2+2*n );
      rsumS[1] = load( grid+n32+i2+3*n );
      
      store( swork+n32+i2    , mul( sub( rsumN[0], rsumS[0] ), rwght ) );
      store( swork+n32+i2+  n, mul( sub( rsumN[1], rsumS[1] ), rwght ) );
      store( swork+n32+i2+2*n, mul( add( rsumN[0], rsumS[0] ), rwcsx ) );
      store( swork+n32+i2+3*n, mul( add( rsumN[1], rsumS[1] ), rwcsx ) );
      
    }
    
}