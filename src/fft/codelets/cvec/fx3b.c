#include "../../../math/cvec.h"

extern inline void fxzm3b_c( const int n,
                             const int l,
                             double *restrict x )

#if defined( avx ) || defined( avx512 )
{
    
    const mmreg rc31 = setdbl( -0.50000000000000000000 );
    const mmreg rc32 = setdbl( +0.86602540378443864676 );
    
    mmreg rxre[4], rxim[4];
    
    for ( int i = 0; i < l; i++ ) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
          
          rxre[2] = load( x+iv1+iv+  2*i*n+2*l*n );
          rxim[2] = load( x+iv1+iv+n+2*i*n+2*l*n );
          rxre[3] = load( x+iv1+iv+  2*i*n+4*l*n );
          rxim[3] = load( x+iv1+iv+n+2*i*n+4*l*n );
          
          rxre[0] = add( rxre[2], rxre[3] );
          rxim[0] = add( rxim[2], rxim[3] );
          rxre[1] = sub( rxre[2], rxre[3] );
          rxim[1] = sub( rxim[2], rxim[3] );
          
          rxre[1] = mul( rc32, rxre[1] );
          rxim[1] = mul( rc32, rxim[1] );
          rxre[2] = mul( rc31, rxre[0] );
          rxim[2] = mul( rc31, rxim[0] );
          
          rxre[3] = load( x+iv1+iv+  2*i*n );
          rxim[3] = load( x+iv1+iv+n+2*i*n );
          
          rxre[0] = add( rxre[0], rxre[3] );
          rxim[0] = add( rxim[0], rxim[3] );
          rxre[2] = add( rxre[2], rxre[3] );
          rxim[2] = add( rxim[2], rxim[3] );
          
          store( x+iv1+iv+  2*i*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n, rxim[0] );
          
          rxre[0] = sub( rxre[2], rxim[1] );
          rxim[0] = add( rxim[2], rxre[1] );
          rxre[3] = add( rxre[2], rxim[1] );
          rxim[3] = sub( rxim[2], rxre[1] );
          
          store( x+iv1+iv+  2*i*n+2*l*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n+2*l*n, rxim[0] );
          store( x+iv1+iv+  2*i*n+4*l*n, rxre[3] );
          store( x+iv1+iv+n+2*i*n+4*l*n, rxim[3] );
          
        }
      }
    }
    
}
#elif defined (fma)
{
    
    const mmreg rc31 = setdbl( -0.50000000000000000000 );
    const mmreg rc32 = setdbl( +0.86602540378443864676 );
    const mmreg rm00 = setdbl( -0.0 );
    
    mmreg rxre[4], rxim[4];
    
    for ( int i = 0; i < l; i++ ) {
      for ( int iv = 0; iv < n; iv+=16 ) {
        for ( int iv1 = 0; iv1 < 16; iv1+=incr ) {
          
          rxre[2] = load( x+iv1+iv+  2*i*n+2*l*n );
          rxim[2] = load( x+iv1+iv+n+2*i*n+2*l*n );
          rxre[3] = load( x+iv1+iv+  2*i*n+4*l*n );
          rxim[3] = load( x+iv1+iv+n+2*i*n+4*l*n );
          
          rxre[0] = add( rxre[2], rxre[3] );
          rxim[0] = add( rxim[2], rxim[3] );
          rxre[1] = sub( rxre[2], rxre[3] );
          rxim[1] = sub( rxim[2], rxim[3] );
          
          rxre[3] = load( x+iv1+iv+  2*i*n );
          rxim[3] = load( x+iv1+iv+n+2*i*n );
          
          rxre[2] = _mm256_fmadd_pd( rc31, rxre[0], rxre[3] );
          rxim[2] = _mm256_fmadd_pd( rc31, rxim[0], rxim[3] );
          
          rxre[0] = add( rxre[0], rxre[3] );
          rxim[0] = add( rxim[0], rxim[3] );
          
          store( x+iv1+iv+  2*i*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n, rxim[0] );
          
          rxre[0] = _mm256_fmsub_pd( rc32, rxim[1], rxre[2] );
          rxim[0] = _mm256_fmadd_pd( rc32, rxre[1], rxim[2] );
          rxre[3] = _mm256_fmadd_pd( rc32, rxim[1], rxre[2] );
          rxim[3] = _mm256_fmsub_pd( rc32, rxre[1], rxim[2] );
          
          rxre[0] = _mm256_xor_pd( rxre[0], rm00 );
          rxim[3] = _mm256_xor_pd( rxim[3], rm00 );
          
          store( x+iv1+iv+  2*i*n+2*l*n, rxre[0] );
          store( x+iv1+iv+n+2*i*n+2*l*n, rxim[0] );
          store( x+iv1+iv+  2*i*n+4*l*n, rxre[3] );
          store( x+iv1+iv+n+2*i*n+4*l*n, rxim[3] );
          
        }
      }
    }
    
}
#endif