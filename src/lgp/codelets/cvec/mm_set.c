#include "../../../math/cvec.h"

extern inline void mm_set_c( const int ma,                  // identifier for m=0 case
                             const int n,                   // howmany roots (step)
                             const double cff,              // recursion coeffs
                             const double *restrict cosx,   // roots
                             const double *restrict sinx,   // sqrt(1-roots**2)
                             double *restrict pmm,          // Lege polys m=0
                             double *restrict pmj1,         // Lege polys for previous step holder
                             double *restrict pmj ) {       // Lege polys
    
    // avx variables for coefficients and polynomials
    mmreg rcff, rpmm;
    
    // set pmj1
    rcff = setzero();
    
    for ( int i2 = 0; i2 < n; i2+=8 ) {
      for ( int i1 = 0; i1 < 8; i1+=incr ) {
      
        store( pmj1+i2+i1, rcff );
        
      }
    }
    
    // set pmm, pmj
    rcff = setdbl( cff );
    
    switch ( ma ) {
      
      case 1:
      {
        
        for ( int i2 = 0; i2 < n; i2+=8 ) {
          for ( int i1 = 0; i1 < 8; i1+=incr ) {
            
            store( pmm+i2+i1, rcff );
            store( pmj+i2+i1, dvv( rcff, load( cosx+i2+i1 ) ) );
            
          }
        }
        
      }
      break;
      
      default:
      {
        
        for ( int i2 = 0; i2 < n; i2+=8 ) {
          for ( int i1 = 0; i1 < 8; i1+=incr ) {
            
            rpmm = mul( rcff, load( sinx+i2+i1 ) );
            rpmm = mul( rpmm, load( pmm +i2+i1 ) );
            
            store( pmm+i2+i1, rpmm );
            store( pmj+i2+i1, dvv( rpmm, load( cosx+i2+i1 ) ) );
            
          }
        }
        
      }
      break;
      
    }
    
}