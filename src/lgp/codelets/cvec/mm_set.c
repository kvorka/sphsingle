#include "../../../math/cvec.h"

extern inline void mm_set_c( const int n,                   // howmany roots (step)
                             const double *restrict cosx,   // roots
                             const double *restrict pmm,    // Lege polys m=0
                             double *restrict pmj1,         // Lege polys for previous step holder
                             double *restrict pmj ) {       // Lege polys
    
    // set pmj1
    for ( int i2 = 0; i2 < n; i2+=incr ) { store( pmj1+i2, setzero() ); }
    
    // set pmj
    for ( int i2 = 0; i2 < n; i2+=incr ) { store( pmj+i2, dvv( load( pmm+i2 ), load( cosx+i2 ) ) ); }
    
}