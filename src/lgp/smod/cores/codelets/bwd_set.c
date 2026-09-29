#include "../../../../math/cvec.h"

extern inline __attribute__((always_inline))
void bwd_set_c( const int n,
                const double *restrict cc,
                const double *restrict pmm,
                      double *restrict pmj1,
                      double *restrict pmj,
                      double *restrict swork )

{
    
    // Constants
    const __td r00 = _t_setzero_pd();
    const __td rc0 = _t_set1_pd( *( cc + 0 ) );
    const __td rc1 = _t_set1_pd( *( cc + 1 ) );
    const __td rc2 = _t_set1_pd( *( cc + 2 ) );
    const __td rc3 = _t_set1_pd( *( cc + 3 ) );
    
    // Memory references
    const double *restrict pm = pmm;
          double *restrict p1 = pmj1;
          double *restrict pj = pmj;
          double *restrict ps = swork;
    
    // Registers to be used
    __td rp0, rp1, rs0, rs1, rs2, rs3;
    
    // Main cycle
    for ( int i = 0; i < n; i++ ) {
        
        _t_store_pd( p1 + vlen0, r00 );
        _t_store_pd( p1 + vlen1, r00 );
        
        rp0 = _t_load_pd( pm + vlen0 );
        rp1 = _t_load_pd( pm + vlen1 );
        
        _t_store_pd( pj + vlen0, rp0 );
        _t_store_pd( pj + vlen1, rp1 );
        
        rs0 = _t_mul_pd( rp0, rc0 );
        rs1 = _t_mul_pd( rp0, rc1 );
        
        _t_store_pd( ps + vlen0, rs0 );
        _t_store_pd( ps + vlen1, rs1 );
        
        rs2 = _t_mul_pd( rp0, rc2 );
        rs3 = _t_mul_pd( rp0, rc3 );
        
        _t_store_pd( ps + vlen2, rs2 );
        _t_store_pd( ps + vlen3, rs3 );
        
        rs0 = _t_mul_pd( rp1, rc0 );
        rs1 = _t_mul_pd( rp1, rc1 );
        
        _t_store_pd( ps + vlen4, rs0 );
        _t_store_pd( ps + vlen5, rs1 );
        
        rs2 = _t_mul_pd( rp1, rc2 );
        rs3 = _t_mul_pd( rp1, rc3 );
        
        _t_store_pd( ps + vlen6, rs2 );
        _t_store_pd( ps + vlen7, rs3 );
        
        p1 += vlen2;
        pj += vlen2;
        pm += vlen2;
        ps += vlen8;
        
    }
    
}