#pragma once
#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void bwd_set_c( const int m,
                const int n,
                const double *restrict fmj,
                const double *restrict cosx,
                const double *restrict sinx,
                const double *restrict cc,
                      double *restrict pmm,
                      double *restrict pmj1,
                      double *restrict pmj,
                      double *restrict swork )

{
    
    // Constants
    const __td r00 = _t_setzero_pd();
    const __td rff = _t_set1_pd( *( fmj    ) );
    const __td rc0 = _t_set1_pd( *( cc + 0 ) );
    const __td rc1 = _t_set1_pd( *( cc + 1 ) );
    const __td rc2 = _t_set1_pd( *( cc + 2 ) );
    const __td rc3 = _t_set1_pd( *( cc + 3 ) );
    
    // Memory references
    const double *restrict pcx = cosx;
    const double *restrict psx = sinx;
          double *restrict pm  = pmm;
          double *restrict p1  = pmj1;
          double *restrict pj  = pmj;
          double *restrict psw = swork;
    
    // Registers to be used
    __td rm0, rm1, rp0, rp1, rs0, rs1, rs2, rs3;
    
    // Main cycle
    switch ( m ) {
        
        case 0:
            
            for ( int i = 0; i < n; i++ ) {
                
                _t_store_pd( p1 + vlen0, r00 );
                _t_store_pd( p1 + vlen1, r00 );
                
                rs2 = _t_load_pd( pcx + vlen0 );
                rs3 = _t_load_pd( pcx + vlen1 );
                
                _t_store_pd( pm + vlen0, rff );
                _t_store_pd( pm + vlen1, rff );
                
                rp0 = _t_div_pd( rff, rs2 );
                rp1 = _t_div_pd( rff, rs3 );
                
                _t_store_pd( pj + vlen0, rp0 );
                _t_store_pd( pj + vlen1, rp1 );
                
                rs0 = _t_mul_pd( rp0, rc0 );
                rs1 = _t_mul_pd( rp0, rc1 );
                
                _t_store_pd( psw + vlen0, rs0 );
                _t_store_pd( psw + vlen1, rs1 );
                
                rs2 = _t_mul_pd( rp0, rc2 );
                rs3 = _t_mul_pd( rp0, rc3 );
                
                _t_store_pd( psw + vlen2, rs2 );
                _t_store_pd( psw + vlen3, rs3 );
                
                rs0 = _t_mul_pd( rp1, rc0 );
                rs1 = _t_mul_pd( rp1, rc1 );
                
                _t_store_pd( psw + vlen4, rs0 );
                _t_store_pd( psw + vlen5, rs1 );
                
                rs2 = _t_mul_pd( rp1, rc2 );
                rs3 = _t_mul_pd( rp1, rc3 );
                
                _t_store_pd( psw + vlen6, rs2 );
                _t_store_pd( psw + vlen7, rs3 );
                
                p1  += vlen2;
                pj  += vlen2;
                pm  += vlen2;
                pcx += vlen2;
                psw += vlen8;
                
            }
            
        break;
        
        default:
            
            for ( int i = 0; i < n; i++ ) {
                
                _t_store_pd( p1 + vlen0, r00 );
                _t_store_pd( p1 + vlen1, r00 );
                
                rm0 = _t_load_pd( pm + vlen0 );
                rm1 = _t_load_pd( pm + vlen1 );
                
                rm0 = _t_mul_pd( rm0, rff );
                rm1 = _t_mul_pd( rm1, rff );
                
                rs0 = _t_load_pd( psx + vlen0 );
                rs1 = _t_load_pd( psx + vlen1 );
                
                rm0 = _t_mul_pd( rs0, rm0 );
                rm1 = _t_mul_pd( rs1, rm1 );
                
                rs2 = _t_load_pd( pcx + vlen0 );
                rs3 = _t_load_pd( pcx + vlen1 );
                
                _t_store_pd( pm + vlen0, rm0 );
                _t_store_pd( pm + vlen1, rm1 );
                
                rp0 = _t_div_pd( rm0, rs2 );
                rp1 = _t_div_pd( rm1, rs3 );
                
                _t_store_pd( pj + vlen0, rp0 );
                _t_store_pd( pj + vlen1, rp1 );
                
                rs0 = _t_mul_pd( rp0, rc0 );
                rs1 = _t_mul_pd( rp0, rc1 );
                
                _t_store_pd( psw + vlen0, rs0 );
                _t_store_pd( psw + vlen1, rs1 );
                
                rs2 = _t_mul_pd( rp0, rc2 );
                rs3 = _t_mul_pd( rp0, rc3 );
                
                _t_store_pd( psw + vlen2, rs2 );
                _t_store_pd( psw + vlen3, rs3 );
                
                rs0 = _t_mul_pd( rp1, rc0 );
                rs1 = _t_mul_pd( rp1, rc1 );
                
                _t_store_pd( psw + vlen4, rs0 );
                _t_store_pd( psw + vlen5, rs1 );
                
                rs2 = _t_mul_pd( rp1, rc2 );
                rs3 = _t_mul_pd( rp1, rc3 );
                
                _t_store_pd( psw + vlen6, rs2 );
                _t_store_pd( psw + vlen7, rs3 );
                
                p1  += vlen2;
                pj  += vlen2;
                pm  += vlen2;
                psx += vlen2;
                pcx += vlen2;
                psw += vlen8;
                
            }
            
        break;
        
    }
    
}