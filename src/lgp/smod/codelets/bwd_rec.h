#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void bwd_rec_c( const int n,
                const double *restrict cc,
                const double *restrict fmj,
                const double *restrict cosx2,
                const double *restrict pmj1,
                      double *restrict pmj,
                      double *restrict swork )

{
    
    // Constants
    const __td rc0 = _t_set1_pd( *( cc + 0 ) );
    const __td rc1 = _t_set1_pd( *( cc + 1 ) );
    const __td rc2 = _t_set1_pd( *( cc + 2 ) );
    const __td rc3 = _t_set1_pd( *( cc + 3 ) );
    
    const __td rf0 = _t_set1_pd( *( fmj + 0 ) );
    const __td rf1 = _t_set1_pd( *( fmj + 1 ) );
    
    // Memory references
    const double *restrict p1  = pmj1;
    const double *restrict px2 = cosx2;
          double *restrict pj  = pmj;
          double *restrict psw = swork;
    
    // Registers to be used
    __td rx0, rx1, rp0, rp1, rs0, rs1, rs2, rs3;
    
    // Main cycle
    for ( int i = 0; i < n; i++ ) {
        
        rx0 = _t_load_pd( px2 + vlen0 );
        rx1 = _t_load_pd( px2 + vlen1 );
        rp0 = _t_load_pd( pj  + vlen0 );
        rp1 = _t_load_pd( pj  + vlen1 );
        
        #if defined (__FMA__)
        rx0 = _t_fmsub_pd( rf0, rx0, rf1 );
        rx1 = _t_fmsub_pd( rf0, rx1, rf1 );
        
        rs0 = _t_load_pd( p1 + vlen0 );
        rs1 = _t_load_pd( p1 + vlen1 );
        
        rp0 = _t_fmsub_pd( rx0, rs0, rp0 );
        rp1 = _t_fmsub_pd( rx1, rs1, rp1 );
        #else
        rx0 = _t_mul_pd( rf0, rx0 );
        rx1 = _t_mul_pd( rf0, rx1 );
        
        rs0 = _t_load_pd( p1 + vlen0 );
        rs1 = _t_load_pd( p1 + vlen1 );
        
        rx0 = _t_sub_pd( rx0, rf1 );
        rx1 = _t_sub_pd( rx1, rf1 );
        
        rx0 = _t_mul_pd( rx0, rs0 );
        rx1 = _t_mul_pd( rx1, rs1 );
        
        rp0 = _t_sub_pd( rx0, rp0 );
        rp1 = _t_sub_pd( rx1, rp1 );
        #endif
        
        _t_store_pd( pj + vlen0, rp0 );
        _t_store_pd( pj + vlen1, rp1 );
        
        #if defined (__FMA__)
        rs0 = _t_load_pd( psw + vlen0 );
        rs1 = _t_load_pd( psw + vlen1 );
        
        rs0 = _t_fmadd_pd( rp0, rc0, rs0 );
        rs1 = _t_fmadd_pd( rp0, rc1, rs1 );
        
        rs2 = _t_load_pd( psw + vlen2 );
        rs3 = _t_load_pd( psw + vlen3 );
        
        _t_store_pd( psw + vlen0, rs0 );
        _t_store_pd( psw + vlen1, rs1 );
        
        rs2 = _t_fmadd_pd( rp0, rc2, rs2 );
        rs3 = _t_fmadd_pd( rp0, rc3, rs3 );
        
        rs0 = _t_load_pd( psw + vlen4 );
        rs1 = _t_load_pd( psw + vlen5 );
        
        _t_store_pd( psw + vlen2, rs2 );
        _t_store_pd( psw + vlen3, rs3 );
        
        rs0 = _t_fmadd_pd( rp1, rc0, rs0 );
        rs1 = _t_fmadd_pd( rp1, rc1, rs1 );
        
        rs2 = _t_load_pd( psw + vlen6 );
        rs3 = _t_load_pd( psw + vlen7 );
        
        _t_store_pd( psw + vlen4, rs0 );
        _t_store_pd( psw + vlen5, rs1 );
        
        rs2 = _t_fmadd_pd( rp1, rc2, rs2 );
        rs3 = _t_fmadd_pd( rp1, rc3, rs3 );
        
        _t_store_pd( psw + vlen6, rs2 );
        _t_store_pd( psw + vlen7, rs3 );
        #else
        rs0 = _t_load_pd( psw + vlen0 );
        rs1 = _t_load_pd( psw + vlen1 );
        
        rx0 = _t_mul_pd( rp0, rc0 );
        rx1 = _t_mul_pd( rp0, rc1 );
        
        rs0 = _t_add_pd( rs0, rx0 );
        rs1 = _t_add_pd( rs1, rx1 );
        
        _t_store_pd( psw + vlen0, rs0 );
        _t_store_pd( psw + vlen1, rs1 );
        
        rs2 = _t_load_pd( psw + vlen2 );
        rs3 = _t_load_pd( psw + vlen3 );
        
        rx0 = _t_mul_pd( rp0, rc2 );
        rx1 = _t_mul_pd( rp0, rc3 );
        
        rs2 = _t_add_pd( rs2, rx0 );
        rs3 = _t_add_pd( rs3, rx1 );
        
        _t_store_pd( psw + vlen2, rs2 );
        _t_store_pd( psw + vlen3, rs3 );
        
        rs0 = _t_load_pd( psw + vlen4 );
        rs1 = _t_load_pd( psw + vlen5 );
        
        rx0 = _t_mul_pd( rp1, rc0 );
        rx1 = _t_mul_pd( rp1, rc1 );
        
        rs0 = _t_add_pd( rs0, rx0 );
        rs1 = _t_add_pd( rs1, rx1 );
        
        _t_store_pd( psw + vlen4, rs0 );
        _t_store_pd( psw + vlen5, rs1 );
        
        rs2 = _t_load_pd( psw + vlen6 );
        rs3 = _t_load_pd( psw + vlen7 );
        
        rx0 = _t_mul_pd( rp1, rc2 );
        rx1 = _t_mul_pd( rp1, rc3 );
        
        rs2 = _t_add_pd( rs2, rx0 );
        rs3 = _t_add_pd( rs3, rx1 );
        
        _t_store_pd( psw + vlen6, rs2 );
        _t_store_pd( psw + vlen7, rs3 );
        #endif
        
        pj  += vlen2;
        p1  += vlen2;
        px2 += vlen2;
        psw += vlen8;
        
    }
    
}