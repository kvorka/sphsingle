#pragma once
#include "../../../math/cvec.h"

static inline __attribute__((always_inline))
void fwd_shf_c( const int n,
                const double *restrict wght,
                const double *restrict cosx,
                const double *restrict grid,
                      double *restrict swork )

{
    
    // Memory references
    const double *restrict px  = cosx;
    const double *restrict pw  = wght;
    const double *restrict pg0 = grid + n * vlen2 * 0;
    const double *restrict pg1 = grid + n * vlen2 * 1;
    const double *restrict pg2 = grid + n * vlen2 * 2;
    const double *restrict pg3 = grid + n * vlen2 * 3;
          double *restrict psw = swork;
    
    // Registers to be used
    __td rf1, rf2, rf3, rf4, rs0, rs1, rs2, rs3, rg0, rg1, rg2, rg3;
    
    // Main cycle
    for ( int i = 0; i < n; i++ ) {
        
        rf1 = _t_load_pd( pw + vlen0 );
        rf3 = _t_load_pd( pw + vlen1 );
        rf2 = _t_load_pd( px + vlen0 );
        rf4 = _t_load_pd( px + vlen1 );
        
        rf2 = _t_mul_pd( rf1, rf2 );
        rf4 = _t_mul_pd( rf3, rf4 );
        
        rg0 = _t_load_pd( pg0 );
        rg1 = _t_load_pd( pg1 );
        
        rs0 = _t_sub_pd( rg0, rg1 );
        rs2 = _t_add_pd( rg0, rg1 );
        
        rg2 = _t_load_pd( pg2 );
        rg3 = _t_load_pd( pg3 );
        
        rs0 = _t_mul_pd( rf1, rs0 );
        rs2 = _t_mul_pd( rf2, rs2 );
        
        rs1 = _t_sub_pd( rg2, rg3 );
        rs3 = _t_add_pd( rg2, rg3 );
        
        rs1 = _t_mul_pd( rf1, rs1 );
        rs3 = _t_mul_pd( rf2, rs3 );
        
        _t_store_pd( psw + vlen0, rs0 );
        _t_store_pd( psw + vlen1, rs1 );
        _t_store_pd( psw + vlen2, rs2 );
        _t_store_pd( psw + vlen3, rs3 );
        
        rg0 = _t_load_pd( pg0 + vlen1 );
        rg1 = _t_load_pd( pg1 + vlen1 );
        
        rs0 = _t_sub_pd( rg0, rg1 );
        rs2 = _t_add_pd( rg0, rg1 );
        
        rg2 = _t_load_pd( pg2 + vlen1 );
        rg3 = _t_load_pd( pg3 + vlen1 );
        
        rs0 = _t_mul_pd( rf3, rs0 );
        rs2 = _t_mul_pd( rf4, rs2 );
        
        rs1 = _t_sub_pd( rg2, rg3 );
        rs3 = _t_add_pd( rg2, rg3 );
        
        rs1 = _t_mul_pd( rf3, rs1 );
        rs3 = _t_mul_pd( rf4, rs3 );
        
        _t_store_pd( psw + vlen4, rs0 );
        _t_store_pd( psw + vlen5, rs1 );
        _t_store_pd( psw + vlen6, rs2 );
        _t_store_pd( psw + vlen7, rs3 );
        
        px  += vlen2;
        pw  += vlen2;
        pg0 += vlen2;
        pg1 += vlen2;
        pg2 += vlen2;
        pg3 += vlen2;
        psw += vlen8;
        
    }
    
}