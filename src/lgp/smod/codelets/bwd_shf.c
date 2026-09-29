#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void bwd_shf_c( const int n,
                const double *restrict cosx,
                const double *restrict swork,
                      double *restrict grid )

{
    
    // Memory references
    const double *restrict px  = cosx;
    const double *restrict psw = swork;
          double *restrict pg0 = grid + n * vlen2 * 0;
          double *restrict pg1 = grid + n * vlen2 * 1;
          double *restrict pg2 = grid + n * vlen2 * 2;
          double *restrict pg3 = grid + n * vlen2 * 3;
    
    // Registers to be used
    __td rc0, rc1, rs0, rs1, rs2, rs3, rg0, rg1, rg2, rg3;
    
    // Main cycle
    for ( int i = 0; i < n; i++ ) {
        
        rc0 = _t_load_pd( px + vlen0 );
        rc1 = _t_load_pd( px + vlen1 );
        
        #if defined (__FMA__)
        rs0 = _t_load_pd( psw + vlen0 );
        rs1 = _t_load_pd( psw + vlen1 );
        rs2 = _t_load_pd( psw + vlen2 );
        rs3 = _t_load_pd( psw + vlen3 );
        
        rg0 = _t_fmadd_pd( rs2, rc0, rs0 );
        rg1 = _t_fmsub_pd( rs2, rc0, rs0 );
        
        _t_store_pd( pg0, rg0 );
        _t_store_pd( pg1, rg1 );
        
        rg2 = _t_fmadd_pd( rs3, rc0, rs1 );
        rg3 = _t_fmsub_pd( rs3, rc0, rs1 );
        
        _t_store_pd( pg2, rg2 );
        _t_store_pd( pg3, rg3 );
        
        rs0 = _t_load_pd( psw + vlen4 );
        rs1 = _t_load_pd( psw + vlen5 );
        rs2 = _t_load_pd( psw + vlen6 );
        rs3 = _t_load_pd( psw + vlen7 );
        
        rg0 = _t_fmadd_pd( rs2, rc1, rs0 );
        rg1 = _t_fmsub_pd( rs2, rc1, rs0 );
        
        _t_store_pd( pg0 + vlen1, rg0 );
        _t_store_pd( pg1 + vlen1, rg1 );
        
        rg2 = _t_fmadd_pd( rs3, rc1, rs1 );
        rg3 = _t_fmsub_pd( rs3, rc1, rs1 );
        
        _t_store_pd( pg2 + vlen1, rg2 );
        _t_store_pd( pg3 + vlen1, rg3 );
        #else
        rs2 = _t_load_pd( psw + vlen2 );
        rs3 = _t_load_pd( psw + vlen3 );
        
        rg0 = _t_mul_pd( rs2, rc0 );
        rg1 = _t_mul_pd( rs2, rc0 );
        
        rs0 = _t_load_pd( psw + vlen0 );
        rs1 = _t_load_pd( psw + vlen1 );
        
        rg2 = _t_mul_pd( rs3, rc0 );
        rg3 = _t_mul_pd( rs3, rc0 );
        
        rg0 = _t_add_pd( rg0, rs0 );
        rg1 = _t_sub_pd( rg1, rs0 );
        
        _t_store_pd( pg0, rg0 );
        _t_store_pd( pg1, rg1 );
        
        rg2 = _t_add_pd( rg2, rs1 );
        rg3 = _t_sub_pd( rg3, rs1 );
        
        _t_store_pd( pg2, rg2 );
        _t_store_pd( pg3, rg3 );
        
        rs2 = _t_load_pd( psw + vlen6 );
        rs3 = _t_load_pd( psw + vlen7 );
        
        rg0 = _t_mul_pd( rs2, rc1 );
        rg1 = _t_mul_pd( rs2, rc1 );
        
        rs0 = _t_load_pd( psw + vlen4 );
        rs1 = _t_load_pd( psw + vlen5 );
        
        rg2 = _t_mul_pd( rs3, rc1 );
        rg3 = _t_mul_pd( rs3, rc1 );
        
        rg0 = _t_add_pd( rg0, rs0 );
        rg1 = _t_sub_pd( rg1, rs0 );
        
        _t_store_pd( pg0 + vlen1, rg0 );
        _t_store_pd( pg1 + vlen1, rg1 );
        
        rg2 = _t_add_pd( rg2, rs1 );
        rg3 = _t_sub_pd( rg3, rs1 );
        
        _t_store_pd( pg2 + vlen1, rg2 );
        _t_store_pd( pg3 + vlen1, rg3 );
        #endif
        
        px  += vlen2;
        pg0 += vlen2;
        pg1 += vlen2;
        pg2 += vlen2;
        pg3 += vlen2;
        psw += vlen8;
        
    }
    
}