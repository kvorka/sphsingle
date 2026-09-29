#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void fwd_set_c( const int n,
                const double *restrict swork,
                const double *restrict pmm,
                      double *restrict pmj1,
                      double *restrict pmj,
                      double *restrict cr )

{
    
    // Constants
    const __td r00 = _t_setzero_pd();
    
    // Memory references
    const double *restrict psw = swork;
    const double *restrict pm  = pmm;
          double *restrict p1  = pmj1;
          double *restrict pj  = pmj;
    
    // Accumulators
    __td rc0 = _t_setzero_pd();
    __td rc1 = _t_setzero_pd();
    __td rc2 = _t_setzero_pd();
    __td rc3 = _t_setzero_pd();
    
    // Registers to be used
    __td rp0, rp1, rs0, rs1, rs2, rs3;
    __m256d reg0, reg1, reg2, reg3;
    
    // Main cycle
    for ( int i = 0; i < n; i++ ) {
        
        _t_store_pd( p1 + vlen0, r00 );
        _t_store_pd( p1 + vlen1, r00 );
        
        rp0 = _t_load_pd( pm + vlen0 );
        rp1 = _t_load_pd( pm + vlen1 );
        
        _t_store_pd( pj + vlen0, rp0 );
        _t_store_pd( pj + vlen1, rp1 );
        
        #if defined (__FMA__)
        rs0 = _t_load_pd( psw + vlen0 );
        rs1 = _t_load_pd( psw + vlen1 );
        
        rc0 = _t_fmadd_pd( rp0, rs0, rc0 );
        rc1 = _t_fmadd_pd( rp0, rs1, rc1 );
        
        rs2 = _t_load_pd( psw + vlen2 );
        rs3 = _t_load_pd( psw + vlen3 );
        
        rc2 = _t_fmadd_pd( rp0, rs2, rc2 );
        rc3 = _t_fmadd_pd( rp0, rs3, rc3 );
        
        rs0 = _t_load_pd( psw + vlen4 );
        rs1 = _t_load_pd( psw + vlen5 );
        
        rc0 = _t_fmadd_pd( rp1, rs0, rc0 );
        rc1 = _t_fmadd_pd( rp1, rs1, rc1 );
        
        rs2 = _t_load_pd( psw + vlen6 );
        rs3 = _t_load_pd( psw + vlen7 );
        
        rc2 = _t_fmadd_pd( rp1, rs2, rc2 );
        rc3 = _t_fmadd_pd( rp1, rs3, rc3 );
        #else
        rs0 = _t_load_pd( psw + vlen0 );
        rs1 = _t_load_pd( psw + vlen1 );
        
        rs0 = _t_mul_pd( rp0, rs0 );
        rs1 = _t_mul_pd( rp0, rs1 );
        
        rs2 = _t_load_pd( psw + vlen2 );
        rs3 = _t_load_pd( psw + vlen3 );
        
        rc0 = _t_add_pd( rs0, rc0 );
        rc1 = _t_add_pd( rs1, rc1 );
        
        rs2 = _t_mul_pd( rp0, rs2 );
        rs3 = _t_mul_pd( rp0, rs3 );
        
        rc2 = _t_add_pd( rs2, rc2 );
        rc3 = _t_add_pd( rs3, rc3 );
        
        rs0 = _t_load_pd( psw + vlen4 );
        rs1 = _t_load_pd( psw + vlen5 );
        
        rs0 = _t_mul_pd( rp1, rs0 );
        rs1 = _t_mul_pd( rp1, rs1 );
        
        rs2 = _t_load_pd( psw + vlen6 );
        rs3 = _t_load_pd( psw + vlen7 );
        
        rc0 = _t_add_pd( rs0, rc0 );
        rc1 = _t_add_pd( rs1, rc1 );
        
        rs2 = _t_mul_pd( rp1, rs2 );
        rs3 = _t_mul_pd( rp1, rs3 );
        
        rc2 = _t_add_pd( rs2, rc2 );
        rc3 = _t_add_pd( rs3, rc3 );
        #endif
        
        p1  += vlen2;
        pj  += vlen2;
        pm  += vlen2;
        psw += vlen8;
        
    }
    
    // Horizontal sum into cr
    rs0 = _t_unpacklo_pd( rc0, rc1 );
    rs1 = _t_unpackhi_pd( rc0, rc1 );
    rs2 = _t_unpacklo_pd( rc2, rc3 );
    rs3 = _t_unpackhi_pd( rc2, rc3 );
    
    rs0 = _t_add_pd( rs0, rs1 );
    rs2 = _t_add_pd( rs2, rs3 );
    
    #if !defined (__AVX512F__)
    reg0 = rs0;
    reg2 = rs2;
    #else
    reg1 = _mm512_extractf64x4_pd( rs0, 1 );
    reg3 = _mm512_extractf64x4_pd( rs2, 1 );
    
    reg0 = _mm256_add_pd( _mm512_castpd512_pd256( rs0 ), reg1 );
    reg2 = _mm256_add_pd( _mm512_castpd512_pd256( rs2 ), reg3 );
    #endif
    
    reg1 = _mm256_permute2f128_pd( reg0, reg2, 0x31 );
    reg3 = _mm256_permute2f128_pd( reg0, reg2, 0x20 );
    
    reg0 = _mm256_add_pd( reg1, reg3 );
    
    _mm256_storeu_pd( cr, reg0 );
    
}