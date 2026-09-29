#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void fwd_rec_c( const int n,
                const double *restrict swork,
                const double *restrict fmj,
                const double *restrict cosx2,
                const double *restrict pmj1,
                      double *restrict pmj,
                      double *restrict cr )

{
    
    // Constants
    const __td rf0 = _t_set1_pd( *( fmj + 0 ) );
    const __td rf1 = _t_set1_pd( *( fmj + 1 ) );
    
    // Memory references
    const double *restrict psw = swork;
    const double *restrict p1  = pmj1;
    const double *restrict px2 = cosx2;
          double *restrict pj  = pmj;
    
    // Accumulators
    __td rc0 = _t_setzero_pd();
    __td rc1 = _t_setzero_pd();
    __td rc2 = _t_setzero_pd();
    __td rc3 = _t_setzero_pd();
    
    // Registers to be used
    __td rx0, rx1, rp0, rp1, rs0, rs1, rs2, rs3;
    __m256d reg0, reg1, reg2, reg3;
    
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
        
        rc0 = _t_add_pd( rs0, rc0 );
        rc1 = _t_add_pd( rs1, rc1 );
        
        rs2 = _t_load_pd( psw + vlen2 );
        rs3 = _t_load_pd( psw + vlen3 );
        
        rs2 = _t_mul_pd( rp0, rs2 );
        rs3 = _t_mul_pd( rp0, rs3 );
        
        rc2 = _t_add_pd( rs2, rc2 );
        rc3 = _t_add_pd( rs3, rc3 );
        
        rs0 = _t_load_pd( psw + vlen4 );
        rs1 = _t_load_pd( psw + vlen5 );
        
        rs0 = _t_mul_pd( rp1, rs0 );
        rs1 = _t_mul_pd( rp1, rs1 );
        
        rc0 = _t_add_pd( rs0, rc0 );
        rc1 = _t_add_pd( rs1, rc1 );
        
        rs2 = _t_load_pd( psw + vlen6 );
        rs3 = _t_load_pd( psw + vlen7 );
        
        rs2 = _t_mul_pd( rp1, rs2 );
        rs3 = _t_mul_pd( rp1, rs3 );
        
        rc2 = _t_add_pd( rs2, rc2 );
        rc3 = _t_add_pd( rs3, rc3 );
        #endif
        
        pj  += vlen2;
        p1  += vlen2;
        px2 += vlen2;
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