#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void bwd_c2r_c( const int jmax,
                const double complex *restrict cjm,
                      double *restrict cab )

{
    
    // Casting memory references
    const double *restrict prjm = ( const double * ) cjm;
    
    // Memory references to be used
    const double *pr;
          double *pcab;
    
    // Indexing variables
    int j, idx, step;
    
    // Registers to be used
    __m128d s1, s2, s3, s4, s5, s6, s7, s8;
    
    // Cycle over orders
    for ( int m = 0; m <= jmax; m++ ) {
        
        // Stepping
        idx  = m * ( m + 1 );
        step = 2 * ( m + 1 );
        
        // Base pointers
        pcab = cab  - idx + ( step - 2 ) * ( jmax + 2 );
        pr   = prjm + idx + ( step - 2 );
        
        // Cycle over degrees
        for ( j = m; j <= jmax-8; j += 8 ) {
            
            s1 = _mm_loadu_pd( pr + 0 * step );
            s2 = _mm_loadu_pd( pr + 1 * step );
            
            _mm_storeu_pd( pcab + 0, s1 );
            _mm_storeu_pd( pcab + 2, s2 );
            
            s3 = _mm_loadu_pd( pr + 2 * step + 2 );
            s4 = _mm_loadu_pd( pr + 3 * step + 6 );
            
            _mm_storeu_pd( pcab + 4, s3 );
            _mm_storeu_pd( pcab + 6, s4 );
            
            s5 = _mm_loadu_pd( pr + 4 * step + 12 );
            s6 = _mm_loadu_pd( pr + 5 * step + 20 );
            
            _mm_storeu_pd( pcab +  8, s5 );
            _mm_storeu_pd( pcab + 10, s6 );
            
            s7 = _mm_loadu_pd( pr + 6 * step + 30 );
            s8 = _mm_loadu_pd( pr + 7 * step + 42 );
            
            _mm_storeu_pd( pcab + 12, s7 );
            _mm_storeu_pd( pcab + 14, s8 );
            
            pr   +=  8 * step + 56;
            pcab += 16;
            step += 16;
            
        }
        
        // Remainder degree cycle
        for ( ; j <= jmax; j++ ) {
            
            _mm_storeu_pd( pcab, _mm_loadu_pd( pr ) );
            
            pr   += step;
            pcab += 2;
            step += 2;
            
        }
        
    }
    
}

extern inline __attribute__((always_inline))
void bwd_rxd_c( const int jmax,
                const double *restrict emj,
                const double *restrict amj,
                const double *restrict icab,
                      double *restrict ocab )

{
    
    // Memory references
    const double *pemj = emj-2;
    const double *pamj = amj-1;
    const double *pci  = icab-2;
          double *pco  = ocab-4;
    
    // Registers to be used
    const __m128d s00 = _mm_setzero_pd();
          __m128d se1, se2, sa1, sc1, sc2, sc3;
    
    // Cycle over orders
    for ( int m = 0; m <= jmax-1; m++ ) {
        
        // j == m
        pci  += 2;
        pco  += 4;
        pemj += 2;
        pamj += 1;
        
        se1 = _mm_set1_pd( *( pemj + 1 ) );
        sa1 = _mm_set1_pd( *( pamj     ) );
        
        sc2 = _mm_loadu_pd( pci + 0 );
        sc1 = _mm_loadu_pd( pci + 2 );
        
        se1 = _mm_mul_pd( sa1, se1 );
        
        sc2 = _mm_mul_pd( sa1, sc2 );
        sc3 = _mm_mul_pd( se1, sc1 );
        
        _mm_storeu_pd( pco + 0, sc3 );
        _mm_storeu_pd( pco + 2, sc2 );
        
        // Cycle over j > m
        #pragma omp unroll partial (2)
        for ( int j = 1; j <= ( jmax-1-m ) / 2; j++ ) {
            
            pci  += 4;
            pco  += 4;
            pemj += 2;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj + 0 ) );
            se2 = _mm_set1_pd( *( pemj + 1 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            
            sc2 = _mm_loadu_pd( pci + 0 );
            sc3 = _mm_loadu_pd( pci + 2 );
            
            sc1 = _mm_fmadd_pd( se2, sc3, sc1 );
            sc2 = _mm_mul_pd( sa1, sc2 );
            
            sc1 = _mm_mul_pd( sa1, sc1 );
            
            _mm_storeu_pd( pco + 0, sc1 );
            _mm_storeu_pd( pco + 2, sc2 );
            
            sc1 = sc3;
            
        }
        
        // j == jmax
        if ( ( jmax - m ) % 2 == 0 ) {
            
            pci  += 4;
            pco  += 4;
            pemj += 2;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj + 0 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );
            
            se1 = _mm_mul_pd( sa1, se1 );
            sc2 = _mm_loadu_pd( pci + 0 );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_mul_pd( sa1, sc2 );
            
            _mm_storeu_pd( pco + 0, sc1 );
            _mm_storeu_pd( pco + 2, sc2 );
            
        } else {
            
            pci  += 2;
            pco  += 4;
            pemj += 1;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj + 1 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );
            
            se1 = _mm_mul_pd( sa1, se1 );
            sc2 = _mm_loadu_pd( pci + 0 );
            
            sc2 = _mm_mul_pd( se1, sc2 );
            
            _mm_storeu_pd( pco + 0, sc2 );
            _mm_storeu_pd( pco + 2, s00 );
            
        }
        
    }
    
    // m == jmax
    sa1 = _mm_set1_pd( *( pamj+1 ) );
    sc1 = _mm_loadu_pd( pci + 2 );
    
    sc1 = _mm_mul_pd( sa1, sc1 );
    
    _mm_storeu_pd( pco + 4, s00 );
    _mm_storeu_pd( pco + 6, sc1 );
    
}