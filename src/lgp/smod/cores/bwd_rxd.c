#include "../../../math/cvec.h"

extern inline void bwd_rxd_c( const int jmax,
                             const double complex *restrict cjm,
                             const double *restrict emj,
                             const double *restrict amj,
                                   double *restrict ocab )

{
    
    // Casting memory references
    const double *restrict prjm = ( const double * ) cjm;
    
    // Memory references
    const double *pemj = emj-2;
    const double *pamj = amj-1;
    const double *pr;
          double *pco  = ocab-4;
    
    // Indexing variables
    int step;
    
    // Registers to be used
    const __m128d s00 = _mm_setzero_pd();
          __m128d se1, se2, sa1, sc1, sc2, sc3, scT;
    
    // Cycle over orders
    for ( int m = 0; m <= jmax-1; m++ ) {
        
        // Stepping and base pointer setup
        step = 2 * m + 2;
        pr   = prjm + m * ( m + 3);
        
        // j == m
        pco  += 4;
        pemj += 2;
        pamj += 1;
        
        se1 = _mm_set1_pd( *( pemj + 1 ) );
        sa1 = _mm_set1_pd( *( pamj     ) );
        
        sc2 = _mm_loadu_pd( pr + 0 * step );
        sc1 = _mm_loadu_pd( pr + 1 * step );
        
        se1 = _mm_mul_pd( sa1, se1 );
        
        sc2 = _mm_mul_pd( sa1, sc2 );
        sc3 = _mm_mul_pd( se1, sc1 );
        
        _mm_storeu_pd( pco + 0, sc3 );
        _mm_storeu_pd( pco + 2, sc2 );
        
        pr   += 2 * step + 2;
        step += 4;
        
        // Cycle over j > m
        #pragma omp unroll partial (2)
        for ( int j = 1; j <= ( jmax-1-m ) / 2; j++ ) {
            
            pco  += 4;
            pemj += 2;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj + 0 ) );
            se2 = _mm_set1_pd( *( pemj + 1 ) );
            
            sa1 = _mm_set1_pd( *( pamj     ) );
            sc1 = _mm_mul_pd( se1, sc1 );
            
            sc2 = _mm_loadu_pd( pr        );
            sc3 = _mm_loadu_pd( pr + step );
            
            #if defined (__FMA__)
            sc1 = _mm_fmadd_pd( se2, sc3, sc1 );
            sc2 = _mm_mul_pd( sa1, sc2 );
            #else
            scT = _mm_mul_pd( se2, sc3 );
            sc2 = _mm_mul_pd( sa1, sc2 );
            
            sc1 = _mm_add_pd( sc1, scT );
            #endif
            
            sc1 = _mm_mul_pd( sa1, sc1 );
            
            _mm_storeu_pd( pco + 0, sc1 );
            _mm_storeu_pd( pco + 2, sc2 );
            
            sc1 = sc3;
            
            pr   += 2 * step + 2;
            step += 4;
            
        }
        
        // j == jmax
        if ( ( jmax - m ) % 2 == 0 ) {
            
            pco  += 4;
            pemj += 2;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj + 0 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );
            
            se1 = _mm_mul_pd( sa1, se1 );
            sc2 = _mm_loadu_pd( pr );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_mul_pd( sa1, sc2 );
            
            _mm_storeu_pd( pco + 0, sc1 );
            _mm_storeu_pd( pco + 2, sc2 );
            
        } else {
            
            pco  += 4;
            pemj += 1;
            pamj += 1;
            
            sa1 = _mm_set1_pd( *( pamj     ) );
            se1 = _mm_set1_pd( *( pemj + 1 ) );
            
            sc1 = _mm_mul_pd( sa1, sc1 );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            
            _mm_storeu_pd( pco + 0, sc1 );
            _mm_storeu_pd( pco + 2, s00 );
            
        }
        
    }
    
    // m == jmax
    pr = prjm + jmax * ( jmax + 3 );
    
    pco  += 4;
    pamj += 1;

    sa1 = _mm_set1_pd( *( pamj ) );
    sc1 = _mm_loadu_pd( pr );
    
    sc1 = _mm_mul_pd( sa1, sc1 );
    
    _mm_storeu_pd( pco + 0, s00 );
    _mm_storeu_pd( pco + 2, sc1 );
    
}