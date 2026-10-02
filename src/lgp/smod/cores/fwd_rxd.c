#include "../../../math/cvec.h"

extern inline void fwd_rxd_c( const int jmax,
                              const double *restrict emj,
                              const double *restrict amj,
                              const double *restrict icab,
                                    double complex *restrict cjm )
{
    // Pointer to double representation of complex array
    double *restrict prjm = ( double * ) cjm;

    // Memory references
    const double *pemj = emj  - 2;
    const double *pamj = amj  - 1;
    const double *pci  = icab - 4;
          double *pr;
    
    // Indexing variables
    int step;
    
    // Registers to be used
    __m128d se1, se2, sa1, sc1, sc2, sc3, scT;

    // Order m == 0
    {
        
        // Zero register
        const __m128d s00 = _mm_setzero_pd();
        
        // Stepping and base pointer setup
        step = 2;
        pr   = prjm;
        
        // j == 0
        pci  += 4;
        pemj += 2;
        pamj += 1;
        
        sa1 = _mm_set1_pd( *( pamj ) );
        
        sc1 = _mm_loadu_pd( pci + 0 );
        sc2 = _mm_loadu_pd( pci + 2 );
        
        sc1 = _mm_mul_pd( sa1, sc1 );
        sc2 = _mm_mul_pd( sa1, sc2 );
        
        sc2 = _mm_unpacklo_pd( sc2, s00 );
        
        _mm_storeu_pd( pr, sc2 );
        
        pr   += step;
        step += 2;

        // Cycle over j > 0
        #pragma omp unroll partial (2)
        for ( int j = 1; j <= ( jmax - 1 ) / 2; j++ ) {
            
            pci  += 4;
            pemj += 2;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj - 1 ) );
            se2 = _mm_set1_pd( *( pemj + 0 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_loadu_pd( pci + 0 );
            sc3 = _mm_loadu_pd( pci + 2 );
            
            sc2 = _mm_mul_pd( sa1, sc2 );
            sc3 = _mm_mul_pd( sa1, sc3 );
            
            #if defined (__FMA__)
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            #else
            scT = _mm_mul_pd( se2, sc2 );
            sc1 = _mm_add_pd( sc1, scT );
            #endif
            
            sc1 = _mm_unpacklo_pd( sc1, s00 );
            sc3 = _mm_unpacklo_pd( sc3, s00 );
            
            _mm_storeu_pd( pr,        sc1 );
            _mm_storeu_pd( pr + step, sc3 );
            
            pr   += 2 * step + 2;
            step += 4;

            sc1 = sc2;
            
        }

        // j == jmax
        if ( jmax % 2 == 0 ) {
            
            pci  += 4;
            pemj += 2;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj - 1 ) );
            se2 = _mm_set1_pd( *( pemj + 0 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );
            
            sc2 = _mm_loadu_pd( pci + 0 );
            sc3 = _mm_loadu_pd( pci + 2 );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_mul_pd( sa1, sc2 );
            
            #if defined (__FMA__)
            sc3 = _mm_mul_pd( sa1, sc3 );
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            #else
            sc3 = _mm_mul_pd( sa1, sc3 );
            scT = _mm_mul_pd( se2, sc2 );
            
            sc1 = _mm_add_pd( sc1, scT );
            #endif
            
            sc3 = _mm_unpacklo_pd( sc3, s00 );
            sc1 = _mm_unpacklo_pd( sc1, s00 );
            
            _mm_storeu_pd( pr,        sc1 );
            _mm_storeu_pd( pr + step, sc3 );
            
        } else {
            
            pci  += 4;
            pemj += 1;
            pamj += 1;
            
            se1 = _mm_set1_pd( *( pemj + 0 ) );
            se2 = _mm_set1_pd( *( pemj + 1 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );
            
            sc2 = _mm_loadu_pd( pci );
            sc1 = _mm_mul_pd( se1, sc1 );
            
            #if defined (__FMA__)
            sc2 = _mm_mul_pd( sa1, sc2 );
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            #else
            sc2 = _mm_mul_pd( sa1, sc2 );
            scT = _mm_mul_pd( se2, sc2 );
            
            sc1 = _mm_add_pd( sc1, scT );
            #endif
            
            sc1 = _mm_unpacklo_pd( sc1, s00 );
            
            _mm_storeu_pd( pr, sc1 );
            
        }
        
    }
    
    // Orders m > 0
    for ( int m = 1; m <= jmax - 1; m++ ) {
        
        // Stepping and base pointer setup
        step = 2 * m + 2;
        pr   = prjm + m * ( m + 3);
        
        // j == m
        pci  += 4;
        pemj += 2;
        pamj += 1;
        
        sa1 = _mm_set1_pd( *( pamj ) );
        
        sc1 = _mm_loadu_pd( pci + 0 );
        sc2 = _mm_loadu_pd( pci + 2 );
        
        sc1 = _mm_mul_pd( sa1, sc1 );
        sc2 = _mm_mul_pd( sa1, sc2 );
        
        _mm_storeu_pd( pr, sc2 );
        
        pr   += step;
        step += 2;
        
        // Cycle over j > m
        #pragma omp unroll partial (2)
        for ( int j = 1; j <= ( jmax - 1 - m ) / 2; j++ ) {

            pci  += 4;
            pemj += 2;
            pamj += 1;

            se1 = _mm_set1_pd( *( pemj - 1 ) );
            se2 = _mm_set1_pd( *( pemj + 0 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );

            sc2 = _mm_loadu_pd( pci + 0 );
            sc3 = _mm_loadu_pd( pci + 2 );

            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_mul_pd( sa1, sc2 );

            #if defined (__FMA__)
            sc3 = _mm_mul_pd( sa1, sc3 );
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            #else
            sc3 = _mm_mul_pd( sa1, sc3 );
            scT = _mm_mul_pd( se2, sc2 );
            
            sc1 = _mm_add_pd( sc1, scT );
            #endif

            _mm_storeu_pd( pr,        sc1 );
            _mm_storeu_pd( pr + step, sc3 );

            sc1 = sc2;
            
            pr   += 2 * step + 2;
            step += 4;

        }

        // j == jmax
        if ( ( jmax - m ) % 2 == 0 ) {

            pci  += 4;
            pemj += 2;
            pamj += 1;

            se1 = _mm_set1_pd( *( pemj - 1 ) );
            se2 = _mm_set1_pd( *( pemj + 0 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );

            sc2 = _mm_loadu_pd( pci + 0 );
            sc3 = _mm_loadu_pd( pci + 2 );

            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_mul_pd( sa1, sc2 );

            #if defined (__FMA__)
            sc3 = _mm_mul_pd( sa1, sc3 );
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            #else
            sc3 = _mm_mul_pd( sa1, sc3 );
            scT = _mm_mul_pd( se2, sc2 );
            
            sc1 = _mm_add_pd( sc1, scT );
            #endif

            _mm_storeu_pd( pr,        sc1 );
            _mm_storeu_pd( pr + step, sc3 );

        } else {

            pci  += 4;
            pemj += 1;
            pamj += 1;

            se1 = _mm_set1_pd( *( pemj + 0 ) );
            se2 = _mm_set1_pd( *( pemj + 1 ) );
            sa1 = _mm_set1_pd( *( pamj     ) );

            sc2 = _mm_loadu_pd( pci );
            sc1 = _mm_mul_pd( se1, sc1 );

            #if defined (__FMA__)
            sc2 = _mm_mul_pd( sa1, sc2 );
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            #else
            sc2 = _mm_mul_pd( sa1, sc2 );
            scT = _mm_mul_pd( se2, sc2 );
            
            sc1 = _mm_add_pd( sc1, scT );
            #endif

            _mm_storeu_pd( pr, sc1 );
            
        }
        
    }

    // Order m == jmax
    {
        
        pr = prjm + jmax * ( jmax + 3 );
        
        sa1 = _mm_set1_pd( *( pamj + 1 ) );
        sc1 = _mm_loadu_pd( pci + 6 );
        
        sc1 = _mm_mul_pd( sa1, sc1 );
        
        _mm_storeu_pd( pr, sc1 );
        
    }
    
}