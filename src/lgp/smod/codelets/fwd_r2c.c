#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void fwd_rxd_c( const int jmax,
                const double *restrict emj,
                const double *restrict icab,
                      double *restrict ocab )

{
    
    // Memory references
    const double *pemj = emj  - 2;
    const double *pci  = icab - 4;
          double *pco  = ocab - 2;
    
    // Registers to be used
    __m128d se1, se2, sc1, sc2, sc3;
    
    // Cycle over orders
    for ( int m = 0; m <= jmax - 1; m++ ) {
        
        // j == m
        pci  += 4;
        pco  += 2;
        pemj += 2;
        
        sc1 = _mm_loadu_pd( pci + 0 );
        sc2 = _mm_loadu_pd( pci + 2 );
        
        _mm_storeu_pd( pco, sc2 );
        
        // Cycle over j > m
        #pragma omp unroll partial (2)
        for ( int j = 1; j <= ( jmax - 1 - m ) / 2; j++ ) {
            
            pci  += 4;
            pco  += 4;
            pemj += 2;
            
            se1 = _mm_set1_pd( *( pemj - 1 ) );
            se2 = _mm_set1_pd( *( pemj + 0 ) );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_loadu_pd( pci ); 
            
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            sc3 = _mm_loadu_pd( pci + 2 );
            
            _mm_storeu_pd( pco - 2, sc1 );
            _mm_storeu_pd( pco + 0, sc3 );
            
            sc1 = sc2;
            
        }
        
        // j == jmax
        if ( ( jmax - m ) % 2 == 0 ) {
            
            pci  += 4;
            pco  += 4;
            pemj += 2;
            
            se1 = _mm_set1_pd( *( pemj - 1 ) );
            se2 = _mm_set1_pd( *( pemj + 0 ) );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_loadu_pd( pci );
            
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            sc3 = _mm_loadu_pd( pci + 2 );
            
            _mm_storeu_pd( pco - 2, sc1 );
            _mm_storeu_pd( pco + 0, sc3 );
            
        } else {
            
            pci  += 4;
            pco  += 2;
            pemj += 1;
            
            se1 = _mm_set1_pd( *( pemj + 0 ) );
            se2 = _mm_set1_pd( *( pemj + 1 ) );
            
            sc1 = _mm_mul_pd( se1, sc1 );
            sc2 = _mm_loadu_pd( pci );
            
            sc1 = _mm_fmadd_pd( se2, sc2, sc1 );
            
            _mm_storeu_pd( pco, sc1 );
            
        }
        
    }
    
    // m == jmax
    _mm_storeu_pd( pco + 2, _mm_loadu_pd( pci + 6 ) );
    
}

extern inline __attribute__((always_inline))
void fwd_r2c_c( const int jmax,
                const double *restrict cab,
                      double complex *restrict cjm )
{
    // Casting memory references
    double *restrict prjm = ( double * ) cjm;
    
    // Memory references to be used
    const double *pcab;
          double *pr;
    
    // Indexing variables
    int m = 0, j, idx, step;
    
    // Registers to be used
    __m128d s1, s2, s3, s4, s5, s6, s7, s8;
    
    // Order zero
    {
        
        // Stepping
        idx  = m * ( m + 1 );
        step = 2 * ( m + 1 );
        
        // Base pointers
        pcab = cab;
        pr   = prjm + ( step - 2 ) + idx;
        
        // Cycle over degrees
        for ( j = m; j <= jmax - 8; j += 8 ) {
            
            s1 = _mm_set_pd( 0.0, pcab[0]  );
            s2 = _mm_set_pd( 0.0, pcab[2]  );
            
            _mm_storeu_pd( pr + 0 * step, s1 );
            _mm_storeu_pd( pr + 1 * step, s2 );
            
            s3 = _mm_set_pd( 0.0, pcab[4]  );
            s4 = _mm_set_pd( 0.0, pcab[6]  );
            
            _mm_storeu_pd( pr + 2 * step + 2, s3 );
            _mm_storeu_pd( pr + 3 * step + 6, s4 );
            
            s5 = _mm_set_pd( 0.0, pcab[8]  );
            s6 = _mm_set_pd( 0.0, pcab[10] );
            
            _mm_storeu_pd( pr + 4 * step + 12, s5 );
            _mm_storeu_pd( pr + 5 * step + 20, s6 );
            
            s7 = _mm_set_pd( 0.0, pcab[12]  );
            s8 = _mm_set_pd( 0.0, pcab[14] );
            
            _mm_storeu_pd( pr + 6 * step + 30, s7 );
            _mm_storeu_pd( pr + 7 * step + 42, s8 );
            
            pr   += 8 * step + 56;
            pcab += 16;
            step += 16;
            
        }
        
        // Remainder cycle 
        for ( ; j <= jmax; j++ ) {
            
            _mm_storeu_pd( pr, _mm_set_pd( 0.0, *pcab ) );
            
            pr   += step;
            pcab += 2;
            step += 2;
            
        }
        
    }
    
    // Non-zero orders
    for ( m = 1; m <= jmax; m++ ) {
        
        // Stepping
        idx  = m * ( m + 1 );
        step = 2 * ( m + 1 );
        
        // Base pointers
        pcab = cab  - idx + ( step - 2 ) * ( jmax + 2 ); 
        pr   = prjm + idx + ( step - 2 );
        
        // Cycle over degrees
        for ( j = m; j <= jmax - 8; j += 8 ) {
            
            s1 = _mm_loadu_pd( pcab + 0 );
            s2 = _mm_loadu_pd( pcab + 2 );
            
            _mm_storeu_pd( pr + 0 * step, s1 );
            _mm_storeu_pd( pr + 1 * step, s2 );
            
            s3 = _mm_loadu_pd( pcab + 4 );
            s4 = _mm_loadu_pd( pcab + 6 );
            
            _mm_storeu_pd( pr + 2 * step + 2, s3 );
            _mm_storeu_pd( pr + 3 * step + 6, s4 );
            
            s5 = _mm_loadu_pd( pcab +  8 );
            s6 = _mm_loadu_pd( pcab + 10 );
            
            _mm_storeu_pd( pr + 4 * step + 12, s5 );
            _mm_storeu_pd( pr + 5 * step + 20, s6 );
            
            s7 = _mm_loadu_pd( pcab + 12 );
            s8 = _mm_loadu_pd( pcab + 14 );
            
            _mm_storeu_pd( pr + 6 * step + 30, s7 );
            _mm_storeu_pd( pr + 7 * step + 42, s8 );
            
            pr   +=  8 * step + 56;
            pcab += 16;
            step += 16;
            
        }
        
        // Remainder degree cycle
        for ( ; j <= jmax; j++ ) {
            
            _mm_storeu_pd( pr, _mm_loadu_pd( pcab ) );
            
            pr   += step;
            pcab += 2;
            step += 2;
            
        }
        
    }
    
}