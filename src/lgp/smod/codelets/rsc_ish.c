#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void is_rescale_c( const int n,
                   const double *restrict amj,
                         double *restrict rcab )

{
    
    // Memory references
    const double *restrict pamj = amj;
          double *restrict pcab = rcab;
    
    // Registers to be used
    __m256d rc0, rc1, rc2, rc3,
            ra0, ra1, ra2, ra3;
    
    // Iterator
    int i = 0;
    
    // Main cycle
    for ( ; i <= n-4; i += 4 ) {
        
        ra0 = _mm256_set1_pd( *( pamj + 0 ) );
        ra1 = _mm256_set1_pd( *( pamj + 1 ) );
        
        rc0 = _mm256_loadu_pd( pcab + 0 );
        rc1 = _mm256_loadu_pd( pcab + 4 );
        
        ra2 = _mm256_set1_pd( *( pamj + 2 ) );
        ra3 = _mm256_set1_pd( *( pamj + 3 ) );
        
        rc0 = _mm256_mul_pd( ra0, rc0 );
        rc1 = _mm256_mul_pd( ra1, rc1 );
        
        rc2 = _mm256_loadu_pd( pcab +  8 );
        rc3 = _mm256_loadu_pd( pcab + 12 );
        
        _mm256_storeu_pd( pcab + 0, rc0 );
        _mm256_storeu_pd( pcab + 4, rc1 );
        
        rc2 = _mm256_mul_pd( ra2, rc2 );
        rc3 = _mm256_mul_pd( ra3, rc3 );
        
        _mm256_storeu_pd( pcab +  8, rc2 );
        _mm256_storeu_pd( pcab + 12, rc3 );
        
        pamj +=  4;
        pcab += 16;
        
    }
    
    // Remainder cycle
    for ( ; i < n; i++ ) {
        
        ra0 = _mm256_set1_pd( *pamj );
        rc0 = _mm256_loadu_pd( pcab );
        
        rc0 = _mm256_mul_pd( ra0, rc0 );
        
        _mm256_storeu_pd( pcab, rc0 );
        
        pamj += 1;
        pcab += 4;
        
    }
    
}