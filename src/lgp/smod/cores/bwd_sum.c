#include "../codelets/bwd_set.h"
#include "../codelets/bwd_rec.h"
#include "../codelets/bwd_shf.h"

extern inline __attribute__((always_inline))
void bwd_sum_m_c( const int n1,
                  const int m,
                  const int nma,
                  const double *restrict fmj,
                  const double *restrict cosx,
                  const double *restrict sinx,
                  const double *restrict cosx2,
                        double *restrict pmm,
                        double *restrict pmj1,
                        double *restrict pmj,
                  const double *restrict cc,
                        double *restrict swork,
                        double *restrict grid )

{
    
    // Memory references
    const double *restrict pfmj = fmj;
    const double *restrict pcc  = cc;
    
    // Iterator
    int ima;
    
    /* Starting from degree j equal to order m, we need to forward the recursion 
    for pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need
    to set the initial value of swork to cc * pmj (first member of the sum). */
    bwd_set_c( m, n1, pfmj+1, cosx, sinx, pcc, pmm, pmj1, pmj, swork );
    
    pfmj += 2;
    pcc  += 4;
    
    /* Following with the recursion for degrees m+1 to jmax. We need to repointer our
    polynomials, follow with recursion and add cc * pmj to our swork accumulator. */
    for ( ima = 1; ima <= nma-1; ima += 2 ) {
        
        bwd_rec_c( n1, pcc+0, pfmj+0, cosx2, pmj,  pmj1, swork );
        bwd_rec_c( n1, pcc+4, pfmj+2, cosx2, pmj1, pmj,  swork );
        
        pfmj += 4;
        pcc  += 8;
        
    }
    
    if ( ima == nma ) { bwd_rec_c( n1, pcc, pfmj, cosx2, pmj, pmj1, swork ); }
    
    /* As we are done with computing the summation, we need to reshufle the data 
    from packed sum to south/north and real/imaginary parts for upcomming FFT.. */
    bwd_shf_c( n1, cosx, swork, grid );
    
}

extern inline __attribute__((always_inline))
void bwd_sum_jmax_c( const int n1,
                     const double *restrict fmj,
                     const double *restrict cosx,
                     const double *restrict sinx,
                           double *restrict pmm,
                           double *restrict pmj1,
                           double *restrict pmj,
                     const double *restrict cc,
                           double *restrict swork,
                           double *restrict grid )

{
    
    /* Starting from degree j equal to order m, we need to forward the recursion 
    for pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need
    to set the initial value of swork to cc * pmj (first member of the sum). */
    bwd_set_c( 1, n1, fmj, cosx, sinx, cc, pmm, pmj1, pmj, swork );
    
    /* As we are done with computing the summation, we need to reshufle the data 
    from packed sum to south/north and real/imaginary parts for upcomming FFT.. */
    bwd_shf_c( n1, cosx, swork, grid );
    
}