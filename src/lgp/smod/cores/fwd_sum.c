#include "../codelets/fwd_set.h"
#include "../codelets/fwd_rec.h"
#include "../codelets/fwd_shf.h"

extern inline __attribute__((always_inline))
void fwd_sum_m_c( const int n1,
                  const int m,
                  const int nma,
                  const double *restrict fmj,
                  const double *restrict cosx,
                  const double *restrict sinx,
                  const double *restrict cosx2,
                  const double *restrict wght,
                        double *restrict pmm,
                        double *restrict pmj1,
                        double *restrict pmj,
                        double *restrict swork,
                        double *restrict cr,
                        double *restrict grid )

{
    
    // Memory references
    const double *restrict pfmj = fmj;
          double *restrict pcr  = cr;
    
    // Iterator
    int ima;
    
    /* After the FFT, we need to shuffle the packing north/south 
    and real/imaginary into packing suitable for summation. */
    fwd_shf_c( n1, wght, cosx, grid, swork );
    
    /* Starting from degree j equal to order m, we need to forward the recursion 
    for pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need
    to set the initial value of swork to cc * pmj (first member of the sum). */
    fwd_set_c( m, n1, pfmj+1, cosx, sinx, swork, pmm, pmj1, pmj, pcr );
    
    pfmj += 2;
    pcr  += 4;
    
    /* Following with the recursion for degrees m+1 to jmax. We need to repointer our
    polynomials, follow with recursion and add cc * pmj to our swork accumulator. */
    for ( ima = 1; ima <= nma-1; ima += 2 ) {
        
        fwd_rec_c( n1, swork, pfmj+0, cosx2, pmj,  pmj1, pcr+0 );
        fwd_rec_c( n1, swork, pfmj+2, cosx2, pmj1, pmj,  pcr+4 );
        
        pfmj += 4;
        pcr  += 8;
        
    }
    
    if ( ima == nma ) { fwd_rec_c( n1, swork, pfmj, cosx2, pmj,  pmj1, pcr ); }
    
}

extern inline __attribute__((always_inline))
void fwd_sum_jmax_c( const int n1,
                     const double *restrict fmj,
                     const double *restrict cosx,
                     const double *restrict sinx,
                     const double *restrict wght,
                           double *restrict pmm,
                           double *restrict pmj1,
                           double *restrict pmj,
                           double *restrict swork,
                           double *restrict cr,
                           double *restrict grid )

{
    
    /* After the FFT, we need to shuffle the packing north/south 
    and real/imaginary into packing suitable for summation. */
    fwd_shf_c( n1, wght, cosx, grid, swork );
    
    /* Starting from degree j equal to order m, we need to forward the recursion 
    for pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need
    to set the initial value of swork to cc * pmj (first member of the sum). */
    fwd_set_c( 1, n1, fmj, cosx, sinx, swork, pmm, pmj1, pmj, cr );
    
}