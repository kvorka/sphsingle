#include "../../../math/cvec.h"

extern inline void bwd_sum_m_c( const int n1,
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
    
    double *pmj2;
    
    bwd_set_c( m, n1, fmj+1, cosx, sinx, cc, pmm, pmj1, pmj, swork );
    
    for ( int ima = 1; ima <= nma; ima++ ) {
        
        pmj2 = pmj1;
        pmj1 = pmj;
        pmj  = pmj2;
        
        bwd_rec_c( n1, cc+4*ima, fmj+2*ima, cosx2, pmj1, pmj, swork );
        
    }
    
    bwd_shf_sub( n1, cosx, swork, grid );
    
}