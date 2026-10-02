#include "../../../math/cvec.h"

extern inline __attribute__((always_inline))
void fxrc0_c( const int m,
                    double *restrict arr1,
                    double *restrict arr2 )

{
    
    // Memory references
    double *restrict parr1 = arr1;
    double *restrict parr2 = arr2;
    
    // Registers to be used
    __td r00, r01, r02, r03, r04, r05,
         r06, r07, r08, r09, r10, r11;
    
    // Main cycle
    for ( int i = 0; i < m; i++ ) {
        
        r00 = _t_load_pd( parr1 + vlen0 );
        r01 = _t_load_pd( parr1 + vlen1 );
        
        r02 = _t_load_pd( parr2 + vlen0 );
        r03 = _t_load_pd( parr2 + vlen1 );
        
        r04 = _t_add_pd( r00, r02 );
        r05 = _t_add_pd( r01, r03 );
        
        r08 = _t_load_pd( parr1 + vlen2 );
        r09 = _t_load_pd( parr1 + vlen3 );
        
        r06 = _t_sub_pd( r00, r02 );
        r07 = _t_sub_pd( r01, r03 );
        
        r10 = _t_load_pd( parr2 + vlen2 );
        r11 = _t_load_pd( parr2 + vlen3 );
        
        _t_store_pd( parr1 + vlen0, r04 );
        _t_store_pd( parr1 + vlen1, r05 );
        
        r00 = _t_add_pd( r08, r10 );
        r01 = _t_add_pd( r09, r11 );
        
        _t_store_pd( parr2 + vlen0, r06 );
        _t_store_pd( parr2 + vlen1, r07 );
        
        r02 = _t_sub_pd( r08, r10 );
        r03 = _t_sub_pd( r09, r11 );
        
        _t_store_pd( parr1 + vlen2, r00 );
        _t_store_pd( parr1 + vlen3, r01 );
        
        parr1 += vlen4;
        
        _t_store_pd( parr2 + vlen2, r02 );
        _t_store_pd( parr2 + vlen3, r03 );
        
        parr2 += vlen4;
        
    }
    
}

extern inline __attribute__((always_inline))
void fxr2c_c( const int m,
              const double *restrict t,
                    double *restrict x11,
                    double *restrict x12,
                    double *restrict x21,
                    double *restrict x22 )

{
    
    // Memory references
    double *restrict px11 = x11;
    double *restrict px12 = x12;
    double *restrict px21 = x21;
    double *restrict px22 = x22;
    
    // Register constants
    const __td rt1 = _t_set1_pd( *( t + 0 ) );
    const __td rt2 = _t_set1_pd( *( t + 1 ) );
    const __td rhf = _t_set1_pd( 0.5 );
    
    // Registers to be used
    __td rx1, rx2, rx3, rx4, rx5, rx6, rx7, rx8, rA1re, rS1im, rA2re, rS2im;
    
    // Main cycle
    for ( int i = 0; i < 4*m; i += 2 ) {
        
        rx1 = _t_load_pd( px11 + vlen0 );
        rx2 = _t_load_pd( px21 + vlen0 );
        
        rA1re = _t_add_pd( rx1, rx2 );
        rx2   = _t_sub_pd( rx1, rx2 );
        
        rx3 = _t_load_pd( px12 + vlen0 );
        rx4 = _t_load_pd( px22 + vlen0 );
        
        rS1im = _t_sub_pd( rx4, rx3 );
        rx3   = _t_add_pd( rx4, rx3 );
        
        rx5 = _t_load_pd( px11 + vlen1 );
        rx6 = _t_load_pd( px21 + vlen1 );
        
        rA2re = _t_add_pd( rx5, rx6 );
        rx6   = _t_sub_pd( rx5, rx6 );
        
        rx7 = _t_load_pd( px12 + vlen1 );
        rx8 = _t_load_pd( px22 + vlen1 );
        
        rS2im = _t_sub_pd( rx8, rx7 );
        rx7   = _t_add_pd( rx8, rx7 );
        
        #if defined (__FMA__)
        rx1 = _t_fmadd_pd(  rt2, rx2, rA1re );
        rx4 = _t_fnmadd_pd( rt2, rx3, rS1im );
        rx5 = _t_fmadd_pd(  rt2, rx6, rA2re );
        rx8 = _t_fnmadd_pd( rt2, rx7, rS2im );
        
        rx1 = _t_fmadd_pd( rt1, rx3, rx1 );
        rx2 = _t_fmadd_pd( rt1, rx2, rx4 );
        rx5 = _t_fmadd_pd( rt1, rx7, rx5 );
        rx6 = _t_fmadd_pd( rt1, rx6, rx8 );
        #else
        rx1 = _t_mul_pd( rt2, rx2 );
        rx4 = _t_mul_pd( rt2, rx3 );
        rx5 = _t_mul_pd( rt2, rx6 );
        rx8 = _t_mul_pd( rt2, rx7 );
        
        rx3 = _t_mul_pd( rt1, rx3 );
        rx2 = _t_mul_pd( rt1, rx2 );
        rx7 = _t_mul_pd( rt1, rx7 );
        rx6 = _t_mul_pd( rt1, rx6 );
        
        rx1 = _t_add_pd( rA1re, rx1 );
        rx4 = _t_sub_pd( rS1im, rx4 );
        rx5 = _t_add_pd( rA2re, rx5 );
        rx8 = _t_sub_pd( rS2im, rx8 );
        
        rx1 = _t_add_pd( rx1, rx3 );
        rx2 = _t_add_pd( rx2, rx4 );
        rx5 = _t_add_pd( rx5, rx7 );
        rx6 = _t_add_pd( rx6, rx8 );
        #endif
        
        rx1 = _t_mul_pd( rhf, rx1 );
        rx2 = _t_mul_pd( rhf, rx2 );
        
        _t_store_pd( px11 + vlen0, rx1 );
        _t_store_pd( px12 + vlen0, rx2 );
        
        rx5 = _t_mul_pd( rhf, rx5 );
        rx6 = _t_mul_pd( rhf, rx6 );
        
        _t_store_pd( px11 + vlen1, rx5 );
        _t_store_pd( px12 + vlen1, rx6 );
        
        px11 += vlen2;
        px12 += vlen2;
        
        rx3 = _t_sub_pd( rA1re, rx1 );
        rx4 = _t_sub_pd( rx2, rS1im );
        
        _t_store_pd( px21 + vlen0, rx3 );
        _t_store_pd( px22 + vlen0, rx4 );
        
        rx7 = _t_sub_pd( rA2re, rx5 );
        rx8 = _t_sub_pd( rx6, rS2im );
        
        _t_store_pd( px21 + vlen1, rx7 );
        _t_store_pd( px22 + vlen1, rx8 );
        
        px21 += vlen2;
        px22 += vlen2;
        
    }
    
}

extern inline __attribute__((always_inline))
void fxc2r_c( const int m,
              const double *restrict t,
                    double *restrict x11,
                    double *restrict x12,
                    double *restrict x21,
                    double *restrict x22 )

{
    
    // Constants
    const double t1 = t[0];
    const double t2 = t[1];
    
    // Temporal variables
    double x1, x2, x3, x4, addre, subre, addim, subim;
    
    // Main loop
    #pragma omp unroll (vlen2) simd uniform (t1,t2) aligned (x11,x12,x21,x22:alignement)
    for ( int i = 0; i < vlen4 * m; i++ ) {
        
        x1 = x11[i];
        x2 = x21[i];
        
        addre = x1 + x2;
        subre = x1 - x2;
        
        x3 = x12[i];
        x4 = x22[i];
        
        addim = x3 + x4;
        subim = x3 - x4;
        
        x1 = addre - subre * t2 - addim * t1;
        x2 = subim - addim * t2 + subre * t1;
        
        x11[i] = x1;
        x12[i] = x2;
        
        x3 = -x1 + 2 * addre;
        x4 = +x2 - 2 * subim;
        
        x21[i] = x3;
        x22[i] = x4;
        
    }
    
}