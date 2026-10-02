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
        
        addim = x4 + x3;
        subim = x4 - x3;
        
        x1 = ( addre + subre * t2 + addim * t1 ) / 2;
        x2 = ( subim - addim * t2 + subre * t1 ) / 2;
        
        x11[i] = x1;
        x12[i] = x2;
        
        x3 = -x1 + addre;
        x4 = +x2 - subim;
        
        x21[i] = x3;
        x22[i] = x4;
        
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