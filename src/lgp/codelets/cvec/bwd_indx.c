#include <stdlib.h>
#include <immintrin.h>

extern inline void bwd_indx_c( const int jmax,                // cutoff degree
                               const double *restrict emj,    // factors
                               const double *restrict icab,   // input to rescale
                               double *restrict ocab) {       // output
  
  __m128d rcab[4];
  
  int m, ma, mj;
  
  m = 0;
  {
    
    { // j == m
      
      ma = 0;
      mj = 0;
      
      rcab[0] = _mm_load_pd( icab   );
      rcab[3] = _mm_load_pd( icab+2 );
      
      rcab[2] = _mm_mul_pd( _mm_set1_pd( emj[mj+1] ), rcab[3] );
      
      _mm_store_pd( ocab,   rcab[2] );
      _mm_store_pd( ocab+2, rcab[0] );
    
    }
    
    for ( int j = 0; j < (jmax-1)/2; j++ ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj  ] ), rcab[3] );
      
      rcab[1] = _mm_load_pd( icab+2*mj   );
      rcab[3] = _mm_load_pd( icab+2*mj+2 );
      
      rcab[2] = _mm_mul_pd( _mm_set1_pd( emj[mj+1] ), rcab[3] );
      
      rcab[0] = _mm_add_pd( rcab[0], rcab[2] );
      
      _mm_store_pd( ocab+4*ma  , rcab[0] );
      _mm_store_pd( ocab+4*ma+2, rcab[1] );
      
    }
    
    if ( jmax%2 == 0 ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj] ), rcab[3] );
      rcab[1] = _mm_load_pd( icab+2*mj   );
      
      _mm_store_pd( ocab+4*ma  , rcab[0] );
      _mm_store_pd( ocab+4*ma+2, rcab[1] );
      
    } else {
      
      ma = ma+1;
      mj = mj+1;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj+1] ), _mm_load_pd( icab+2*mj ) );
      
      _mm_store_pd( ocab+4*ma, rcab[0] );
      
    }
    
  }
  
  for ( m = 1; m < jmax; m++ ) {
    
    { // j == m
      
      ma = ma+1;
      mj = mj+1;
      
      rcab[0] = _mm_load_pd( icab+2*mj   );
      rcab[3] = _mm_load_pd( icab+2*mj+2 );
      
      rcab[1] = _mm_mul_pd( _mm_set1_pd( emj[mj+m+1] ), rcab[3] );
      
      _mm_store_pd( ocab+4*ma,   rcab[1] );
      _mm_store_pd( ocab+4*ma+2, rcab[0] );
      
    }
    
    for ( int j = 0; j < (jmax-1-m)/2; j++ ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj+m  ] ), rcab[3] );
      
      rcab[1] = _mm_load_pd( icab+2*mj   );
      rcab[3] = _mm_load_pd( icab+2*mj+2 );
      
      rcab[2] = _mm_mul_pd( _mm_set1_pd( emj[mj+m+1] ), rcab[3] );
      
      rcab[0] = _mm_add_pd( rcab[0], rcab[2] );
      
      _mm_store_pd( ocab+4*ma  , rcab[0] );
      _mm_store_pd( ocab+4*ma+2, rcab[1] );
      
    }
    
    if ( (jmax-m)%2 == 0 ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj+m] ), rcab[3] );
      rcab[1] = _mm_load_pd( icab+2*mj   );
      
      _mm_store_pd( ocab+4*ma  , rcab[0] );
      _mm_store_pd( ocab+4*ma+2, rcab[1] );
      
    } else {
      
      ma = ma+1;
      mj = mj+1;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj+m+1] ), _mm_load_pd( icab+2*mj ) );
      
      _mm_store_pd( ocab+4*ma, rcab[0] );
      
    }
    
  }
  
  m = jmax;
  {
    
    ma = ma+1;
    mj = mj+1;
    
    _mm_store_pd( ocab+4*ma+2, _mm_load_pd( icab+2*mj ) );
    
  }
  
}