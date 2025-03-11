#include <stdlib.h>
#include <immintrin.h>

extern inline void fwd_indx_c( const int jmax,                // cutoff degree
                               const double *restrict emj,    // factors
                               const double *restrict ocab,   // input to rescale
                               double *restrict icab) {       // output
  
  __m128d rcab[4];
  
  int m, ma, mj;
  
  m = 0;
  {
    
    { // j == m
      
      ma = 0;
      mj = 0;
      
      rcab[3] = _mm_load_pd( ocab );
      
      _mm_store_pd( icab, _mm_load_pd( ocab+2 ) );
      
    }
    
    for ( int j = 0; j < (jmax-1)/2; j++ ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj-1] ), rcab[3] );
      
      rcab[3] = _mm_load_pd( ocab+4*ma   );
      rcab[2] = _mm_load_pd( ocab+4*ma+2 );
      
      rcab[1] = _mm_mul_pd( _mm_set1_pd( emj[mj  ] ), rcab[3] );
      
      rcab[1] = _mm_add_pd( rcab[0], rcab[1] );
      
      _mm_store_pd( icab+2*mj-2, rcab[1] );
      _mm_store_pd( icab+2*mj  , rcab[2] );
      
    }
    
    if ( jmax%2 == 0 ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[1] = _mm_load_pd( ocab+4*ma   );
      rcab[2] = _mm_load_pd( ocab+4*ma+2 );
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj-1] ), rcab[3] );
      rcab[1] = _mm_mul_pd( _mm_set1_pd( emj[mj  ] ), rcab[1] );
      
      rcab[1] = _mm_add_pd( rcab[0], rcab[1] );
      
      _mm_store_pd( icab+2*mj-2, rcab[1] );
      _mm_store_pd( icab+2*mj  , rcab[2] );
      
    } else {
      
      ma = ma+1;
      mj = mj+1;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj  ] ), rcab[3] );
      rcab[1] = _mm_mul_pd( _mm_set1_pd( emj[mj+1] ), _mm_load_pd( ocab+4*ma ) );
      
      rcab[1] = _mm_add_pd( rcab[0], rcab[1] );
      
      _mm_store_pd( icab+2*mj, rcab[1] );
      
    }
    
  }
  
  for ( m = 1; m < jmax; m++ ) {
    
    { // j == m
      
      ma = ma+1;
      mj = mj+1;
      
      rcab[3] = _mm_load_pd( ocab+4*ma );
      
      _mm_store_pd( icab+2*mj, _mm_load_pd( ocab+4*ma+2 ) );
      
    }
    
    for ( int j = 0; j < (jmax-m-1)/2; j++ ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj+m-1] ), rcab[3] );
      
      rcab[3] = _mm_load_pd( ocab+4*ma   );
      rcab[2] = _mm_load_pd( ocab+4*ma+2 );
      
      rcab[1] = _mm_mul_pd( _mm_set1_pd( emj[mj+m  ] ), rcab[3] );
      
      rcab[1] = _mm_add_pd( rcab[0], rcab[1] );
      
      _mm_store_pd( icab+2*mj-2, rcab[1] );
      _mm_store_pd( icab+2*mj  , rcab[2] );
      
    }
    
    if ( (jmax-m)%2 == 0 ) {
      
      ma = ma+1;
      mj = mj+2;
      
      rcab[1] = _mm_load_pd( ocab+4*ma   );
      rcab[2] = _mm_load_pd( ocab+4*ma+2 );
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj+m-1] ), rcab[3] );
      rcab[1] = _mm_mul_pd( _mm_set1_pd( emj[mj+m  ] ), rcab[1] );
      
      rcab[1] = _mm_add_pd( rcab[0], rcab[1] );
      
      _mm_store_pd( icab+2*mj-2, rcab[1] );
      _mm_store_pd( icab+2*mj  , rcab[2] );
      
    } else {
      
      ma = ma+1;
      mj = mj+1;
      
      rcab[1] = _mm_load_pd( ocab+4*ma   );
      
      rcab[0] = _mm_mul_pd( _mm_set1_pd( emj[mj+m  ] ), rcab[3] );
      rcab[1] = _mm_mul_pd( _mm_set1_pd( emj[mj+m+1] ), rcab[1] );
      
      rcab[1] = _mm_add_pd( rcab[0], rcab[1] );
      
      _mm_store_pd( icab+2*mj, rcab[1] );
      
    }
    
  }
  
  m = jmax;
  {
    
    ma = ma+1;
    mj = mj+1;
    
    _mm_store_pd( icab+2*mj, _mm_load_pd( ocab+4*ma+2 ) );
    
  }
  
}