#include <stdlib.h>
#include <immintrin.h>

#if defined( avx )
  
  typedef __m256d mmreg;
  
  #ifndef incr
  #define incr 4
  #endif
  
  inline __m256d load( const double *arr) { return _mm256_load_pd( arr ); }
  inline void store( double *arr, const __m256d reg ) { _mm256_store_pd( arr, reg ); }
  
  inline __m256d broadcast( const double *num ) { return _mm256_broadcast_sd( num ); }
  inline __m256d setzero() { return _mm256_setzero_pd(); }
  inline __m256d setdbl( const double cff ) { return _mm256_set1_pd( cff ); }
  
  inline __m256d add( const __m256d reg1, const __m256d reg2 ) { return _mm256_add_pd( reg1, reg2 ); }
  inline __m256d sub( const __m256d reg1, const __m256d reg2 ) { return _mm256_sub_pd( reg1, reg2 ); }
  inline __m256d mul( const __m256d reg1, const __m256d reg2 ) { return _mm256_mul_pd( reg1, reg2 ); }
  inline __m256d dvv( const __m256d reg1, const __m256d reg2 ) { return _mm256_div_pd( reg1, reg2 ); }
  
  inline double hadd( const __m256d reg ) {
    
    __m128d rsum = _mm_add_pd( _mm256_extractf128_pd( reg, 0 ), _mm256_extractf128_pd( reg, 1 ) );
    return _mm_cvtsd_f64( _mm_add_pd( rsum, _mm_unpackhi_pd( rsum, rsum ) ) );
    
  }
  
#elif defined( avx512 )
  
  typedef __m512d mmreg;
  
  #ifndef incr
  #define incr 8
  #endif
  
  inline __m512d load( const double *arr) { return _mm512_load_pd( arr ); }
  inline void store( double *arr, const __m512d reg ) { _mm512_store_pd( arr, reg ); }
  
  inline __m512d broadcast( const double *num ) { return _mm512_broadcast_sd( num ); }
  
  inline __m512d add( const __m512d reg1, const __m512d reg2) { return _mm512_add_pd( reg1, reg2 ); }
  inline __m512d sub( const __m512d reg1, const __m512d reg2) { return _mm512_sub_pd( reg1, reg2 ); }
  inline __m512d mul( const __m512d reg1, const __m512d reg2) { return _mm512_mul_pd( reg1, reg2 ); }
  
#endif