#include "../../../math/cvec.h"

// Definition of constant needed for fft
#ifndef rc31
#define rc31 setdbl( -0.5 )
#endif

#ifndef rc32
#define rc32 setdbl( +0.86602540378443864676 )
#endif

#ifndef rc51
#define rc51 setdbl( +0.25 )
#endif

#ifndef rc52
#define rc52 setdbl( +0.5590169943749474241 )
#endif

#ifndef rc53
#define rc53 setdbl( +0.6180339887498948482 )
#endif

#ifndef rc54
#define rc54 setdbl( -0.9510565162951535721 )
#endif

// Definition of special register values: 2 and -1
#if defined( fma ) || defined( avx512fma )
  
  #ifndef rtwo
  #define rtwo setdbl(+2.0)
  #endif
  
  #ifndef rm00
  #define rm00 setdbl( -0.0 )
  #endif
  
#endif