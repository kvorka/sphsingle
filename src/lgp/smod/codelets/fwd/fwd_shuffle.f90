submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer        :: i1, i2, n2
    real(kind=dbl) :: w1, w2, cx1, cx2, s11, s21, s31, s41, s12, s22, s32, s42
    
    n2 = ( n1 / 2 ) * 2
    
    do i2 = 1, n2, 2
      !$omp simd
      do i1 = 1, ndbl
        w1 = wght(i1,i2  )
        w2 = wght(i1,i2+1)
        
        cx1 = cosx(i1,i2  )
        cx2 = cosx(i1,i2+1)
        
        s11 = g1(i1,i2)
        s21 = g2(i1,i2)
        s31 = g3(i1,i2)
        s41 = g4(i1,i2)
        
        s12 = g1(i1,i2+1)
        s22 = g2(i1,i2+1)
        s32 = g3(i1,i2+1)
        s42 = g4(i1,i2+1)
        
        g1(i1,i2) = ( s11 - s21 ) * w1
        g3(i1,i2) = ( s11 + s21 ) * w1 * cx1
        g2(i1,i2) = ( s31 - s41 ) * w1
        g4(i1,i2) = ( s31 + s41 ) * w1 * cx1
        
        g1(i1,i2+1) = ( s12 - s22 ) * w2
        g3(i1,i2+1) = ( s12 + s22 ) * w2 * cx2
        g2(i1,i2+1) = ( s32 - s42 ) * w2
        g4(i1,i2+1) = ( s32 + s42 ) * w2 * cx2
      end do
    end do
    
    if ( n2 /= n1 ) then
      !$omp simd
      do i1 = 1, ndbl
        w1 = wght(i1,n1)
        cx1 = cosx(i1,n1)
        
        s11 = g1(i1,n1)
        s21 = g2(i1,n1)
        s31 = g3(i1,n1)
        s41 = g4(i1,n1)
        
        g1(i1,n1) = ( s11 - s21 ) * w1
        g3(i1,n1) = ( s11 + s21 ) * w1 * cx1
        g2(i1,n1) = ( s31 - s41 ) * w1
        g4(i1,n1) = ( s31 + s41 ) * w1 * cx1
      end do
    end if
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle