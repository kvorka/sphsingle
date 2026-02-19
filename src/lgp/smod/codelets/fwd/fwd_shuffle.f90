submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer        :: i1, i2
    real(kind=dbl) :: w, cx, g1, g2, g3, g4
    
    do i2 = 1, n1
      !$omp simd
      do i1 = 1, ndbl
        w  = wght(i1,i2)
        cx = cosx(i1,i2)
        
        g1 = grid(i1,i2,1)
        g2 = grid(i1,i2,2)
        g3 = grid(i1,i2,3)
        g4 = grid(i1,i2,4)
        
        swork(i1,1,i2) = ( g1 - g2 ) * w
        swork(i1,3,i2) = ( g1 + g2 ) * w * cx
        swork(i1,2,i2) = ( g3 - g4 ) * w
        swork(i1,4,i2) = ( g3 + g4 ) * w * cx
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle