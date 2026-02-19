submodule (lege_poly) bwd_shuffle
  implicit none; contains
  
  module procedure bwd_shuffle_sub
    integer        :: i1, i2
    real(kind=dbl) :: cx, s1, s2, s3, s4
    
    do i2 = 1, n1
      !$omp simd
      do i1 = 1, ndbl
        cx = cosx(i1,i2)
        
        s1 = swork(i1,1,i2)
        s2 = swork(i1,3,i2)
        s3 = swork(i1,2,i2)
        s4 = swork(i1,4,i2)
        
        grid(i1,i2,1) = s2 * cx + s1
        grid(i1,i2,2) = s2 * cx - s1
        grid(i1,i2,3) = s4 * cx + s3
        grid(i1,i2,4) = s4 * cx - s3
      end do
    end do
    
  end procedure bwd_shuffle_sub

end submodule bwd_shuffle