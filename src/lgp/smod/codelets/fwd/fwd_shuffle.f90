submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: w(:), cosx(:), g1(:), g2(:), g3(:), g4(:), s1(:), s2(:), s3(:), s4(:)
    
    do i2 = 1, this%n, 16
      cosx(1:16) => this%cosx(i2:i2+15)
      w(1:16)    => this%wght(i2:i2+15)
      
      g1(1:16) => grid(i2:i2+15,1,1)
      g2(1:16) => grid(i2:i2+15,2,1)
      g3(1:16) => grid(i2:i2+15,1,2)
      g4(1:16) => grid(i2:i2+15,2,2)
      
      s1(1:16) => swork(i2:i2+15,1,1)
      s2(1:16) => swork(i2:i2+15,2,1)
      s3(1:16) => swork(i2:i2+15,1,2)
      s4(1:16) => swork(i2:i2+15,2,2)
      
      !$omp simd aligned (w,cosx,g1,g2,g3,g4,s1,s2,s3,s4:alig)
      do i1 = 1, 16
        s1(i1) = ( g1(i1) - g2(i1) ) * w(i1)
        s2(i1) = ( g3(i1) - g4(i1) ) * w(i1)
        s3(i1) = ( g1(i1) + g2(i1) ) * w(i1) * cosx(i1)
        s4(i1) = ( g3(i1) + g4(i1) ) * w(i1) * cosx(i1)
      end do
    end do
    
  end procedure fwd_shuffle_sub
  
end submodule fwd_shuffle