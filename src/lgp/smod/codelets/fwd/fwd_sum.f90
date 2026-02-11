submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: p1(:,:), s1(:,:), s2(:,:), s3(:,:), s4(:,:)
    
    do i2 = 1, this%n32, 32
      p1(1:16,1:2) => pmj(i2:i2+31)
      s1(1:16,1:2) => swork(i2:i2+31,1)
      s2(1:16,1:2) => swork(i2:i2+31,2)
      s3(1:16,1:2) => swork(i2:i2+31,3)
      s4(1:16,1:2) => swork(i2:i2+31,4)
      
      !$omp simd aligned (p1,s1,s2,s3,s4:alig)
      do i1 = 1, 16
        cr(1) = cr(1) + p1(i1,1) * s1(i1,1) + p1(i1,2) * s1(i1,2)
        cr(2) = cr(2) + p1(i1,1) * s2(i1,1) + p1(i1,2) * s2(i1,2)
        cr(3) = cr(3) + p1(i1,1) * s3(i1,1) + p1(i1,2) * s3(i1,2)
        cr(4) = cr(4) + p1(i1,1) * s4(i1,1) + p1(i1,2) * s4(i1,2)
      end do
    end do
    
    do i2 = this%n32+1, this%n, 16
      p1(1:16,1:1) => pmj(i2:i2+15)
      s1(1:16,1:1) => swork(i2:i2+15,1)
      s2(1:16,1:1) => swork(i2:i2+15,2)
      s3(1:16,1:1) => swork(i2:i2+15,3)
      s4(1:16,1:1) => swork(i2:i2+15,4)
      
      !$omp simd aligned (p1,s1,s2,s3,s4:alig)
      do i1 = 1, 16
        cr(1) = cr(1) + p1(i1,1) * s1(i1,1)
        cr(2) = cr(2) + p1(i1,1) * s2(i1,1)
        cr(3) = cr(3) + p1(i1,1) * s3(i1,1)
        cr(4) = cr(4) + p1(i1,1) * s4(i1,1)
      end do
    end do
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
