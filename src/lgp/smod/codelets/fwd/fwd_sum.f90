submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer                             :: i1, i2, n64
    real(kind=dbl)                      :: c1, c2, c3, c4
    real(kind=dbl), pointer, contiguous :: p(:), s1(:), s2(:), s3(:), s4(:)
    
    n64 = (this%n/64)*64
    
    c1 = cr(1)
    c2 = cr(2)
    c3 = cr(3)
    c4 = cr(4)
    
    do i2 = 1, n64, 64
      p(1:64)  => pmj(i2:i2+63)
      
      s1(1:64) => swork(i2:i2+63,1)
      s2(1:64) => swork(i2:i2+63,2)
      s3(1:64) => swork(i2:i2+63,3)
      s4(1:64) => swork(i2:i2+63,4)
      
      !$omp simd aligned (p,s1,s2,s3,s4:alig)
      do i1 = 1, 64
        c1 = c1 + p(i1) * s1(i1)
        c2 = c2 + p(i1) * s2(i1)
        c3 = c3 + p(i1) * s3(i1)
        c4 = c4 + p(i1) * s4(i1)
      end do
    end do
    
    do i2 = n64+1, this%n, 16
      p(1:16)  => pmj(i2:i2+15)
      
      s1(1:16) => swork(i2:i2+15,1)
      s2(1:16) => swork(i2:i2+15,2)
      s3(1:16) => swork(i2:i2+15,3)
      s4(1:16) => swork(i2:i2+15,4)
      
      !$omp simd aligned (p,s1,s2,s3,s4:alig)
      do i1 = 1, 16
        c1 = c1 + p(i1) * s1(i1)
        c2 = c2 + p(i1) * s2(i1)
        c3 = c3 + p(i1) * s3(i1)
        c4 = c4 + p(i1) * s4(i1)
      end do
    end do
    
    cr(1) = c1
    cr(2) = c2
    cr(3) = c3
    cr(4) = c4
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
