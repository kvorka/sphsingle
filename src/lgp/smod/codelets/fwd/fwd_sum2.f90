submodule (lege_poly) fwd_sum2
  implicit none; contains
  
  module procedure fwd_sum2_sub
    integer        :: i1, i2
    real(kind=dbl) :: cff1, cff2
    
    cff1 = fmj(1)
    cff2 = fmj(2)
    
    !$omp simd
    do i1 = 1, ndbl
      acc(i1,1) = zero
      acc(i1,2) = zero
      acc(i1,3) = zero
      acc(i1,4) = zero
    end do
    
    do i2 = 1, n1
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,i2) = ( cff1 * cosx2(i1,i2) - cff2 ) * pmj1(i1,i2) - pmj(i1,i2)
        
        acc(i1,1) = acc(i1,1) + pmj(i1,i2) * swork(i1,1,i2)
        acc(i1,2) = acc(i1,2) + pmj(i1,i2) * swork(i1,2,i2)
        acc(i1,3) = acc(i1,3) + pmj(i1,i2) * swork(i1,3,i2)
        acc(i1,4) = acc(i1,4) + pmj(i1,i2) * swork(i1,4,i2)
      end do
    end do
    
    !$omp simd
    do i1 = 1, ndbl
      cr(1) = cr(1) + acc(i1,1)
      cr(2) = cr(2) + acc(i1,2)
      cr(3) = cr(3) + acc(i1,3)
      cr(4) = cr(4) + acc(i1,4)
    end do
    
  end procedure fwd_sum2_sub
  
end submodule fwd_sum2