submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer :: i1, i2
    
    !$omp simd
    do i1 = 1, ndbl
      acc(i1,1) = 0._dbl
      acc(i1,2) = 0._dbl
      acc(i1,3) = 0._dbl
      acc(i1,4) = 0._dbl
    end do
    
    do i2 = 1, this%n_dbl, 2
      !$omp simd
      do i1 = 1, ndbl
        acc(i1,1) = acc(i1,1) + pmj(i1,i2) * swork(i1,i2,1) + pmj(i1,i2+1) * swork(i1,i2+1,1)
        acc(i1,2) = acc(i1,2) + pmj(i1,i2) * swork(i1,i2,2) + pmj(i1,i2+1) * swork(i1,i2+1,2)
        acc(i1,3) = acc(i1,3) + pmj(i1,i2) * swork(i1,i2,3) + pmj(i1,i2+1) * swork(i1,i2+1,3)
        acc(i1,4) = acc(i1,4) + pmj(i1,i2) * swork(i1,i2,4) + pmj(i1,i2+1) * swork(i1,i2+1,4)
      end do
    end do
    
    !$omp simd
    do i1 = 1, ndbl
      cr(1) = cr(1) + acc(i1,1)
      cr(2) = cr(2) + acc(i1,2)
      cr(3) = cr(3) + acc(i1,3)
      cr(4) = cr(4) + acc(i1,4)
    end do
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
