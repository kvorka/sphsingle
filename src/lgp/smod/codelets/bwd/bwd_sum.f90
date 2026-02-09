submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_sub
    integer :: i1, i2
    
    do i2 = 1, this%n32, 32
      !$omp simd
      do i1 = 0, 31
        swork(i1+i2,1) = swork(i1+i2,1) + pmj(i1+i2) * cc(1)
        swork(i1+i2,2) = swork(i1+i2,2) + pmj(i1+i2) * cc(2)
        swork(i1+i2,3) = swork(i1+i2,3) + pmj(i1+i2) * cc(3)
        swork(i1+i2,4) = swork(i1+i2,4) + pmj(i1+i2) * cc(4)
      end do
    end do
    
    do i2 = this%n32+1, this%n, 16
      !$omp simd
      do i1 = 0, 15
        swork(i1+i2,1) = swork(i1+i2,1) + pmj(i1+i2) * cc(1)
        swork(i1+i2,2) = swork(i1+i2,2) + pmj(i1+i2) * cc(2)
        swork(i1+i2,3) = swork(i1+i2,3) + pmj(i1+i2) * cc(3)
        swork(i1+i2,4) = swork(i1+i2,4) + pmj(i1+i2) * cc(4)
      end do
    end do
    
  end procedure bwd_sum_sub
  
end submodule bwd_sum