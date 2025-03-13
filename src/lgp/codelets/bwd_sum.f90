submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_c
    integer :: i1, i2
    
    do i1 = 1, (n/32)*32, 32
      !$omp simd
      do i2 = 0, 31
        swork(i2+i1,1) = swork(i2+i1,1) + pmj(i2+i1) * cc(1)
        swork(i2+i1,2) = swork(i2+i1,2) + pmj(i2+i1) * cc(2)
        swork(i2+i1,3) = swork(i2+i1,3) + pmj(i2+i1) * cc(3)
        swork(i2+i1,4) = swork(i2+i1,4) + pmj(i2+i1) * cc(4)
      end do
    end do
    
    do i1 = (n/32)*32+1, n, 8
      !$omp simd
      do i2 = 0, 7
        swork(i2+i1,1) = swork(i2+i1,1) + pmj(i2+i1) * cc(1)
        swork(i2+i1,2) = swork(i2+i1,2) + pmj(i2+i1) * cc(2)
        swork(i2+i1,3) = swork(i2+i1,3) + pmj(i2+i1) * cc(3)
        swork(i2+i1,4) = swork(i2+i1,4) + pmj(i2+i1) * cc(4)
      end do
    end do
    
  end procedure bwd_sum_c

end submodule bwd_sum