submodule (lege_poly) bwd_sum
  implicit none; contains
  
  module procedure bwd_sum_c
    integer :: i1
    
    !$omp simd
    do i1 = 1, n
      swork(i1,1) = swork(i1,1) + pmj(i1) * cc(1)
      swork(i1,2) = swork(i1,2) + pmj(i1) * cc(2)
      swork(i1,3) = swork(i1,3) + pmj(i1) * cc(3)
      swork(i1,4) = swork(i1,4) + pmj(i1) * cc(4)
    end do
    
  end procedure bwd_sum_c

end submodule bwd_sum