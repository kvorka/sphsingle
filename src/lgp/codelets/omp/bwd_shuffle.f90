submodule (lege_poly) bwd_shuffle
  implicit none; contains
  
  module procedure bwd_shuffle_c
    integer :: i1
    
    !$omp simd
    do i1 = 1, n
      grid(i1,1,1) = swork(i1,1,2) * cosx(i1) + swork(i1,1,1)
      grid(i1,2,1) = swork(i1,1,2) * cosx(i1) - swork(i1,1,1)
      grid(i1,1,2) = swork(i1,2,2) * cosx(i1) + swork(i1,2,1)
      grid(i1,2,2) = swork(i1,2,2) * cosx(i1) - swork(i1,2,1)
    end do
    
  end procedure bwd_shuffle_c

end submodule bwd_shuffle