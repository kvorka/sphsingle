submodule (lege_poly) fwd_shuffle
  implicit none; contains
  
  module procedure fwd_shuffle_c
    integer :: i1
    
    !$omp simd
    do i1 = 1, n
      swork(i1,1,1) = ( sumN(i1,1) - sumS(i1,1) ) * w(i1)
      swork(i1,2,1) = ( sumN(i1,2) - sumS(i1,2) ) * w(i1)
      swork(i1,1,2) = ( sumN(i1,1) + sumS(i1,1) ) * w(i1) * cosx(i1)
      swork(i1,2,2) = ( sumN(i1,2) + sumS(i1,2) ) * w(i1) * cosx(i1)
    end do
    
  end procedure fwd_shuffle_c
  
end submodule fwd_shuffle