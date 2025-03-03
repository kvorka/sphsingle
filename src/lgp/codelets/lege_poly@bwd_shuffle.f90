submodule (lege_poly) bwd_shuffle
  implicit none; contains
  
  module procedure bwd_shuffle_sub
    integer :: i1
    
    do i1 = 1, step
      sumN(i1,1) = swork(i1,1,2) * cosx(i1) + swork(i1,1,1)
      sumN(i1,2) = swork(i1,2,2) * cosx(i1) + swork(i1,2,1)
      sumS(i1,1) = swork(i1,1,2) * cosx(i1) - swork(i1,1,1)
      sumS(i1,2) = swork(i1,2,2) * cosx(i1) - swork(i1,2,1)
    end do
    
  end procedure bwd_shuffle_sub

end submodule bwd_shuffle