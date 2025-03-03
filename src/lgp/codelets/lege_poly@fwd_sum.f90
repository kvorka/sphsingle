submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer :: i1
    
    do i1 = 1, step
      cr(1) = cr(1) + pmj(i1) * swork(i1,1)
      cr(2) = cr(2) + pmj(i1) * swork(i1,2)
      cr(3) = cr(3) + pmj(i1) * swork(i1,3)
      cr(4) = cr(4) + pmj(i1) * swork(i1,4)
    end do
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
