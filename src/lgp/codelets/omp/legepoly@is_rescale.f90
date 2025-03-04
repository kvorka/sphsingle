submodule (lege_poly) is_rescale
  implicit none; contains
  
  module procedure is_rescale_c
    integer :: i1, i2
    
    do i2 = 1, n
      !$omp simd
      do i1 = 1, 4
        rcab(i1,i2) = cff(i2) * rcab(i1,i2)
      end do
    end do
    
  end procedure is_rescale_c
  
end submodule is_rescale