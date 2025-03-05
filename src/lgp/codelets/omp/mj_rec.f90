submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mj_rec_c
    integer :: i2
    
    !$omp simd
    do i2 = 1, n
      pmj(i2)  = ( cff(1) * cosx2(i2) - cff(2) ) * pmj1(i2) - pmj(i2)
    end do
    
  end procedure mj_rec_c
  
end submodule poly_mj