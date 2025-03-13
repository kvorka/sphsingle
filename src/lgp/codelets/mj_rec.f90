submodule (lege_poly) poly_mj
  implicit none; contains
  
  module procedure mj_rec_c
    integer :: i1, i2
    
    do i1 = 1, (n/32)*32, 32
      !$omp simd
      do i2 = 0, 31
        pmj(i2+i1)  = ( cff(1) * cosx2(i2+i1) - cff(2) ) * pmj1(i2+i1) - pmj(i2+i1)
      end do
    end do
    
    do i1 = (n/32)*32+1, n, 8
      !$omp simd
      do i2 = 0, 7
        pmj(i2+i1)  = ( cff(1) * cosx2(i2+i1) - cff(2) ) * pmj1(i2+i1) - pmj(i2+i1)
      end do
    end do
    
  end procedure mj_rec_c
  
end submodule poly_mj