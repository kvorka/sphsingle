submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mm_set_c
    integer :: i1, i2
    
    call zero_rarray_c( n, pmj1 )
    
    do i2 = 1, (n/32)*32, 32
      !$omp simd
      do i1 = 0, 31
        pmj(i1+i2)  = pmm(i1+i2)
      end do
    end do
    
    do i2 = (n/32)*32+1, n, 8
      !$omp simd
      do i1 = 0, 7
        pmj(i1+i2)  = pmm(i1+i2)
      end do
    end do
    
  end procedure mm_set_c
  
end submodule poly_mm