submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mm_set_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: pmm(:,:)
    
    pmm(1:16,1:this%n_16) => this%pmm(1:this%n,m+1)
    
    do i2 = 1, this%n_16
      !$omp simd aligned (pmm:alig)
      do i1 = 1, 16
        pmj1(i1,i2) = 0._dbl
        pmj(i1,i2)  = pmm(i1,i2)
      end do
    end do
    
  end procedure mm_set_sub
  
end submodule poly_mm