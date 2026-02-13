submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mm_set_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: pmm(:,:)
    
    pmm(1:ndbl,1:this%n_dbl) => this%pmm(1:this%n,m+1)
    
    do i2 = 1, this%n_dbl, 2
      !$omp simd aligned (pmm:alig)
      do i1 = 1, ndbl
        pmj1(i1,i2  ) = zero
        pmj1(i1,i2+1) = zero

        pmj(i1,i2  ) = pmm(i1,i2  )
        pmj(i1,i2+1) = pmm(i1,i2+1)
      end do
    end do
    
  end procedure mm_set_sub
  
end submodule poly_mm