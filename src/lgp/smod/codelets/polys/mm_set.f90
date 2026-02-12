submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mm_set_sub
    integer                             :: i1, i2
    real(kind=dbl), pointer, contiguous :: pmm(:,:,:)
    
    pmm(1:ndbl,1:2,1:this%n_step) => this%pmm(1:this%n,m+1)
    
    do i2 = 1, this%nUnroll3, 3
      !$omp simd aligned (pmm:alig)
      do i1 = 1, ndbl
        pmj1(i1,1,i2  ) = 0._dbl
        pmj1(i1,2,i2  ) = 0._dbl
        pmj1(i1,1,i2+1) = 0._dbl
        pmj1(i1,2,i2+1) = 0._dbl
        pmj1(i1,1,i2+2) = 0._dbl
        pmj1(i1,2,i2+2) = 0._dbl
        
        pmj(i1,1,i2  ) = pmm(i1,1,i2  )
        pmj(i1,2,i2  ) = pmm(i1,2,i2  )
        pmj(i1,1,i2+1) = pmm(i1,1,i2+1)
        pmj(i1,2,i2+1) = pmm(i1,2,i2+1)
        pmj(i1,1,i2+2) = pmm(i1,1,i2+2)
        pmj(i1,2,i2+2) = pmm(i1,2,i2+2)
      end do
    end do
    
    do i2 = this%nUnroll3+1, this%n_step
      !$omp simd aligned (pmm:alig)
      do i1 = 1, ndbl
        pmj1(i1,1,i2) = 0._dbl
        pmj1(i1,2,i2) = 0._dbl
        
        pmj(i1,1,i2) = pmm(i1,1,i2)
        pmj(i1,2,i2) = pmm(i1,2,i2)
      end do
    end do
    
  end procedure mm_set_sub
  
end submodule poly_mm