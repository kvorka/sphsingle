submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mm_set_sub
    integer                             :: i1, i2, n64
    real(kind=dbl), pointer, contiguous :: p_pmm(:), p_pmj(:), p_pmj1(:)
    
    n64 = (this%n/64)*64
    
    do i2 = 1, n64, 64
      p_pmm(1:64)   => this%pmm(i2:i2+63,m+1)
      p_pmj(1:64)   => pmj(i2:i2+63)
      p_pmj1(1:64)  => pmj1(i2:i2+63)
      
      !$omp simd aligned (p_pmm,p_pmj,p_pmj1:alig)
      do i1 = 1, 64
        p_pmj1(i1) = 0._dbl
        p_pmj(i1)  = p_pmm(i1)
      end do
    end do
    
    do i2 = n64+1, this%n, 16
      p_pmm(1:16)   => this%pmm(i2:i2+15,m+1)
      p_pmj(1:16)   => pmj(i2:i2+15)
      p_pmj1(1:16)  => pmj1(i2:i2+15)
      
      !$omp simd aligned (p_pmm,p_pmj,p_pmj1:alig)
      do i1 = 1, 16
        p_pmj1(i1) = 0._dbl
        p_pmj(i1)  = p_pmm(i1)
      end do
    end do
    
  end procedure mm_set_sub
  
end submodule poly_mm