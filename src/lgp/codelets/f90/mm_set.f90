submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mm_set_c
    integer :: i2
    
    select case (ma)
      case (1)
        do i2 = 1, n
          pmm(i2) = cff
        end do
      
      case default
        do i2 = 1, n
          pmm(i2) = cff * sinx(i2) * pmm(i2)
        end do
    end select
    
    call zero_rarray_c( n, pmj1 )
    
    do i2 = 1, n
      pmj(i2)  = pmm(i2) / cosx(i2)
    end do
    
  end procedure mm_set_c
  
end submodule poly_mm