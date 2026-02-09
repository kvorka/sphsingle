submodule (lege_poly) polyinit
  implicit none; contains
  
  module procedure compute_pmm_sub
    integer :: m, ma, i2
    
    call alloc_aligned_2d_sub( this%n, this%jmax+1, this%c_pmm, this%pmm )
    
    do m = 0, this%jmax
      ma = this%mamj(m)
      
      select case (ma)
        case (1)
          do i2 = 1, this%n
            this%pmm(i2,m+1) = this%fmj(2,ma)
          end do
      
      case default
        do i2 = 1, this%n
          this%pmm(i2,m+1) = this%fmj(2,ma) * sqrt( 1-this%cosx(i2)**2 ) * this%pmm(i2,m)
        end do
        
      end select
    end do
    
    do m = 0, this%jmax
      do i2 = 1, this%n
        this%pmm(i2,m+1) = this%pmm(i2,m+1) / this%cosx(i2)
      end do
    end do
    
  end procedure compute_pmm_sub
  
end submodule polyinit