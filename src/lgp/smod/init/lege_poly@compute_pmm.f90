submodule (lege_poly) compute_pmm
  implicit none; contains
  
  module procedure compute_pmm_sub
    integer :: m, ma, i1, i2
    
    this%c_pmm = malloc( alig, this%n * (this%jmax+1) * size_d )
    call c_f_pointer( this%c_pmm, this%pmm, [ndbl,this%n_dbl,this%jmax+1] )
    
    do m = 0, this%jmax
      ma = this%mamj(m)
      
      select case (ma)
        case (1)
          do i2 = 1, this%n_dbl
            !$omp simd
            do i1 = 1, ndbl
              this%pmm(i1,i2,m+1) = this%fmj(2,ma)
            end do
          end do
      
      case default
        do i2 = 1, this%n_dbl
          !$omp simd
          do i1 = 1, ndbl
            this%pmm(i1,i2,m+1) = this%fmj(2,ma) * sqrt( 1-this%cosx(i1,i2)**2 ) * this%pmm(i1,i2,m)
          end do
        end do
        
      end select
    end do
    
    do m = 0, this%jmax
      do i2 = 1, this%n_dbl
        !$omp simd
        do i1 = 1, ndbl
          this%pmm(i1,i2,m+1) = this%pmm(i1,i2,m+1) / this%cosx(i1,i2)
        end do
      end do
    end do
    
  end procedure compute_pmm_sub
  
end submodule compute_pmm