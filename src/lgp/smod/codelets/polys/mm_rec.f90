submodule (lege_poly) poly_mm
  implicit none; contains
  
  module procedure mm_rec_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: fac
    real(kind=dbl), pointer, contiguous :: scx(:,:)
    
    fac = this%fmj(2,ma)
    
    if ( mod(this%n_dbl,4) == 0 ) then
      
      select case (ma)
        case (1)
          do i2 = 1, this%n_dbl4, 4
            !$omp simd
            do i1 = 1, ndbl
              pmm(i1,i2  ) = fac
              pmm(i1,i2+1) = fac
              pmm(i1,i2+2) = fac
              pmm(i1,i2+3) = fac
            end do
          end do
          
        case default
          scx(1:ndbl,1:this%n_dbl) => this%sinx
          
          do i2 = 1, this%n_dbl4, 4
            !$omp simd aligned (scx:alig)
            do i1 = 1, ndbl
              pmm(i1,i2  ) = fac * scx(i1,i2  ) * pmm(i1,i2  )
              pmm(i1,i2+1) = fac * scx(i1,i2+1) * pmm(i1,i2+1)
              pmm(i1,i2+2) = fac * scx(i1,i2+2) * pmm(i1,i2+2)
              pmm(i1,i2+3) = fac * scx(i1,i2+3) * pmm(i1,i2+3)
            end do
          end do
      end select
      
      scx(1:ndbl,1:this%n_dbl) => this%cosx
      
      do i2 = 1, this%n_dbl4, 4
        !$omp simd aligned (scx:alig)
        do i1 = 1, ndbl
          pmj(i1,i2  ) = pmm(i1,i2  ) / scx(i1,i2  )
          pmj(i1,i2+1) = pmm(i1,i2+1) / scx(i1,i2+1)
          pmj(i1,i2+2) = pmm(i1,i2+2) / scx(i1,i2+2)
          pmj(i1,i2+3) = pmm(i1,i2+3) / scx(i1,i2+3)
        end do
      end do
      
    else
      
      select case (ma)
        case (1)
          do i2 = 1, this%n_dbl4, 4
            !$omp simd
            do i1 = 1, ndbl
              pmm(i1,i2  ) = fac
              pmm(i1,i2+1) = fac
              pmm(i1,i2+2) = fac
              pmm(i1,i2+3) = fac
            end do
          end do
          
          !$omp simd
          do i1 = 1, ndbl
            pmm(i1,this%n_dbl4+1) = fac
            pmm(i1,this%n_dbl4+2) = fac
          end do
          
        case default
          scx(1:ndbl,1:this%n_dbl) => this%sinx
          
          do i2 = 1, this%n_dbl4, 4
            !$omp simd aligned (scx:alig)
            do i1 = 1, ndbl
              pmm(i1,i2  ) = fac * scx(i1,i2  ) * pmm(i1,i2  )
              pmm(i1,i2+1) = fac * scx(i1,i2+1) * pmm(i1,i2+1)
              pmm(i1,i2+2) = fac * scx(i1,i2+2) * pmm(i1,i2+2)
              pmm(i1,i2+3) = fac * scx(i1,i2+3) * pmm(i1,i2+3)
            end do
          end do
          
          !$omp simd
          do i1 = 1, ndbl
            pmm(i1,this%n_dbl4+1) = fac * scx(i1,this%n_dbl4+1) * pmm(i1,this%n_dbl4+1)
            pmm(i1,this%n_dbl4+2) = fac * scx(i1,this%n_dbl4+2) * pmm(i1,this%n_dbl4+2)
          end do
      end select
      
      scx(1:ndbl,1:this%n_dbl) => this%cosx
      
      do i2 = 1, this%n_dbl4, 4
        !$omp simd aligned (scx:alig)
        do i1 = 1, ndbl
          pmj(i1,i2  ) = pmm(i1,i2  ) / scx(i1,i2  )
          pmj(i1,i2+1) = pmm(i1,i2+1) / scx(i1,i2+1)
          pmj(i1,i2+2) = pmm(i1,i2+2) / scx(i1,i2+2)
          pmj(i1,i2+3) = pmm(i1,i2+3) / scx(i1,i2+3)
        end do
      end do
      
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,this%n_dbl4+1) = pmm(i1,this%n_dbl4+1) / scx(i1,this%n_dbl4+1)
        pmj(i1,this%n_dbl4+2) = pmm(i1,this%n_dbl4+2) / scx(i1,this%n_dbl4+2)
      end do
      
    end if
    
    call zero_rarray_sub( this%n, pmj1 )
    
  end procedure mm_rec_sub
  
end submodule poly_mm