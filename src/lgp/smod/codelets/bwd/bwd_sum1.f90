submodule (lege_poly) bwd_sum1
  implicit none; contains
  
  module procedure bwd_sum1_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: fac, c1, c2, c3, c4
    real(kind=dbl), pointer, contiguous :: scx(:,:,:)
    
    fac = this%fmj(2,ma)
    
    select case (ma)
      case (1)
        do i2 = 1, this%n_dbl_2
          !$omp simd
          do i1 = 1, ndbl
            pmm(i1,1,i2) = fac
            pmm(i1,2,i2) = fac
            
            pmj1(i1,1,i2) = zero
            pmj1(i1,2,i2) = zero
          end do
        end do
        
      case default
        scx(1:ndbl,1:2,1:this%n_dbl_2) => this%sinx
        
        do i2 = 1, this%n_dbl_2
          !$omp simd aligned (scx:alig)
          do i1 = 1, ndbl
            pmm(i1,1,i2) = fac * scx(i1,1,i2) * pmm(i1,1,i2)
            pmm(i1,2,i2) = fac * scx(i1,2,i2) * pmm(i1,2,i2)
            
            pmj1(i1,1,i2) = zero
            pmj1(i1,2,i2) = zero
          end do
        end do
    end select
    
    scx(1:ndbl,1:2,1:this%n_dbl_2) => this%cosx
    
    c1 = cc(1)
    c2 = cc(2)
    c3 = cc(3)
    c4 = cc(4)
    
    do i2 = 1, this%n_dbl_2
      !$omp simd
      do i1 = 1, ndbl
        pmj(i1,1,i2) = pmm(i1,1,i2) / scx(i1,1,i2)
        pmj(i1,2,i2) = pmm(i1,2,i2) / scx(i1,2,i2)
        
        swork(i1,1,i2) = pmj(i1,1,i2) * c1
        swork(i1,5,i2) = pmj(i1,2,i2) * c1
        swork(i1,2,i2) = pmj(i1,1,i2) * c2
        swork(i1,6,i2) = pmj(i1,2,i2) * c2
        swork(i1,3,i2) = pmj(i1,1,i2) * c3
        swork(i1,7,i2) = pmj(i1,2,i2) * c3
        swork(i1,4,i2) = pmj(i1,1,i2) * c4
        swork(i1,8,i2) = pmj(i1,2,i2) * c4
      end do
    end do
    
  end procedure bwd_sum1_sub
  
end submodule bwd_sum1