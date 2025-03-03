submodule (lateral_grid) harmsy
  implicit none; contains
  
  module procedure harmsy_sub
    integer                     :: itheta, i2, i1
    real(kind=dbl), allocatable :: cosx(:), sinx(:), cosx2(:), pmm(:), pmj1(:), pmj(:)
    real(kind=dbl), allocatable :: sumN(:), sumS(:), swork(:), rcc(:)
    
    !Transform to suitable real input
    call this%lgp%allocate_lgp_arr_sub( rcc )
    call this%lgp%index_bwd_sub( cin, rcc )
    
    !Allocating work memory
    allocate( swork(4*step), sumN(step*this%fourtrans%n), sumS(step*this%fourtrans%n), &
            & pmm(step), pmj1(step), pmj(step), cosx(step), sinx(step), cosx2(step)    )
    
    !Cycle over latitudes :: calculating step at once
    do itheta = 1, (this%lgp%n/step)*step, step
      do i1 = 1, step
        cosx(i1)  = this%lgp%rw(itheta+i1-1,1)
        sinx(i1)  = this%lgp%rw(itheta+i1-1,2)
        cosx2(i1) = this%lgp%rw(itheta+i1-1,3)
      end do
      
      call zero_rarray_sub( step*this%fourtrans%n, sumN )
      call zero_rarray_sub( step*this%fourtrans%n, sumS )
      
      call this%lgp%bwd_legesum_sub( rcc, sumN, sumS, cosx, sinx, cosx2, pmm, pmj1, pmj, swork )
      
      call this%fourtrans%fft_c2r_sub( sumN )
      call this%fourtrans%fft_c2r_sub( sumS )
      
      do i2 = 0, this%fourtrans%n-1
        do i1 = 1, step
          grid(itheta+i1-1,i2+1,1) = sumN(i1+i2*step)
          grid(itheta+i1-1,i2+1,2) = sumS(i1+i2*step)
        end do
      end do
    end do
    
    !Cleaning
    deallocate( swork, sumN, sumS, pmm, pmj1, pmj, cosx, sinx, cosx2, rcc )
    
  end procedure harmsy_sub
  
end submodule harmsy
