submodule (lateral_grid) harman
  implicit none; contains
  
  module procedure harman_sub
    integer                     :: itheta, i1, i2
    real(kind=dbl), allocatable :: cosx(:), sinx(:), cosx2(:), wght(:), pmm(:), pmj1(:), pmj(:)
    real(kind=dbl), allocatable :: sumN(:), sumS(:), swork(:), rcr(:)
    
    !Allocate input array
    call this%lgp%allocate_lgp_arr_sub( rcr )
    
    !Allocating work memory
    allocate( swork(4*step), sumN(step*this%fourtrans%n), sumS(step*this%fourtrans%n), pmm(step), &
            & pmj1(step), pmj(step), cosx(step), sinx(step), cosx2(step), wght(step)              )
    
    !Cycle over latitudes :: computing step at once
    do itheta = 1, (this%lgp%n/step)*step, step
      do i1 = 1, step
        cosx(i1)  = this%lgp%rw(itheta+i1-1,1)
        sinx(i1)  = this%lgp%rw(itheta+i1-1,2)
        cosx2(i1) = this%lgp%rw(itheta+i1-1,3)
        wght(i1)  = this%lgp%rw(itheta+i1-1,4)
      end do
      
      do i2 = 0, this%fourtrans%n-1
        do i1 = 1, step
          sumN(i1+i2*step) = grid(itheta+i1-1,i2+1,1)
          sumS(i1+i2*step) = grid(itheta+i1-1,i2+1,2)
        end do
      end do
      
      call this%fourtrans%fft_r2c_sub( sumN )
      call this%fourtrans%fft_r2c_sub( sumS )
      
      call this%lgp%fwd_legesum_sub( sumN, sumS, rcr, cosx, sinx, cosx2, wght, pmm, pmj1, pmj, swork )
    end do
    
    !Reindex output array
    call this%lgp%index_fwd_sub( rcr, cout )
    
    !Cleaning
    deallocate( swork, sumN, sumS, pmm, pmj1, pmj, cosx, sinx, cosx2, wght, rcr )
    
  end procedure harman_sub
  
end submodule harman
