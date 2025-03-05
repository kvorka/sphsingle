submodule (lege_poly) bwd
  implicit none; contains
  
  module procedure bwd_legesum_sub
    integer                             :: m, j, ma
    real(kind=dbl), pointer, contiguous :: work(:), swork(:), pmm(:), pmj1(:), pmj(:), pmj2(:)
    type(c_ptr)                         :: c_work
    
    call alloc_aligned_sub( default_alig, 7*this%n, c_work, work )
    
    pmm   => work(          1 :   this%n )
    pmj1  => work(   this%n+1 : 2*this%n )
    pmj   => work( 2*this%n+1 : 3*this%n )
    swork => work( 3*this%n+1 : 7*this%n )
    
    ma  = 0
    
    do m = 0, this%jmax
      call zero_rarray_c( 4*this%n, swork )
      
      !j = m
        ma = ma+1
        
        call mm_set_c( ma, this%n, this%fmj(2,ma), this%cosx, this%sinx, pmm, pmj1, pmj )
        call bwd_sum_c( this%n, pmj, cc(1,ma), swork )
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
        
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call mj_rec_c( this%n, this%fmj(1,ma), this%cosx2, pmj1, pmj )
        call bwd_sum_c( this%n, pmj, cc(1,ma), swork )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        ma = ma+1
        
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call mj_rec_c( this%n, this%fmj(1,ma), this%cosx2, pmj1, pmj )
        call bwd_sum_c( this%n, pmj, cc(1,ma), swork )
      end if
      
      call bwd_shuffle_c( this%n, this%cosx, swork, grid(1,m) )
    end do
    
    call free( c_work )
    
  end procedure bwd_legesum_sub
  
end submodule bwd