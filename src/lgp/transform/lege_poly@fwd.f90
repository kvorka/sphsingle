submodule (lege_poly) fwd
  implicit none; contains

  module procedure fwd_legesum_sub
    integer                             :: j, m, ma
    real(kind=dbl), pointer, contiguous :: work(:), swork(:), pmm(:), pmj1(:), pmj(:), pmj2(:)
    type(c_ptr)                         :: c_work
    
    call alloc_aligned_sub( default_alig, 7*this%n, c_work, work )
    
    pmm   => work(          1 :   this%n )
    pmj1  => work(   this%n+1 : 2*this%n )
    pmj   => work( 2*this%n+1 : 3*this%n )
    swork => work( 3*this%n+1 : 7*this%n )
    
    ma  = 0
    
    do m = 0, this%jmax
      call fwd_shuffle_c( this%n, this%cosx, this%wght, sumN(1,m), sumS(1,m), swork )
      
      !j = m
        ma = ma+1
        
        call mm_set_c( ma, this%n, this%fmj(2,ma), this%cosx, this%sinx, pmm, pmj1, pmj )
        call fwd_sum_c( this%n, pmj, swork, cr(1,ma) )
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
        
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call mj_rec_c( this%n, this%fmj(1,ma), this%cosx2, pmj1, pmj )
        call fwd_sum_c( this%n, pmj, swork, cr(1,ma) )
      end do
      
      if ( mod(this%jmax-m,2) /= 0 ) then
        ma = ma+1
        
        pmj2 => pmj1
        pmj1 => pmj
        pmj  => pmj2
        
        call mj_rec_c( this%n, this%fmj(1,ma), this%cosx2, pmj1, pmj )
        call fwd_sum_c( this%n, pmj, swork, cr(1,ma) )
      end if
    end do
    
    call free( c_work )
    
  end procedure fwd_legesum_sub
  
end submodule fwd