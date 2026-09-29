submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_m_sub
    integer                             :: ima
    real(kind=dbl), pointer, contiguous :: pmj2(:)
    
    !!! After the FFT, we need to shuffle the packing north/south and real/imaginary
    !!! into packing suitable for summation.
    call fwd_shf_sub( n1, wght, cosx, grid, swork )
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    call fwd_set_sub( n1, swork, pmm, pmj1, pmj, cr(1,ma1) )
    
    !! Following with the recursion for degrees m+1 to jmax. We need to repointer our
    !! polynomials, follow with recursion and add cc * pmj to our swork accumulator.
    do ima = ma1+1, ma2
      pmj2 => pmj1
      pmj1 => pmj
      pmj  => pmj2
      
      call fwd_rec_sub( n1, swork, fmj(1,ima), cosx2, pmj1, pmj, cr(1,ima) )
    end do
    
  end procedure fwd_sum_m_sub
  
end submodule fwd_sum