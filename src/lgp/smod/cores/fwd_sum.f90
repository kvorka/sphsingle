submodule (lege_poly) fwd_sum_m
  implicit none; contains
  
  module procedure fwd_sum_m_sub
    integer :: ima
    
    !!! After the FFT, we need to shuffle the packing north/south and real/imaginary
    !!! into packing suitable for summation.
    call fwd_shf_sub( n1, wght, cosx, grid, swork )
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    call fwd_set_sub( m, n1, fmj(1), cosx, sinx, swork, pmm, pmj1, pmj, cr )
    
    !! Following with the recursion for degrees m+1 to jmax. We need to repointer our
    !! polynomials, follow with recursion and add cc * pmj to our swork accumulator.
    do ima = 1, nma-1, 2
      call fwd_rec_sub( n1, swork, fmj(2*ima  ), cosx2, pmj,  pmj1, cr(4*ima  ) )
      call fwd_rec_sub( n1, swork, fmj(2*ima+2), cosx2, pmj1, pmj,  cr(4*ima+4) )
    end do
    
    if ( ima == nma ) then
      call fwd_rec_sub( n1, swork, fmj(2*nma), cosx2, pmj, pmj1, cr(4*nma) )
    end if
    
  end procedure fwd_sum_m_sub
  
  module procedure fwd_sum_jmax_sub
    
    !!! After the FFT, we need to shuffle the packing north/south and real/imaginary
    !!! into packing suitable for summation.
    call fwd_shf_sub( n1, wght, cosx, grid, swork )
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    call fwd_set_sub( 1, n1, fmj, cosx, sinx, swork, pmm, pmj1, pmj, cr )
    
  end procedure fwd_sum_jmax_sub
  
end submodule fwd_sum_m