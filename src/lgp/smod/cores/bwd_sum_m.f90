submodule (lege_poly) bwd_sum_m
  implicit none; contains
  
  module procedure bwd_sum_m_sub
    integer                             :: ima
    real(kind=dbl), pointer, contiguous :: pmj2(:)
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    call bwd_set_sub( m, n1, fmj(1), cosx, sinx, cc, pmm, pmj1, pmj, swork )
    
    !! Following with the recursion for degrees m+1 to jmax. We need to repointer our
    !! polynomials, follow with recursion and add cc * pmj to our swork accumulator.
    do ima = 1, nma
      pmj2 => pmj1
      pmj1 => pmj
      pmj  => pmj2
      
      call bwd_rec_sub( n1, cc(4*ima), fmj(2*ima), cosx2, pmj1, pmj, swork )
    end do
    
    !! As we are done with computing the summation, we need to reshufle the data from
    !! packed sum to south/north and real/imaginary parts for upcomming FFT..
    call bwd_shf_sub( n1, cosx, swork, grid )
    
  end procedure bwd_sum_m_sub
  
  module procedure bwd_sum_jmax_sub
    
    !! Starting from degree j equal to order m, we need to forward the recursion for
    !! pmm, prepare the recursion for pmj by setting pmj1 to zero, and we need to
    !! set the initial value of swork to cc * pmj (first member of the sum).
    call bwd_set_sub( 1, n1, fmj, cosx, sinx, cc, pmm, pmj1, pmj, swork )
    
    !! As we are done with computing the summation, we need to reshufle the data from
    !! packed sum to south/north and real/imaginary parts for upcomming FFT..
    call bwd_shf_sub( n1, cosx, swork, grid )
    
  end procedure bwd_sum_jmax_sub
  
end submodule bwd_sum_m