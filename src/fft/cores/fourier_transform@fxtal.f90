submodule (fourier_transform) fxtal
  implicit none; contains
  
  module procedure fxztal
    integer :: it, ip, isd, l, k1, icdd
    
    l   = this%n/2
    it  = this%it(l-1)
    isd = 0
    k1  = 1
    ip  = mod(it,4)+2
    
    select case (ip)
      case (4)
        call fxzm4b_c( howmany, l/4, x )
      case (2)
        call fxzm2b_c( howmany, l/2, x )
      case (3)
        call fxzm3b_c( howmany, l/3, x )
      case (5)
        call fxzm5b_c( howmany, l/5, x )
    end select
    
    do icdd = 2, this%it(this%n/2)
      l   = l / ip
      it  = it / 4
      isd = isd + k1 * (ip-1)
      k1  = k1 * ip
      ip  = mod(it,4)+2
      
      select case (ip)
        case (4)
          call fxzm4a_c( howmany, k1, l/4, x, this%t(1+2*isd) )
        case (2)
          call fxzm2a_c( howmany, k1, l/2, x, this%t(1+2*isd) )
        case (3)
          call fxzm3a_c( howmany, k1, l/3, x, this%t(1+2*isd) )
        case (5)
          call fxzm5a_c( howmany, k1, l/5, x, this%t(1+2*isd) )
      end select
    end do
    
  end procedure fxztal
  
end submodule fxtal
