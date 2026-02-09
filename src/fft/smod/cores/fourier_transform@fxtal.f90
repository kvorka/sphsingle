submodule (fourier_transform) fxtal
  implicit none; contains
  
  module procedure fxztal
    integer :: i, ip, isd, l, k1, icdd
    
    l   = n/2
    i   = it(1)
    isd = 0
    k1  = 1
    ip  = mod(i,4)+2
    
    select case (ip)
      case (4)
        call fxzm4b( m, l, x )
      case (2)
        call fxzm2b( m, l, x )
      case (3)
        call fxzm3b( m, l, x )
      case (5)
        call fxzm5b( m, l, x )
    end select
    
    do icdd = 2, it(2)
      l   = l / ip
      i   = i / 4
      isd = isd + k1 * (ip-1)
      k1  = k1 * ip
      ip  = mod(i,4)+2
      
      select case (ip)
        case (4)
          call fxzm4a( m, k1, l, x, t(2*isd) )
        case (2)
          call fxzm2a( m, k1, l, x, t(2*isd) )
        case (3)
          call fxzm3a( m, k1, l, x, t(2*isd) )
        case (5)
          call fxzm5a( m, k1, l, x, t(2*isd) )
      end select
    end do
    
  end procedure fxztal
  
end submodule fxtal