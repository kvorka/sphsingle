submodule (math) lege
  implicit none; contains
  
  module procedure lege_fn
    real(kind=qbl) :: p1, p2, fac
    integer        :: i
    
    p1      = 1._qbl
    lege_fn = x
    
    do i = 2, deg
      fac     = 2 - 1._qbl / i
      
      p2      = fac * ( lege_fn * x - p1 ) + p1
      p1      = lege_fn
      lege_fn = p2
    end do
    
  end procedure lege_fn
  
end submodule lege