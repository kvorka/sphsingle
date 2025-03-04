submodule (lege_poly) roots
  implicit none; contains
  
   real(kind=qbl) function lege_fn(deg, x)
    integer,        intent(in) :: deg
    real(kind=qbl), intent(in) :: x
    real(kind=qbl)             :: p1, p2, fac
    integer                    :: i
    
    p1      = 1._qbl
    lege_fn = x
    
    do i = 2, deg
      fac     = 2 - 1._qbl / i
      
      p2      = fac * ( lege_fn * x - p1 ) + p1
      p1      = lege_fn
      lege_fn = p2
    end do
    
  end function lege_fn
  
  module procedure find_roots_sub
    integer                     :: i, n, ncnt
    real(kind=qbl)              :: xincr, x1, fx1, x2, fx2, root, froot
    real(kind=qbl), allocatable :: xclose(:)
    
    !!**********************************************************************!!
    !!* Close to roots array holder and holder arrays.                     *!!
    !!**********************************************************************!!
    call alloc_aligned_sub( default_alig, this%n, this%c_cosx,  this%cosx  )
    call alloc_aligned_sub( default_alig, this%n, this%c_sinx,  this%sinx  )
    call alloc_aligned_sub( default_alig, this%n, this%c_cosx2, this%cosx2 )
    call alloc_aligned_sub( default_alig, this%n, this%c_wght,  this%wght  )
    
    !!**********************************************************************!!
    !!* Seek for efficient stepping to use within the bisection method and *!!
    !!* starting points [xclose,xclose+xincr].                             *!!
    !!**********************************************************************!!
    allocate( xclose(this%n) )
    
    n = this%n**2 / 4
    
    do
      n     = 6 * n / 5
      xincr = 1._qbl / n
      ncnt  = 0
      
      !$omp parallel do private (fx1, fx2)
      do i = 1, n
        fx1 = lege_fn( 2*this%n, (i-1) * xincr )
        fx2 = lege_fn( 2*this%n, (i  ) * xincr )
        
        if ( fx1 * fx2 < 0._qbl ) then
          !$omp critical
          ncnt         = ncnt+1
          xclose(ncnt) = (i-1) * xincr
          !$omp end critical
        end if
      end do
      !$omp end parallel do
      
      if ( ncnt == this%n ) then
        exit
      else
        do i = 1, ncnt
          xclose(i) = 0._qbl
        end do
      end if
    end do
    
    !!**********************************************************************!!
    !!* Bisection                                                          *!!
    !!**********************************************************************!!
    !$omp parallel do private (x1,fx1,x2,fx2,root,froot)
    do i = 1, this%n
      x1  = xclose(i)
      fx1 = lege_fn(2*this%n, x1)
      
      x2  = x1+xincr
      fx2 = lege_fn(2*this%n, x2)
      
      do
        root  = ( x1 + x2 ) / 2
        froot = lege_fn(2*this%n, root)
        
        if ( abs(froot) < qeps ) then
          exit
        else if ( fx1 * froot < 0._qbl ) then
          x2  = root
          fx2 = froot
        else
          x1  = root
          fx1 = froot
        end if
      end do
      
      this%cosx(i)  = root
      this%sinx(i)  = sqrt( 1 - root**2 )
      this%cosx2(i) = root**2
      this%wght(i) = qpi * (1-root**2) / ( this%n * lege_fn(2*this%n-1, root) )**2
    end do
    !$omp end parallel do
    
    !!**********************************************************************!!
    !!* Cleaning.                                                          *!!
    !!**********************************************************************!!
    deallocate( xclose )
    
  end procedure find_roots_sub
  
end submodule roots