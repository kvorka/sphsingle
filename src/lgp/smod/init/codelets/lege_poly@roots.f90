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
    integer                             :: i, i1, i2, n, ncnt
    real(kind=qbl)                      :: xincr, x1, fx1, x2, fx2, x3, fx3, root, froot
    real(kind=qbl), target, allocatable :: xclose(:)
    real(kind=qbl), pointer, contiguous :: p2xclose(:,:)
    
    !!**********************************************************************!!
    !!* Close to roots array holder and holder arrays.                     *!!
    !!**********************************************************************!!
    this%c_cosx = malloc( alig, this%n * size_d )
    call c_f_pointer( this%c_cosx, this%cosx, [ndbl,this%n_dbl] )
    
    this%c_cosx2 = malloc( alig, this%n * size_d )
    call c_f_pointer( this%c_cosx2, this%cosx2, [ndbl,this%n_dbl] )
    
    this%c_wght = malloc( alig, this%n * size_d )
    call c_f_pointer( this%c_wght, this%wght, [ndbl,this%n_dbl] )
    
    !!**********************************************************************!!
    !!* Bracket the positions of the roots.                                *!!
    !!**********************************************************************!!
    allocate( xclose(this%n) ); p2xclose(1:ndbl,1:this%n_dbl) => xclose
    
    n    = this%n**2 / 4
    ncnt = 0
    
    do while ( ncnt < this%n )
      n     = 6 * n / 5
      xincr = 1._qbl / n
      ncnt  = 0
      
      !$omp parallel do private (fx1, fx2)
      do i = 1, n
        fx1 = lege_fn( 2*this%n, (i-1) * xincr )
        fx2 = lege_fn( 2*this%n, (i  ) * xincr )
        
        if ( fx1 * fx2 < qzero ) then
          !$omp critical
          ncnt         = ncnt+1
          xclose(ncnt) = (i-1) * xincr
          !$omp end critical
        end if
      end do
      !$omp end parallel do
    end do
    
    !!**********************************************************************!!
    !!* Riddler                                                            *!!
    !!**********************************************************************!!
    !$omp parallel do private (x1,fx1,x2,fx2,x3,fx3,root,froot)
    do i2 = 1, this%n_dbl
      do i1 = 1, ndbl
        x1  = p2xclose(i1,i2)
        fx1 = lege_fn(2*this%n, x1)
        
        x2  = x1+xincr
        fx2 = lege_fn(2*this%n, x2)
        
        do
          x3  = ( x1 + x2 ) / 2
          fx3 = lege_fn(2*this%n, x3)
          
          root  = x3 + (x3-x1) * sign(1._qbl,fx1-fx2) * fx3 / sqrt( fx3**2 - fx1*fx2 )
          froot = lege_fn(2*this%n, root)
          
          if ( abs(froot) < qeps ) then
            exit
          else if ( fx3 * froot < qzero ) then
            x1  = x3
            fx1 = fx3
            x2  = root
            fx2 = froot
          else if ( fx1 * froot < qzero ) then
            x1  = root
            fx1 = froot
          else if ( fx2 * froot < qzero ) then
            x2  = root
            fx2 = froot
          end if
        end do
        
        this%cosx(i1,i2)  = real( root, kind=dbl )
        this%cosx2(i1,i2) = real( root**2, kind=dbl )
        this%wght(i1,i2)  = real( qpi * (1-root**2) / ( this%n * lege_fn(2*this%n-1, root) )**2, kind=dbl )
      end do
    end do
    !$omp end parallel do
    
    !!**********************************************************************!!
    !!* Cleaning.                                                          *!!
    !!**********************************************************************!!
    p2xclose => null(); deallocate( xclose )
    
  end procedure find_roots_sub
  
end submodule roots
