submodule (lege_poly) init
  implicit none; contains
  
  module procedure init_lege_sub
    integer                              :: j, m, mj, ma, i, i1, i2, n, ncnt
    real(kind=qbl)                       :: xincr, x1, fx1, x2, fx2, x3, fx3, root, froot
    real(kind=qbl),          allocatable :: qamj(:), qemj(:)
    real(kind=qbl), target,  allocatable :: xclose(:)
    real(kind=qbl), pointer, contiguous  :: p2xclose(:,:)
    
    !! Set the constants needed within this class. This includes maximum degree,
    !! maximum compound degree, number of roots needed for G.-L. quadrature (keep
    !! in mind, that only half is needed due to equatorial symmetry/antisymmetry
    !! of associated Legendre polynomials) and number of non-zero outputs after
    !! the summation.
    this%jmax = jmax
    this%jms  = jmax*(jmax+1)/2+jmax+1
    
    this%n     = (3*jmax/2+1)/2+1+ndbl-mod((3*jmax/2+1)/2+1,ndbl)
    this%n_dbl = this%n / ndbl
    
    this%nFreq = 4 * this%n * ( this%jmax+1 )
    
    !! This is a wild counter. I did not want to spent time with maths, so I just
    !! copied the cycles developed during writing this library and added counter.
    !! Roughly, it corresponds to number of quadruplets (odd j, even j, real, imag)
    !! encountered during evaluation of summation over degrees.
    this%nrma = 0
      do m = 0, this%jmax
        this%nrma = this%nrma+1
        
        if ( m < this%jmax) then
          do j = 1, (this%jmax-1-m)/2
            this%nrma = this%nrma+1
          end do
          
          this%nrma = this%nrma+1
        end if
      end do
    
    allocate( this%mamj(0:this%jmax+1) )
    
    ma = 0
        
    do m = 0, this%jmax
      !j = m
        ma = ma+1
        this%mamj(m) = ma
      
      do j = 1, (this%jmax-m)/2
        ma = ma+1
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        ma = ma+1
      end if
    end do
    
    this%mamj(this%jmax+1) = ma+1
    
    !! This seeks for roots of the Legendre polynomials. Then computes everything neeeded,
    !! notably the root (cosx), root squared (cosx2) and associated weight (wght). First,
    !! instances of this class are alocated. Second, the positions of the roots are bracketed.
    !! In the end, the roots are found using the Riddlers method.
    this%c_cosx = malloc( alig, this%n * size_d )
    call c_f_pointer( this%c_cosx, this%cosx, [ndbl,this%n_dbl] )
    
    this%c_cosx2 = malloc( alig, this%n * size_d )
    call c_f_pointer( this%c_cosx2, this%cosx2, [ndbl,this%n_dbl] )
    
    this%c_wght = malloc( alig, this%n * size_d )
    call c_f_pointer( this%c_wght, this%wght, [ndbl,this%n_dbl] )
    
    allocate( xclose(this%n) )
    p2xclose(1:ndbl,1:this%n_dbl) => xclose
    
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
    
    p2xclose => null()
    deallocate( xclose )
    
    !! Computing coefficients needed for recurion and associated work on the coefficients,
    !! like rescaling prior and after the recursion. The math is done in quadruple precision
    !! in order to keep everything as precise as possible, while the results are stored in
    !! double precision to save space and increase speed.
    allocate( this%emj((this%jmax+3)*(this%jmax+2)/2), qemj((this%jmax+3)*(this%jmax+2)/2) )
      
    do m = 0, this%jmax+1
      do j = m, this%jmax+1
        qemj(m*(this%jmax+2)-m*(m+1)/2+j+1)     = sqrt((j**2-m**2)/(4*j**2-1._qbl))
        this%emj(m*(this%jmax+2)-m*(m+1)/2+j+1) = real( sqrt((j**2-m**2)/(4*j**2-1._qbl)), kind=dbl )
      end do
    end do
    
    allocate( this%amj(this%nrma), qamj(this%nrma) ); ma = 0
    
    do m = 0, this%jmax
      !j = m
        mj = m*(this%jmax+2)-(m-2)*(m+1)/2
        ma = ma+1
        
        qamj(ma)     = 1._qbl
        this%amj(ma) = real( qamj(ma), kind=dbl )
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        ma = ma+1
        
        qamj(ma)     = 1._qbl / ( qemj(mj) * qemj(mj-1) ) / qamj(ma-1)
        this%amj(ma) = real( qamj(ma), kind=dbl )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        mj = mj+2
        ma = ma+1
        
        qamj(ma)     = 1._qbl / ( qemj(mj) * qemj(mj-1) ) / qamj(ma-1)
        this%amj(ma) = real( qamj(ma), kind=dbl )
      end if
    end do
    
    allocate( this%fmj(2,this%nrma) ) ; ma = 0
    
    do m = 0, this%jmax
      !j = m
        mj = m*(this%jmax+2)-(m-2)*(m+1)/2
        ma = ma+1
        
        if ( m == 0) then
          this%fmj(1,ma) = 1._dbl
          this%fmj(2,ma) = real( 1._qbl / sqrt(4*qpi), kind=dbl )
        else
          this%fmj(1,ma) = 1._dbl
          this%fmj(2,ma) = real( -sqrt( (2*m+1._qbl) / (2*m) ), kind=dbl )
        end if
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        ma = ma+1
        
        this%fmj(1,ma) = real( qamj(ma-1)**2, kind=dbl )
        this%fmj(2,ma) = real( ( qemj(mj-1)**2 + qemj(mj-2)**2 ) * qamj(ma-1)**2, kind=dbl )
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        mj = mj+2
        ma = ma+1
        
        this%fmj(1,ma) = real( qamj(ma-1)**2, kind=dbl )
        this%fmj(2,ma) = real( ( qemj(mj-1)**2 + qemj(mj-2)**2 ) * qamj(ma-1)**2, kind=dbl )
      end if
    end do
    
    deallocate( qemj, qamj )
    
    !! Precomputing the starting points of recursion for each order m. This helps parallelization
    !! and costs only (jmax+1)*n space. Also increases precision, as it is possible again to use
    !! quadruple precision for math, but save the results in doubles.
    this%c_pmm = malloc( alig, this%n * (this%jmax+1) * size_d )
    call c_f_pointer( this%c_pmm, this%pmm, [ndbl,this%n_dbl,this%jmax+1] )
    
    do m = 0, this%jmax
      ma = this%mamj(m)
      
      select case (ma)
        case (1)
          do i2 = 1, this%n_dbl
            !$omp simd
            do i1 = 1, ndbl
              this%pmm(i1,i2,m+1) = this%fmj(2,ma)
            end do
          end do
      
      case default
        do i2 = 1, this%n_dbl
          !$omp simd
          do i1 = 1, ndbl
            this%pmm(i1,i2,m+1) = this%fmj(2,ma) * sqrt( 1-this%cosx(i1,i2)**2 ) * this%pmm(i1,i2,m)
          end do
        end do
        
      end select
    end do
    
    do m = 0, this%jmax
      do i2 = 1, this%n_dbl
        !$omp simd
        do i1 = 1, ndbl
          this%pmm(i1,i2,m+1) = this%pmm(i1,i2,m+1) / this%cosx(i1,i2)
        end do
      end do
    end do
    
    !! The last bit of operation is rescaling the weights with whatever factor does the FFT
    !! include or does not include.
    this%wght = this%wght / real(wfac, kind=dbl)
    
  end procedure init_lege_sub
  
end submodule init