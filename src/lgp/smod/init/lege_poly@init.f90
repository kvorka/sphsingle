submodule (lege_poly) init
  implicit none; contains
  
  module procedure init_lege_sub
    integer                     :: j, m, mj, ma, i, i1, i2
    real(kind=qbl)              :: x1, fx1, x2, fx2, x3, fx3, root, froot
    real(kind=qbl), allocatable :: qamj(:), qemj(:), qfmj(:,:), qroots(:), qpmm(:,:,:)
    
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
    !! instances of this class are alocated. Second, the roots are found by Riddlers method
    !! with initial bracketing proposed by Stjeltjes and Markoff. The weights are additionaly
    !! rescaled by some factor comming from the FFT.
    call alloc_aligned_sub( this%n, this%c_cosx,  this%cosx  )
    call alloc_aligned_sub( this%n, this%c_cosx2, this%cosx2 )
    call alloc_aligned_sub( this%n, this%c_wght,  this%wght  )
    
    allocate( qroots(this%n) )
    
    !$omp parallel do private (x1,fx1,x2,fx2,x3,fx3,root,froot)
    do i2 = 1, this%n
      x1  = cos( (i2-0.5_qbl) * qpi / (2*this%n) )
      fx1 = lege_fn(2*this%n, x1)
      
      x2  = cos( i2 * qpi / (2*this%n+1) )
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
      
      qroots(i2)     = root
      this%cosx(i2)  = real( root, kind=dbl )
      this%cosx2(i2) = real( root**2, kind=dbl )
      this%wght(i2)  = real( qpi * (1-root**2) / ( this%n * lege_fn(2*this%n-1, root) )**2, kind=dbl ) / real( wfac, kind=dbl )
    end do
    !$omp end parallel do
    
    !! Computing coefficients needed for rearranging even/odd degrees before/after transforms. 
    !! The math is done in quadruple precision in order to keep everything as precise as possible, 
    !! while the results are stored in double precision to save space and increase speed.
    allocate( qemj((this%jmax+3)*(this%jmax+2)/2) )
      
    do m = 0, this%jmax+1
      do j = m, this%jmax+1
        qemj(m*(this%jmax+2)-m*(m+1)/2+j+1) = sqrt( ( j**2 - m**2 ) / ( 4*j**2 - 1._qbl ) )
      end do
    end do
    
    !! Computing coefficients needed for rescaling even/odd degrees before/after transforms. 
    !! The math is done in quadruple precision in order to keep everything as precise as possible, 
    !! while the results are stored in double precision to save space and increase speed.
    allocate( qamj(this%nrma) ); ma = 0
    
    do m = 0, this%jmax
      !j = m
        mj = m*(this%jmax+2)-(m-2)*(m+1)/2
        ma = ma+1
        
        qamj(ma) = 1._qbl
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        ma = ma+1
        
        qamj(ma) = 1._qbl / ( qemj(mj) * qemj(mj-1) ) / qamj(ma-1)
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        mj = mj+2
        ma = ma+1
        
        qamj(ma) = 1._qbl / ( qemj(mj) * qemj(mj-1) ) / qamj(ma-1)
      end if
    end do
    
    !! Computing coefficients needed for on-the-fly recursion during transforms. The math is done 
    !! in quadruple precision in order to keep everything as precise as possible, while the results 
    !! are stored in double precision to save space and increase speed.
    allocate( qfmj(2,this%nrma) ) ; ma = 0
    
    do m = 0, this%jmax
      !j = m
        mj = m*(this%jmax+2)-(m-2)*(m+1)/2
        ma = ma+1
        
        if ( m == 0) then
          qfmj(1,ma) = 1._qbl
          qfmj(2,ma) = 1._qbl / sqrt(4*qpi)
        else
          qfmj(1,ma) = 1._qbl
          qfmj(2,ma) = -sqrt( (2*m+1._qbl) / (2*m) )
        end if
      
      do j = 1, (this%jmax-m)/2
        mj = mj+2
        ma = ma+1
        
        qfmj(1,ma) =                                     qamj(ma-1)**2
        qfmj(2,ma) = ( qemj(mj-1)**2 + qemj(mj-2)**2 ) * qamj(ma-1)**2
      end do
      
      if ( mod((this%jmax-m),2) /= 0 ) then
        mj = mj+2
        ma = ma+1
        
        qfmj(1,ma) =                                     qamj(ma-1)**2
        qfmj(2,ma) = ( qemj(mj-1)**2 + qemj(mj-2)**2 ) * qamj(ma-1)**2
      end if
    end do
    
    !! Precomputing the starting points of recursion for each order m. This helps parallelization
    !! and costs only (jmax+1)*n space. Also increases precision, as it is possible again to use
    !! quadruple precision for math, but save the results in doubles.
    allocate( qpmm(ndbl,this%n_dbl,this%jmax+1) )
    
    do m = 0, this%jmax
      ma = this%mamj(m)
      
      select case (ma)
        case (1)
          do i2 = 1, this%n_dbl
            do i1 = 1, ndbl
              qpmm(i1,i2,m+1) = qfmj(2,ma)
            end do
          end do
      
      case default
        do i2 = 1, this%n_dbl
          do i1 = 1, ndbl
            qpmm(i1,i2,m+1) = qfmj(2,ma) * sqrt( 1-qroots(i1+(i2-1)*ndbl)**2 ) * qpmm(i1,i2,m)
          end do
        end do
        
      end select
    end do
    
    do m = 0, this%jmax
      do i2 = 1, this%n_dbl
        do i1 = 1, ndbl
          qpmm(i1,i2,m+1) = qpmm(i1,i2,m+1) / qroots(i1+(i2-1)*ndbl)
        end do
      end do
    end do
    
    !! Prepare the array holders from this class and save the results into the double
    !! precision to save space and also the operation counts.
    allocate( this%emj((this%jmax+3)*(this%jmax+2)/2), this%amj(this%nrma), this%fmj(2,this%nrma) )
    
    this%emj = real( qemj, kind=dbl )
    this%amj = real( qamj, kind=dbl )
    this%fmj = real( qfmj, kind=dbl )
    
    this%c_pmm = malloc( alig, this%n * (this%jmax+1) * size_d )
    call c_f_pointer( this%c_pmm, this%pmm, [ndbl,this%n_dbl,this%jmax+1] )
    
    this%pmm = real( qpmm, kind=dbl )
    
    !! Cleaning.
    deallocate( qpmm, qemj, qamj, qfmj, qroots )
    
  end procedure init_lege_sub
  
end submodule init