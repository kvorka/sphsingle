module sphsingle
  use lateral_grid
  use physical_grid
  implicit none; public
  
  !! Only numbers that can be expressed as 2^a.3^b.5^c-3 with a > 0 are supported due to FFT constraints.
  integer, parameter :: addmissible_jmax(48) = [   5,   7,   9,  13,  15,  21,  27,  29,  33,  37,  45,  47,  51,  57,  61,   69, & 
                                               &  77,  87,  93,  97, 105, 117, 125, 141, 147, 157, 159, 177, 189, 197, 213,  237, &
                                               & 247, 253, 267, 285, 297, 317, 321, 357, 381, 397, 429, 447, 477, 497, 997, 1021  ]
  
  type(T_lateralGrid), private   :: sph
  integer,             protected :: nth, nph
  
  contains
  
  subroutine init_sphpack(jcut)
    integer, intent(in) :: jcut
    
    if ( .not. ( any( addmissible_jmax == jcut ) ) ) then
      write(*,*) 'Due to FFT, this value of jmax is prohibited. Check sphsingle.f90 for supported values.'
      stop
    end if
    
    call sph%init_sub( jcut )
    
    nth = 2 * sph%lgp%n
    nph = sph%fft%n
    
  end subroutine init_sphpack
  
  subroutine alloc_grid(grid)
    class(T_grid), intent(out) :: grid
    
    call grid%alloc_sub( nth, nph )
    
  end subroutine alloc_grid
  
  subroutine free_grid(grid)
    class(T_grid), intent(inout) :: grid
    
    call grid%free_sub()
    
  end subroutine free_grid
  
  subroutine harmsy(cajm, grid)
    complex(kind=dbl), intent(in)  :: cajm(*)
    real(kind=dbl),    intent(out) :: grid(*)
    
    call sph%harmsy_sub( cajm, grid )
    
  end subroutine harmsy
  
  subroutine harman(grid, cajm)
    real(kind=dbl),    intent(inout) :: grid(*)
    complex(kind=dbl), intent(out)   :: cajm(*)
    
    call sph%harman_sub( grid, cajm )
    
  end subroutine harman
  
  subroutine clean_sphpack()
    
    call sph%deallocate_sub()
    
  end subroutine clean_sphpack
  
end module sphsingle