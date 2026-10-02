module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                             :: nFreq, jmax, jms, n, n_dbl_2, nrma
    integer,        allocatable         :: mamj(:)
    real(kind=dbl), allocatable         :: emj(:), fmj(:), amj(:)
    real(kind=dbl), pointer, contiguous :: cosx(:), sinx(:), cosx2(:), wght(:)
    type(c_ptr)                         :: c_cosx, c_sinx, c_cosx2, c_wght
    
    contains
    
    procedure, public, pass :: init_sub => init_lege_sub
    procedure, public, pass :: bwd_legesum_sub
    procedure, public, pass :: fwd_legesum_sub
    procedure, public, pass :: deallocate_sub => deallocate_lege_sub
    
  end type T_legep
  
  real(kind=qbl), parameter :: qeps  = 1.0d-28
  real(kind=qbl), parameter :: qpi   = acos(-1._qbl)
  real(kind=qbl), parameter :: qzero = 0._qbl
  
  !! Class subroutines
  interface
    module subroutine init_lege_sub(this, jmax, wfac)
      class(T_legep), intent(inout) :: this
      integer,        intent(in)    :: jmax
      integer,        intent(in)    :: wfac
    end subroutine init_lege_sub
    
    module subroutine bwd_legesum_sub(this, cjm, grid)
      class(T_legep),         intent(in)  :: this
      complex(kind=dbl),      intent(in)  :: cjm(*)
      real(kind=dbl), target, intent(out) :: grid(0:*)
    end subroutine bwd_legesum_sub
    
    module subroutine fwd_legesum_sub(this, grid, cjm)
      class(T_legep),         intent(in)    :: this
      real(kind=dbl), target, intent(inout) :: grid(0:*)
      complex(kind=dbl),      intent(out)   :: cjm(*)
    end subroutine fwd_legesum_sub
    
    module subroutine deallocate_lege_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine deallocate_lege_sub
  end interface
  
  !! Cores
  interface
    module subroutine bwd_rxd_sub(jmax, cjm, emj, amj, ocab) &
    & bind(C, name="bwd_rxd_c")
      integer, value,    intent(in)  :: jmax
      real(kind=dbl),    intent(in)  :: emj(*), amj(*)
      complex(kind=dbl), intent(in)  :: cjm(*)
      real(kind=dbl),    intent(out) :: ocab(*)
    end subroutine bwd_rxd_sub
    
    module subroutine bwd_sum_m_sub(n1, m, nma, fmj, cosx, sinx, cosx2, pmm, pmj1, pmj, cc, swork, grid) &
    & bind(C, name="bwd_sum_m_c")
      integer, value, intent(in)    :: n1, m, nma
      real(kind=dbl), intent(in)    :: fmj(*), cosx(*), sinx(*), cosx2(*), cc(*)
      real(kind=dbl), intent(inout) :: pmm(*)
      real(kind=dbl), intent(out)   :: swork(*), grid(*), pmj(*), pmj1(*)
    end subroutine bwd_sum_m_sub
    
    module subroutine bwd_sum_jmax_sub(n1, fmj, cosx, sinx, pmm, pmj1, pmj, cc, swork, grid) &
    & bind(C, name="bwd_sum_jmax_c")
      integer, value, intent(in)    :: n1
      real(kind=dbl), intent(in)    :: fmj(*), cosx(*), sinx(*), cc(*)
      real(kind=dbl), intent(inout) :: pmm(*), pmj(*), pmj1(*)
      real(kind=dbl), intent(out)   :: swork(*), grid(*)
    end subroutine bwd_sum_jmax_sub
    
    module subroutine fwd_sum_m_sub(n1, m, nma, fmj, cosx, sinx, cosx2, wght, pmm, pmj1, pmj, swork, cr, grid) &
    & bind(C, name="fwd_sum_m_c")
      integer, value, intent(in)    :: n1, m, nma
      real(kind=dbl), intent(in)    :: fmj(*), cosx(*), sinx(*), cosx2(*), wght(*), grid(*)
      real(kind=dbl), intent(inout) :: cr(*), pmm(*)
      real(kind=dbl), intent(out)   :: swork(*), pmj(*), pmj1(*)
    end subroutine fwd_sum_m_sub
    
    module subroutine fwd_sum_jmax_sub(n1, fmj, cosx, sinx, wght, pmm, pmj1, pmj, swork, cr, grid) &
    & bind(C, name="fwd_sum_jmax_c")
      integer, value, intent(in)    :: n1
      real(kind=dbl), intent(in)    :: fmj(*), cosx(*), sinx(*), wght(*), grid(*)
      real(kind=dbl), intent(inout) :: cr(*), pmm(*), pmj(*), pmj1(*)
      real(kind=dbl), intent(out)   :: swork(*)
    end subroutine fwd_sum_jmax_sub
    
    module subroutine fwd_rxd_sub(jmax, emj, amj, icab, cjm) &
    & bind(C, name="fwd_rxd_c")
      integer, value,    intent(in)  :: jmax
      real(kind=dbl),    intent(in)  :: emj(*), amj(*), icab(*)
      complex(kind=dbl), intent(out) :: cjm(*)
    end subroutine fwd_rxd_sub
  end interface
  
end module lege_poly