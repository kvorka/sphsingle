module lege_poly
  use math
  implicit none
  
  type, public :: T_legep
    integer                             :: nFreq, jmax, jms, n, n_dbl_2, nrma
    integer,        allocatable         :: mamj(:)
    real(kind=dbl), allocatable         :: emj(:), fmj(:,:), amj(:)
    real(kind=dbl), pointer, contiguous :: cosx(:), sinx(:), cosx2(:), wght(:)
    type(c_ptr)                         :: c_cosx, c_sinx, c_cosx2, c_wght
    
    contains
    
    procedure, public,  pass :: init_sub       => init_lege_sub
    procedure, public,  pass :: deallocate_sub => deallocate_lege_sub
    
    procedure, public, pass :: index_bwd_sub, bwd_legesum_sub
    procedure, public, pass :: index_fwd_sub, fwd_legesum_sub
    
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
    
    module subroutine deallocate_lege_sub(this)
      class(T_legep), intent(inout) :: this
    end subroutine deallocate_lege_sub
    
    module subroutine index_bwd_sub(this, cjm, rcab)
      class(T_legep),    intent(in)  :: this
      complex(kind=dbl), intent(in)  :: cjm(*)
      real(kind=dbl),    intent(out) :: rcab(*)
    end subroutine index_bwd_sub
    
    module subroutine bwd_legesum_sub(this, cc, grid)
      class(T_legep),         intent(in)  :: this
      real(kind=dbl),         intent(in)  :: cc(4,*)
      real(kind=dbl), target, intent(out) :: grid(0:*)
    end subroutine bwd_legesum_sub
    
    module subroutine index_fwd_sub(this, rcab, cjm)
      class(T_legep),    intent(in)    :: this
      real(kind=dbl),    intent(inout) :: rcab(*)
      complex(kind=dbl), intent(out)   :: cjm(*)
    end subroutine index_fwd_sub
    
    module subroutine fwd_legesum_sub(this, grid, cr)
      class(T_legep),         intent(in)    :: this
      real(kind=dbl), target, intent(inout) :: grid(0:*)
      real(kind=dbl),         intent(out)   :: cr(4,*)
    end subroutine fwd_legesum_sub
  end interface
  
  !! Cores
  interface
    module subroutine bwd_sum_m_sub(n1, ma1, ma2, fmj, cosx, sinx, cosx2, pmm, pmj1, pmj, cc, swork, grid)
      integer,                             intent(in)    :: n1, ma1, ma2
      real(kind=dbl),                      intent(in)    :: fmj(2,ma1:*), cosx(*), sinx(*), cosx2(*), cc(4,ma1:*)
      real(kind=dbl),                      intent(inout) :: pmm(*)
      real(kind=dbl),                      intent(out)   :: swork(*), grid(*)
      real(kind=dbl), pointer, contiguous, intent(inout) :: pmj(:), pmj1(:)
    end subroutine bwd_sum_m_sub
    
    module subroutine bwd_sum_jmax_sub( n1, fmj, cosx, sinx, pmm, pmj1, pmj, cc, swork, grid )
      integer,        intent(in)    :: n1
      real(kind=dbl), intent(in)    :: fmj(*), cosx(*), sinx(*), cc(*)
      real(kind=dbl), intent(inout) :: pmm(*), pmj(*), pmj1(*)
      real(kind=dbl), intent(out)   :: swork(*), grid(*)
    end subroutine bwd_sum_jmax_sub
    
    module subroutine fwd_sum_m_sub(n1, ma1, ma2, fmj, cosx, sinx, cosx2, wght, pmm, pmj1, pmj, swork, cr, grid)
      integer,        intent(in)                         :: n1, ma1, ma2
      real(kind=dbl), intent(in)                         :: fmj(2,ma1:*), cosx(*), sinx(*), cosx2(*), wght(*), grid(*)
      real(kind=dbl), intent(inout)                      :: cr(4,ma1:*), pmm(*)
      real(kind=dbl), intent(out)                        :: swork(*)
      real(kind=dbl), pointer, contiguous, intent(inout) :: pmj(:), pmj1(:)
    end subroutine fwd_sum_m_sub
    
    module subroutine fwd_sum_jmax_sub(n1, fmj, cosx, sinx, wght, pmm, pmj1, pmj, swork, cr, grid)
      integer,        intent(in)    :: n1
      real(kind=dbl), intent(in)    :: fmj(*), cosx(*), sinx(*), wght(*), grid(*)
      real(kind=dbl), intent(inout) :: cr(*), pmm(*), pmj(*), pmj1(*)
      real(kind=dbl), intent(out)   :: swork(*)
    end subroutine fwd_sum_jmax_sub
  end interface
  
  !! Codelets
  interface
    module subroutine bwd_c2r_sub(jmax, cjm, cab) bind(C, name="bwd_c2r_c")
      integer, value,    intent(in)  :: jmax
      complex(kind=dbl), intent(in)  :: cjm(*)
      real(kind=dbl),    intent(out) :: cab(*)
    end subroutine bwd_c2r_sub
    
    module subroutine bwd_rxd_sub(jmax, emj, amj, icab, ocab) bind(C, name="bwd_rxd_c")
      integer, value, intent(in)  :: jmax
      real(kind=dbl), intent(in)  :: emj(*), amj(*), icab(*)
      real(kind=dbl), intent(out) :: ocab(*)
    end subroutine bwd_rxd_sub
    
    module subroutine fwd_rxd_sub(jmax, emj, amj, icab, ocab) bind(C, name="fwd_rxd_c")
      integer, value, intent(in)  :: jmax
      real(kind=dbl), intent(in)  :: emj(*), amj(*), icab(*)
      real(kind=dbl), intent(out) :: ocab(*)
    end subroutine fwd_rxd_sub
    
    module subroutine fwd_r2c_sub(jmax, cab, cjm) bind(C, name="fwd_r2c_c")
      integer, value,    intent(in)  :: jmax
      real(kind=dbl),    intent(in)  :: cab(*)
      complex(kind=dbl), intent(out) :: cjm(*)
    end subroutine fwd_r2c_sub
    
    module subroutine bwd_set_sub(ma, n, fmj, cosx, sinx, cc, pmm, pmj1, pmj, swork) bind(C, name="bwd_set_c")
      integer, value, intent(in)    :: n, ma
      real(kind=dbl), intent(in)    :: fmj(*), cosx(*), sinx(*), cc(*)
      real(kind=dbl), intent(inout) :: pmm(*)
      real(kind=dbl), intent(out)   :: pmj1(*), pmj(*), swork(*)
    end subroutine bwd_set_sub
    
    module subroutine bwd_rec_sub(n, cc, fmj, cosx2, pmj1, pmj, swork) bind(C, name="bwd_rec_c")
      integer, value, intent(in)    :: n
      real(kind=dbl), intent(in)    :: cc(*), fmj(*), cosx2(*), pmj1(*)
      real(kind=dbl), intent(inout) :: pmj(*), swork(*)
    end subroutine bwd_rec_sub
    
    module subroutine bwd_shf_sub(n, cosx, swork, grid) bind(C, name="bwd_shf_c")
      integer, value, intent(in)  :: n
      real(kind=dbl), intent(in)  :: cosx(*), swork(*)
      real(kind=dbl), intent(out) :: grid(*)
    end subroutine bwd_shf_sub
    
    module subroutine fwd_shf_sub(n, wght, cosx, grid, swork) bind(C, name="fwd_shf_c")
      integer, value, intent(in)  :: n
      real(kind=dbl), intent(in)  :: wght(*), cosx(*), grid(*)
      real(kind=dbl), intent(out) :: swork(*)
    end subroutine fwd_shf_sub
    
    module subroutine fwd_set_sub(ma, n, fmj, cosx, sinx, swork, pmm, pmj1, pmj, cr) bind(C, name="fwd_set_c")
      integer, value, intent(in)    :: n, ma
      real(kind=dbl), intent(in)    :: swork(*), fmj(*), cosx(*), sinx(*)
      real(kind=dbl), intent(out)   :: pmj1(*), pmj(*)
      real(kind=dbl), intent(inout) :: cr(*), pmm(*)
    end subroutine fwd_set_sub
    
    module subroutine fwd_rec_sub(n, swork, fmj, cosx2, pmj1, pmj, cr) bind(C, name="fwd_rec_c")
      integer, value, intent(in)    :: n
      real(kind=dbl), intent(in)    :: swork(*), fmj(*), cosx2(*), pmj1(*)
      real(kind=dbl), intent(inout) :: pmj(*), cr(*)
    end subroutine fwd_rec_sub
  end interface
  
end module lege_poly
