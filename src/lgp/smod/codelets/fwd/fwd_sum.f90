submodule (lege_poly) fwd_sum
  implicit none; contains
  
  module procedure fwd_sum_sub
    integer                             :: i1, i2
    real(kind=dbl)                      :: c1, c2, c3, c4
    real(kind=dbl), pointer, contiguous :: p1(:), p2(:), p3(:), p4(:), s11(:), s21(:), s31(:), s41(:), &
                                                                     & s12(:), s22(:), s32(:), s42(:), &
                                                                     & s13(:), s23(:), s33(:), s43(:), &
                                                                     & s14(:), s24(:), s34(:), s44(:)
    
    c1 = cr(1)
    c2 = cr(2)
    c3 = cr(3)
    c4 = cr(4)
    
    do i2 = 1, this%n64, 64
      p1(1:16)  => pmj(i2   :i2+15)
      p2(1:16)  => pmj(i2+16:i2+31)
      p3(1:16)  => pmj(i2+32:i2+47)
      p4(1:16)  => pmj(i2+48:i2+63)
      
      s11(1:16) => swork(i2:i2+15,1)
      s21(1:16) => swork(i2:i2+15,2)
      s31(1:16) => swork(i2:i2+15,3)
      s41(1:16) => swork(i2:i2+15,4)
      
      s12(1:16) => swork(i2+16:i2+31,1)
      s22(1:16) => swork(i2+16:i2+31,2)
      s32(1:16) => swork(i2+16:i2+31,3)
      s42(1:16) => swork(i2+16:i2+31,4)
      
      s13(1:16) => swork(i2+32:i2+47,1)
      s23(1:16) => swork(i2+32:i2+47,2)
      s33(1:16) => swork(i2+32:i2+47,3)
      s43(1:16) => swork(i2+32:i2+47,4)
      
      s14(1:16) => swork(i2+48:i2+63,1)
      s24(1:16) => swork(i2+48:i2+63,2)
      s34(1:16) => swork(i2+48:i2+63,3)
      s44(1:16) => swork(i2+48:i2+63,4)
      
      !$omp simd aligned (p1,p2,p3,p4,s11,s12,s13,s14,s21,s22,s23,s24,s31,s32,s33,s34,s41,s42,s43,s44:alig)
      do i1 = 1, 16
        c1 = c1 + p1(i1) * s11(i1) + p2(i1) * s12(i1) + p3(i1) * s13(i1) + p4(i1) * s14(i1)
        c2 = c2 + p1(i1) * s21(i1) + p2(i1) * s22(i1) + p3(i1) * s23(i1) + p4(i1) * s24(i1)
        c3 = c3 + p1(i1) * s31(i1) + p2(i1) * s32(i1) + p3(i1) * s33(i1) + p4(i1) * s34(i1)
        c4 = c4 + p1(i1) * s41(i1) + p2(i1) * s42(i1) + p3(i1) * s43(i1) + p4(i1) * s44(i1)
      end do
    end do
    
    do i2 = this%n64+1, this%n, 16
      p1(1:16)  => pmj(i2:i2+15)
      
      s11(1:16) => swork(i2:i2+15,1)
      s21(1:16) => swork(i2:i2+15,2)
      s31(1:16) => swork(i2:i2+15,3)
      s41(1:16) => swork(i2:i2+15,4)
      
      !$omp simd aligned (p1,s12,s21,s31,s41:alig)
      do i1 = 1, 16
        c1 = c1 + p1(i1) * s11(i1)
        c2 = c2 + p1(i1) * s21(i1)
        c3 = c3 + p1(i1) * s31(i1)
        c4 = c4 + p1(i1) * s41(i1)
      end do
    end do
    
    cr(1) = c1
    cr(2) = c2
    cr(3) = c3
    cr(4) = c4
    
  end procedure fwd_sum_sub
  
end submodule fwd_sum
