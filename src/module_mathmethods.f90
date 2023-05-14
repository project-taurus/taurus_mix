!==============================================================================!
! MODULE MathMethods                                                           !
!                                                                              !
! This module contains the variables and routines related to mathematical fun- !
! ctions and numerical methods.                                                !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine multiply_weights_diag                                           !
! - subroutine multiply_weights_full                                           !
! - subroutine ClebschGordan                                                   !
! - function kdelta                                                            !
!==============================================================================!
MODULE MathMethods 

use Constants

implicit none
public

CONTAINS 

!------------------------------------------------------------------------------!
! function kdelta                                                              !
!                                                                              !
! Computes the Kronecker delta: \delta_ij = 1 if i = j                         !
!                                         = 0 otherwise                        !
!                                                                              !
! Input: i,j = values to comare                                                !
!                                                                              !
! Output: delta = 0 or 1                                                       !
!------------------------------------------------------------------------------!
function kdelta(i,j) result(delta)

integer, intent(in) :: i, j
integer :: delta

if ( i == j ) then 
  delta = 1    
else
  delta = 0    
endif

end function kdelta

!------------------------------------------------------------------------------!
! subroutine multiply_weights_diag                                             !
!                                                                              !
! Computes the diagonal part of: C = F^T A F, which corresponds to the expect- !
! ation values of the operators A for the different eigenstates of H.          !
!                                                                              !
! Input: F =  weights f                                                        !
!        A =  projected matrix elements of operator A                          !
!                                                                              !
! Output: D = PGCM matrix elements of operator A                               !
!------------------------------------------------------------------------------!
subroutine multiply_weights_diag(F,A,D,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: A, F  
real(r64), dimension(ndim), intent(out) :: D
real(r64), dimension(ndim,ndim) :: B, C
integer :: i 

call dgemm('t','n',ndim,ndim,ndim,one,F,ndim,A,ndim,zero,B,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,B,ndim,F,ndim,zero,C,ndim)

do i = 1, ndim
  D(i) = C(i,i)
enddo

end subroutine multiply_weights_diag

!------------------------------------------------------------------------------!
! subroutine multiply_weights_full                                             !
!                                                                              !
! Computes the matrix elements multiplication: C = F^T A G, which is used to   !
! get the transition probabilities between different states.                   !
!                                                                              !
! Input: F,G =  weights f, g (different states)                                !
!        A =  projected matrix elements of operator A                          !
!                                                                              !
! Output: C = PGCM matrix elements of operator A                               !
!------------------------------------------------------------------------------!
subroutine multiply_weights_full(F,A,G,C,ndim)

integer, intent(in) :: ndim
real(r64), dimension(ndim,ndim), intent(in) :: A, F, G
real(r64), dimension(ndim,ndim), intent(out) :: C
real(r64), dimension(ndim,ndim) :: B 

call dgemm('t','n',ndim,ndim,ndim,one,F,ndim,A,ndim,zero,B,ndim)
call dgemm('n','n',ndim,ndim,ndim,one,B,ndim,G,ndim,zero,C,ndim)

end subroutine multiply_weights_full

!------------------------------------------------------------------------------!
! subroutine ClebschGordan                                                     !
!                                                                              !
! Computes the ClebschGordan for the group SU(2). The algorithm used was taken !
! from technical notes from NASA written by W. F. Ford and R. C. Bruley.       !
! Ref: NASA TN D-6173                                                          !
!                                                                              !
! (j1,m1,j2,m2|j3,m3) = c * g                                                  !
! with                                                                         !
! c = D(j1j2j3) * [(j1-m1)!(j2+m2)!(j1+m1)!(j2-m2)!(j3+m3)!(j3-m3)!]**1/2      !
! g = sqrt(2*j3+1) sum_l (-1)**l [(j1+j2-j3-l)!(j1-m1-l)!(j2+m2-l)!            !
!                                 (j3-j2+m1+l)!(j3-j1-m1+l)!l!]**-1            !
! D(j1j2j3) = [(j1+j2-j3)!(j2+j3-j1)!(j3+j1-j2)!]**1/2 / [(j1+j2+j3+1)!]**1/2  !
!                                                                              !
! Input: j1,j2,j3 = angular momentum x 2                                       !
!        m1,m2,m3 = z-component x 2
!                                                                              !
! Output: cg = Clebsh-Gordan coefficient                                       !
!------------------------------------------------------------------------------!
subroutine ClebschGordan (j1,j2,j3,m1,m2,m3,cg)

integer, intent(in) :: j1, j2, j3, m1, m2, m3
real(r64), intent(out) :: cg 
integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, k1, k2, k3, k4, k5, k6, &
           l, l1, l2
real(r64) :: p, q, h, hl, hlm1

cg = zero

!!! Computes the ingredients for the factor c
n1 = 1 + (j1 + j2 - j3)/2
n2 = 1 + (j2 + j3 - j1)/2
n3 = 1 + (j3 + j1 - j2)/2
n4 = 1 + (j1 - m1)/2
n5 = 1 + (j2 + m2)/2
n6 = 1 + (j1 + m1)/2
n7 = 1 + (j2 - m2)/2
n8 = 1 + (j3 + m3)/2
n9 = 1 + (j3 - m3)/2
n10= n1 + n2 + n3 - 1

if ( (min(n1,n2,n3,n4,n5,n6,n7,n8,n9) < 1) .or. (m1+m2 /= m3) ) return

p =  log_gamma(n1+zero) + log_gamma(n2+zero) + log_gamma(n3+zero) &
   + log_gamma(n4+zero) + log_gamma(n5+zero) + log_gamma(n6+zero) &
   + log_gamma(n7+zero) + log_gamma(n8+zero) + log_gamma(n9+zero) &
   - log_gamma(n10+zero)

!!! Computes the ingredients for the factor g
k1 = n1
k2 = n4
k3 = n5
k4 = n4 - n3
k5 = n5 - n2

l1 = max(0,k4,k5)
l2 = min(k1,k2,k3)

h  = one 
hl = one

do l = l1+1, l2
  hlm1 = hl
  hl = hlm1 * (l - k1) * (l - k2) * (l - k3) / ((l - k4) * (l - k5) * l)
  h = h + hl
enddo

k1 = k1 - l1
k2 = k2 - l1
k3 = k3 - l1
k4 = l1 + 1 - k4 
k5 = l1 + 1 - k5 
k6 = l1 + 1 

q =  log_gamma(k1+zero) + log_gamma(k2+zero) + log_gamma(k3+zero) &
   + log_gamma(k4+zero) + log_gamma(k5+zero) + log_gamma(k6+zero)

!!! Computes the final value combining the two parts.
cg = sqrt(j3 + one) * (-1)**l1 * exp(0.5d0*p - q) * h

end subroutine ClebschGordan

END MODULE MathMethods 
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
