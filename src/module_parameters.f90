!==============================================================================!
! MODULE Parameters                                                            !
!                                                                              !
! This module contains the main parameters used to define the calculations and !
! run the code.                                                                !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_characters_JP                                               !
! - subroutine set_quantnumb_index (inactive)                                  !
!==============================================================================!
MODULE Parameters    

use Constants

implicit none
public

!!! General parameters of the HWG equation
integer :: hwg_phys,   & ! Type of physics (0=spectroscopy)      
           hwg_algo,   & ! Option to choose the the algorithm for HWG equation
           hwg_norm,   & !    "   "  normalize the states before mixing
           hwg_rmev,   & !    "   "  remove states giving large negative ev
           hwg_conv,   & !    "   "  perform convergence analysis (norm ev)
           hwg_read,   & !    "   "  read the weights from a file
           hwg_Z,      & ! Number of active protons Z                               
           hwg_N,      & !    "   "    "    neutrons N
           hwg_A,      & !    "   "    "    nucleons A
           hwg_Zc,     & ! Number of core protons Z                               
           hwg_Nc,     & !    "   "   "   neutrons N
           hwg_Ac,     & !    "   "   "   nucleons A
           hwg_2jmin,  & ! Minimum value for twice the angular momentum 2*J
           hwg_2jmax,  & ! Maximum   "    "    "    "    "        "      "
           hwg_pmin,   & ! Minimum value for the parity P
           hwg_pmax,   & ! Maximum   "    "   "    "    "
           hwg_Edis      ! Maximum excitation energy display in the spectrum

real(r64) :: hwg_echp, & ! Electric charge for protons     
             hwg_echn    !    "        "    "  neutrons   

!!! Basis information
integer :: HOsh_dim ! dimension of the basis
integer, dimension(:), allocatable :: HOsh_n,  & ! quantum number  n
                                      HOsh_l,  & !     "      "    l
                                      HOsh_2j, & !     "      "   2j
                                      HOsh_na    ! label of the shell

!!! Quantum numbers written as nice characters
character(3) :: char_Jleg, char_Kleg 
character(1), dimension(:), allocatable :: char_P
character(5), dimension(:), allocatable :: char_J
character(7), dimension(:,:), allocatable :: char_JP

!!! Switch to read/calculate the different observables
logical :: do_Tl=.false.,  &
           do_E1=.false.,  &
           do_E2=.false.,  &
           do_E3=.false.,  &
           do_M1=.false.,  &
           do_M2=.false.,  &
           do_occn=.false., &
           do_0nu=.false.

!!! Index system (deactivated)
!integer, dimension(:,:), allocatable :: index_2kj, &
!                                        index_2jp

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_characters_JP                                                 !
!                                                                              !
! Sets character arrays containing nice strings that can be used when printing !
! the results. For example 2j=13, p=1 => "13/2^+".                             !
!------------------------------------------------------------------------------!
subroutine set_characters_JP

integer :: j, p, ialloc=0 
character(3) :: temp1

allocate ( char_J(hwg_2jmin:hwg_2jmax), char_P(hwg_pmin:hwg_pmax), &
           char_JP(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of arrays for the characters'   

!!! Angular momentum
do j = hwg_2jmin, hwg_2jmax, 2
  if ( (-1)**j == -1 ) then
    write(temp1,'(1i3)') j
    char_J(j) = adjustr(temp1) // '/2'
  else
    write(temp1,'(1i3)') j/2
    char_J(j) = '  ' // adjustr(temp1) 
  endif
enddo

!!! Parity
do p = hwg_pmin, hwg_pmax, 2
  if ( p == -1 ) then 
    char_P(p) = '-'
  else
    char_P(p) = '+'
  endif
enddo

!!! Angular momentum + parity
do j = hwg_2jmin, hwg_2jmax, 2
  do p = hwg_pmin, hwg_pmax, 2
    char_JP(j,p) = adjustl(char_J(j) // '^' // char_P(p))
  enddo
enddo

!!! Angular momentum J and K for legend of tables
if ( (1)**hwg_2jmin == -1 ) then
  char_Jleg = '2*J'
  char_Kleg = '2*K'
else
  char_Jleg = '  J'
  char_Kleg = '  K'
endif

end subroutine set_characters_JP

!------------------------------------------------------------------------------!
! subroutine set_quantnumb_index                                               !
!                                                                              !
! This routine is not used right now. It can be reactivated if the memory      !
! becomes an issue in the future. (Probably would need some changes).          !
!                                                                              !
! Sets an index system for the quantum bumber to reduce the storage in memory. !
! For example:  HwG_2jmin=2, hwg_2jmax=6, hwg_pmin=-1, hwg_pmax=+1             !
!   2j   p     index                                                           !
!    2  +1  =>   1                                                             !
!    4  +1  =>   2                                                             !
!    6  +1  =>   3                                                             !
!    2  -1  =>  -1                                                             !
!    4  -1  =>  -2                                                             !
!    6  -1  =>  -3                                                             !
!------------------------------------------------------------------------------!
!subroutine set_quantnumb_index
!
!integer :: m, j, k, p, ialloc=0
!
!!!! Computes the numbers of 2*J and P values considered
!hwg_2jdim = (hwg_2jmax - hwg_2jmin)/2 + 1
!hwg_pdim  = (hwg_pmax  - hwg_pmin )/2 + 1
!
!hwg_2jpmin = min(0,hwg_pmin) * hwg_2jdim
!hwg_2jpmax = max(0,hwg_pmax) * hwg_2jdim
!
!!!! Computes the index system
!allocate ( index_2jp(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
!           index_2kj(-hwg_2jmax:hwg_2jmax,0:1), &
!           stat=ialloc )
!if ( ialloc /= 0 ) stop 'Error during allocation of arrays for the indices'   
!
!m = (hwg_2jmin+1)/2
!
!! Index 2jp
!do j = hwg_2jmin, hwg_2jmax, 2
!  do p = hwg_pmin, hwg_pmax, 2
!    index_2jp(j,p) = p *  (1 - m + (j+1)/2)
!  enddo
!enddo
!
!! Index 2k
!index_2kj(:,0) = 1
!
!do k = -hwg_2jmax, hwg_2jmax, 2
!  index_2kj(k,1) = sign(1,k) * (abs(k)+1)/2
!enddo
!
!end subroutine set_quantnumb_index

END MODULE Parameters     
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
