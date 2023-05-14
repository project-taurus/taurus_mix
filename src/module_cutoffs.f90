!==============================================================================!
! MODULE Cutoffs                                                               !
!                                                                              !
! This module contains the variables and routines related to the cut-offs used !
! when reading the states or solving the HWG equations.                        !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_cutoffs                                                     !
! - subroutine print_cutoffs_blockJP                                           !
!==============================================================================!
MODULE Cutoffs     

use Parameters 

implicit none
public

!!! Global cutoffs
integer :: cutoff_ldim ! Max. number of reference states read (0 = all)
real(r64) :: cutoff_over,  & ! Cutoff on the projected overlap 
             cutoff_negev, & !   "    "  the negative eigenvalues
             cutoff_algo,  & !   "    "  the norm eigenvalues  
             cutoff_J,     & !   "    "  expect. val of J/J^2
             cutoff_A        !   "    "  expect. val of N/Z/A

!!! Cutoffs specific to the quantum number blocks (read in inputs)
integer :: cutoff_spec_dim, & ! Number of specific cutoffs (0 = none)
           cutoff_label_dim   !   "    "     "        "    of type "label"
integer, dimension(:), allocatable :: cutoff_spec_2j, & ! J value of cutoff
                                      cutoff_spec_p, &  ! P value of cutoff 
                                      cutoff_spec_lab2k ! 2*K  of cutoff 
integer(i64), dimension(:), allocatable :: cutoff_spec_label ! Label of cutoff
real(r64), dimension(:), allocatable :: cutoff_spec_value   ! Value of cutoff
character(1), dimension(:), allocatable :: cutoff_spec_type ! Type: O,N,J,A   


!!! Cutoffs specific to the quantum number blocks (same as above but now
!!! assigned to a given block)
integer, dimension(:,:,:), allocatable :: cutoff_block_lab2k    
integer(i64), dimension(:,:,:), allocatable :: cutoff_block_label
real(r64), dimension(:,:), allocatable :: cutoff_block_over,  &
                                          cutoff_block_negev, &
                                          cutoff_block_algo,  &
                                          cutoff_block_J,     &
                                          cutoff_block_A

CONTAINS 

!------------------------------------------------------------------------------!
! subroutine set_cutoffs                                                       !
!                                                                              !
! Sets cutoffs for the different quantumn number blocks. The cutoffs are first !
! set using the global values. Then, if needed, specific values are assigned.  !
!------------------------------------------------------------------------------!
subroutine set_cutoffs

integer :: i, j, ndim, n2j, np, kvalue, ialloc=0 
integer(i64) :: lvalue
real(r64) :: xvalue
character(1) :: ctype 

j = 0
ndim = cutoff_label_dim

allocate ( cutoff_block_over(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax),  &
           cutoff_block_negev(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           cutoff_block_algo(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax),  &
           cutoff_block_J(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           cutoff_block_A(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           cutoff_block_label(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           cutoff_block_lab2k(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of arrays for cutoffs'   

!!! Sets all the cutoffs using the global values 
cutoff_block_over = cutoff_over
cutoff_block_negev= cutoff_negev
cutoff_block_algo = cutoff_algo
cutoff_block_J = cutoff_J
cutoff_block_A = cutoff_A
cutoff_block_label = 0
cutoff_block_lab2k = 999

!!! Sets the specific values 
do i = 1, cutoff_spec_dim
  n2j = cutoff_spec_2j(i)
  np = cutoff_spec_p(i)
  ctype = cutoff_spec_type(i)
  xvalue = cutoff_spec_value(i)
  lvalue = cutoff_spec_label(i)
  kvalue = cutoff_spec_lab2k(i)
 
  if ( ctype == 'O' ) then
    cutoff_block_over(n2j,np) = xvalue
  elseif ( ctype == 'N' ) then
    cutoff_block_negev(n2j,np) = xvalue
  elseif ( ctype == 'S' ) then
    cutoff_block_algo(n2j,np) = xvalue
  elseif ( ctype == 'J' ) then
    cutoff_block_J(n2j,np) = xvalue
  elseif ( ctype == 'A' ) then
    cutoff_block_A(n2j,np) = xvalue
  elseif ( ctype == 'L' ) then
    j = j + 1
    cutoff_block_label(j,n2j,np) = lvalue
    cutoff_block_lab2k(j,n2j,np) = kvalue
  endif
enddo

end subroutine set_cutoffs

!------------------------------------------------------------------------------!
! subroutine print_cutoffs_blockJP                                             !
!                                                                              !
! Prints the cutoffs used to reject the states for the block (2j,p).           !
!                                                                              !
! Input: block_2j = value of 2*j for the symmetry block                        !
!        block_p  =   "   "   p   "   "     "       "                          !
!------------------------------------------------------------------------------!
subroutine print_cutoffs_blockJP(block_2j,block_p)

integer, intent(in) :: block_2j, block_p
integer :: i
character(len=*), parameter :: format1 = "(1a15,1x,1es10.3)", &
                               format2 = "(1a15,1x,1i19,1i5)"

!!! Print the properties
print '(/,"Cutoffs used",/,12("="),//, &
       &5x,"Type",10x,"Value",/,26("-"))'

print format1, 'Proj. overlap  ', cutoff_block_over(block_2j,block_p)
if ( hwg_algo == 0 ) then
  print format1, 'Norm eigenval. ', cutoff_block_algo(block_2j,block_p)
else
  print format1, 'QZ eigenvalues ', cutoff_block_algo(block_2j,block_p)
endif
if ( hwg_rmev == 1 ) then
  print format1, '|neg. eigenv.| ', cutoff_block_negev(block_2j,block_p)
endif
print format1, 'Exp. val. Jz/J ', cutoff_block_J(block_2j,block_p)
print format1, 'Exp. val. Z/N/A', cutoff_block_A(block_2j,block_p)

do i = 1, cutoff_label_dim
  if ( cutoff_block_label(i,block_2j,block_p) /= 0 ) then 
    print format2, 'States (lab,2K)', cutoff_block_label(i,block_2j,block_p), &
                                      cutoff_block_lab2k(i,block_2j,block_p)
  endif
enddo

end subroutine print_cutoffs_blockJP    

END MODULE Cutoffs      
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
