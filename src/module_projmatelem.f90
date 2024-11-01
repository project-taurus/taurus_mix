!==============================================================================!
! MODULE ProjMatElem                                                           !
!                                                                              !
! This module contains the variables and routines related to the projected     !
! matrix elements computed by taurus_pav.                                      !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine read_projme_init                                                !
! - subroutine read_projmatelem_blockJP                                        !
! - subroutine read_projmatelem_states                                         !
! - subroutine print_projmatelem_states                                        !
! - subroutine normalize_projmatelem_states                                    !
! - subroutine remove_negev_projmatelem_states                                 !
! - subroutine read_projmatelem_occnumb                                        !
! - subroutine read_projmatelem_transitions_elm                                !
!==============================================================================!
MODULE ProjMatElem       

use Constants
use Cutoffs   
use Parameters  

implicit none
public

!!! Properties of the states defining the mixing
integer :: projme_tdim ! Number of possibles states (labels * quant. numb.)
integer, dimension(-3:0,-1:1,2) :: projme_ldim, & ! Number of labels read
                                   projme_bdim    ! Number of states for (2j,p)
integer(i64), dimension(:), allocatable :: projme_label_read ! Labels read
integer(i64), dimension(:,:,:,:), allocatable :: projme_label ! Labels stored
integer(i32), dimension(:,:,:,:), allocatable :: projme_2kj   ! 2*kj  
integer(i32), dimension(:,:), allocatable :: projme_read, & ! # mat. elem. read
                                             warnings_zero  ! # mat. elem. = 0

!!! Projected matrix elements for the different quantities
real(r64), dimension(:,:), allocatable :: projme_over, & ! overlap
                                          projme_ener, & ! energy
                                          projme_pari, & ! parity
                                          projme_prot, & ! proton  number
                                          projme_neut, & ! neutron   "
                                          projme_nucl, & ! nucleon   "
                                          projme_amj2, & ! ang. mom. J^2
                                          projme_ist2, & ! isospin T^2
                                          projme_r2so    ! isospin T^2

real(r64), dimension(:,:), allocatable :: projme_0nuF, &  ! Fermi (0n2b)
                                          projme_0nuG, &  ! Gamow-Teller (0n2b)
                                          projme_0nuT     ! Tensor (0n2b)

real(r64), dimension(:,:,:), allocatable :: projme_E1, &  ! E1 transition ELM
                                            projme_E2, &  ! E2      "      "
                                            projme_E3, &  ! E3      "      "
                                            projme_M1, &  ! M1      "      "
                                            projme_M2, &  ! M2      "      "
                                            projme_r2p, & ! radius protons
                                            projme_r2n, & !   "    neutrons
                                            projme_r2m    !   "    matter   

real(r64), dimension(:,:,:,:), allocatable :: projme_occn ! occupation numbers 

real(r64), dimension(:,:,:), allocatable :: projme_norm ! normalization

!!! Local arrays 
real(r64), dimension(:,:,:), allocatable, private :: transition

CONTAINS

!------------------------------------------------------------------------------!
! subroutine read_projme_init                                                  !
!                                                                              !
! Reads the file containing the projected states. The first read is used to    !
! determine the number of states possibles by counting the number of diagonal  !
! matrix elements with different labels.                                       !
!                                                                              !
! Remark: The number of states determined in this routines does not take into  !
! account the various cutoffs. It could be done if memory becomes a problem,   !
! but right now it would just add an unnecessary complexity to the routine.    !
!------------------------------------------------------------------------------!
subroutine read_projmatelem_initial

integer :: i, htype, facn, iexit, ialloc=0, ialloctot=0
integer(i64) :: label_l, label_r 
integer(i64), dimension(10000) :: label_read=0
logical :: is_new

!!! Initialization
projme_bdim = 0
projme_ldim = 0
projme_tdim = 0

!!! Reading the states
do 
  read(utst,iostat=iexit) label_l, label_r
  if ( iexit /= 0 ) exit

  !!! Only considers diagonal states
  if ( label_l /= label_r ) cycle

  !!! Detects if the state is a new one and stores its identity   
  is_new = .true.

  do i = 1, projme_ldim(0,0,1)
    if ( label_r == label_read(i) ) is_new = .false.
  enddo

  if ( is_new ) then 
    projme_ldim(0,0,1) = projme_ldim(0,0,1) + 1
    label_read(projme_ldim(0,0,1)) = label_r
  endif

  !!! Exits if the maximum number of states has been read
  if ( projme_ldim(0,0,1) == cutoff_ldim ) exit
enddo

rewind(utst)

!!! Reading the occupation numbers to get the information about the basis
if ( do_occn ) then 
  read(utoc) i, htype, HOsh_dim
  rewind(utoc)

  allocate( HOsh_n(HOsh_dim), HOsh_l(HOsh_dim), HOsh_2j(HOsh_dim), &
            HOsh_na(HOsh_dim), stat=ialloc )
  if ( ialloc /= 0 ) stop 'Error during allocation of shells'

  read(utoc) i, htype, HOsh_dim, HOsh_na
  rewind(utoc)

  !!! Determines the quantum numbers of the shells
  if ( (htype == 1) .or. (htype == 2) ) then
    facn = 1000
  else
    facn = 10000
  endif

  do i = 1, HOsh_dim
    HOsh_n(i)  =  HOsh_na(i) / facn 
    HOsh_l(i)  = (HOsh_na(i) - HOsh_n(i)*facn) / 100
    HOsh_2j(i) =  HOsh_na(i) - HOsh_n(i)*facn - HOsh_l(i)*100
  enddo
endif

!!! Computes the total number of states possibles given the symmetries
projme_tdim = projme_ldim(0,0,1) * (hwg_2jmax + 1)

!!! Allocations: arrays to identify the states
allocate( projme_label_read(projme_ldim(0,0,1)),              &
          projme_label(projme_tdim,-3:0,hwg_pmin:hwg_pmax,1), &
          projme_2kj(projme_tdim,-3:0,hwg_pmin:hwg_pmax,1),   &
          stat=ialloc )
projme_label_read(1:projme_ldim(0,0,1)) = label_read(1:projme_ldim(0,0,1))
projme_label = 0
projme_2kj = 0

!!! Allocations: expectation values of operators and transitions
! The arrays are initialized at zero when reading the matrix elements
! (the dimensions could be reduced if we encounter memory problems)
allocate( projme_over(projme_tdim,projme_tdim), &
          projme_ener(projme_tdim,projme_tdim), &
          projme_pari(projme_tdim,projme_tdim), &
          projme_prot(projme_tdim,projme_tdim), &
          projme_neut(projme_tdim,projme_tdim), &
          projme_nucl(projme_tdim,projme_tdim), &
          projme_amj2(projme_tdim,projme_tdim), &
          projme_ist2(projme_tdim,projme_tdim), &
          projme_r2so(projme_tdim,projme_tdim), &
          projme_read(projme_tdim,projme_tdim), &
          projme_r2p(projme_tdim,projme_tdim,0:0), &
          projme_r2n(projme_tdim,projme_tdim,0:0), &
          projme_r2m(projme_tdim,projme_tdim,0:0), &
          projme_occn(projme_tdim,projme_tdim,HOsh_dim,2), &
          transition(projme_tdim,projme_tdim,-3:0), &
          warnings_zero(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
          stat=ialloc )

ialloctot = ialloctot + ialloc

warnings_zero = 0

if ( do_E1 ) then 
  allocate( projme_E1(projme_tdim,projme_tdim,-1:0), stat=ialloc )
  ialloctot = ialloctot + ialloc
endif

if ( do_E2 ) then 
  allocate( projme_E2(projme_tdim,projme_tdim,-2:0), stat=ialloc )
  ialloctot = ialloctot + ialloc
endif

if ( do_E3 ) then 
  allocate( projme_E3(projme_tdim,projme_tdim,-3:0), stat=ialloc )
  ialloctot = ialloctot + ialloc
endif

if ( do_M1 ) then 
  allocate( projme_M1(projme_tdim,projme_tdim,-1:0), stat=ialloc )
  ialloctot = ialloctot + ialloc
endif

if ( do_M2 ) then 
  allocate( projme_M2(projme_tdim,projme_tdim,-2:0), stat=ialloc )
  ialloctot = ialloctot + ialloc
endif

if ( hwg_norm == 1 ) then 
  allocate( projme_norm(projme_tdim,hwg_pmin:hwg_pmax,-3:0), &
            stat=ialloc )
  ialloctot = ialloctot + ialloc
  projme_norm = zero
endif

if ( ialloctot /= 0 ) stop 'Error during allocation of arrays for projected &
                            &matrix elements'

end subroutine read_projmatelem_initial

!------------------------------------------------------------------------------!
! subroutine read_projmatelem_blockJP                                          !
!                                                                              !
! Reads the file containing the projected states. The second reading is used   !
! to determine the number of states for a given (2j,p).                        !
!                                                                              !
! Input: block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "    "     "       "                         !
!                                                                              !
! Remark: Need to adapt the cutoff on N,Z,A when we will consider iso. proj.   !
!------------------------------------------------------------------------------!
subroutine read_projmatelem_blockJP(block_2j,block_p)

integer, intent(in) :: block_2j, block_p
integer :: i, icut, cutdisp, n2j, n2mj, n2kj, hkj, np, iexit
integer(i64) :: label_l, label_r
integer, dimension(1) :: loclab
real(r64) :: xover, xener, xpari, xprot, xneut, xnucl, xamj2, xamsp, xamsn, &
             xist2, xr2p, xr2n, xr2m, xj, xt
logical :: is_new_state, is_new_label
character(len=*), parameter :: format1 = "(1i19,1i5,1f12.8,1f12.5,1x,4f12.6,& 
                                          &2f10.5)"

!!! Initialization
cutdisp = 0
projme_bdim(0,block_p,1) = 0
projme_ldim(0,block_p,1) = 0
projme_label(:,0,block_p,1) = 0
projme_2kj(:,0,block_p,1) = 0

!!! Reading            
do
  read(utst,iostat=iexit) label_l, label_r, n2j, n2mj, n2kj, np, &
                          xover, xener, xpari, xprot, xneut, xamj2, xamsp, &
                          xamsn, xist2, xr2p, xr2n, xr2m
  xnucl = xprot + xneut
  if ( iexit /= 0 ) exit

  !!! Only considers diagonal states
  if ( label_l /= label_r ) cycle
  if ( n2mj    /= n2kj    ) cycle
 
  !!! Only considers reference states previously recorded   
  loclab = findloc(projme_label_read,label_l)
  if ( loclab(1) == 0 ) cycle
 
  !!! Only considers the states with good quantum numbers
  if ( n2j /= block_2j ) cycle
  if ( np  /= block_p  ) cycle

  !!! Only considers the states that respect the cutoffs
  icut = 0

  do i = 1, cutoff_label_dim
    if ( (label_r == cutoff_block_label(i,n2j,np)) .and. &
         ((n2kj == cutoff_block_lab2k(i,n2j,np)) .or. &
          (99 == cutoff_block_lab2k(i,n2j,np))) ) icut = 1
  enddo

  xj = 0.5d0 * (-1 + sqrt(1+4*abs(xamj2/xover)))
  xt = 0.5d0 * (-1 + sqrt(1+4*abs(xist2/xover)))
  if ( abs(xover) < cutoff_block_over(n2j,np) ) icut = 1 
  if ( xener/xover > cutoff_block_ener(n2j,np) ) icut = 1 
  if ( abs(xj - 0.5d0*block_2j) > cutoff_block_J(n2j,np) ) icut = 1
  if ( abs(xprot/xover - hwg_Zv) > cutoff_block_A(n2j,np) ) icut = 1
  if ( abs(xneut/xover - hwg_Nv) > cutoff_block_A(n2j,np) ) icut = 1
  !if ( abs(xnucl/xover - hwg_Av) > cutoff_block_A(n2j,np) ) cycle 
 
  if ( icut == 1 ) then
    hkj = n2kj
    if ( (-1)**n2kj == 1 ) hkj = n2kj/2
  
    if ( abs(xover) > zero_eps ) then
      if ( cutdisp == 0 ) then
        cutdisp = 1
        print '(/,"States discarded",/,16("="),//, &
              &7x,"Label",9x,1a3,3x,"Overlap",7x,"Energy",9x,"P",11x,"Z",&
              &11x,"N",11x,"A",10x,"J",9x,"T",/,116("-"))',char_Kleg
      endif
      write(uto,format1) label_r, hkj, xover, xener/xover, xpari/xover, &
                         xprot/xover, xneut/xover, xnucl/xover, xj, xt
    endif
    cycle
  endif

  !!! Detects if the state is a new one and stores its identity   
  is_new_state = .true.
  is_new_label = .true.

  do i = 1, projme_bdim(0,block_p,1)
    if ( label_r == projme_label(i,0,block_p,1) ) then
      is_new_label = .false.
      if ( n2kj == projme_2kj(i,0,block_p,1) ) then  
        is_new_state = .false.
        exit
      endif
    endif
  enddo
  
  if ( is_new_label ) then 
    projme_ldim(0,block_p,1) = projme_ldim(0,block_p,1) + 1
  endif
  
  if ( is_new_state ) then 
    projme_bdim(0,block_p,1) = projme_bdim(0,block_p,1) + 1
    projme_label(projme_bdim(0,block_p,1),0,block_p,1) = label_r
    projme_2kj(projme_bdim(0,block_p,1),0,block_p,1) = n2kj   
  endif
enddo

rewind(utst)

end subroutine read_projmatelem_blockJP   

!------------------------------------------------------------------------------!
! subroutine read_projmatelem_states                                           !
!                                                                              !
! Reads the file containing the expectation values of various operators. The   !
! third reading is used to store the information on the states.                !
!                                                                              !
! Remark: could probably be merged with read_projmatelem_blockJP.              !
!                                                                              !
! Input: block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "    "     "       "                         !
!------------------------------------------------------------------------------!
subroutine read_projmatelem_states(block_2j,block_p)

integer, intent(in) :: block_2j, block_p
integer :: i, n2j, n2mj, n2kj, np, idx_l, idx_r, iexit
integer(i64) :: label_l, label_r
integer, dimension(1) :: loclab
real(r64) :: xover, xener, xpari, xprot, xneut, xnucl, xamj2, xamsp, xamsn, &
             xist2, xr2p, xr2n, xr2m, xr2so

!!! Initialization
projme_read = 0
projme_over = zero
projme_ener = zero
projme_pari = zero 
projme_prot = zero
projme_neut = zero
projme_nucl = zero
projme_amj2 = zero
projme_ist2 = zero 
projme_r2p  = zero
projme_r2n  = zero
projme_r2m  = zero
projme_r2so = zero

!!! Reading            
do 
  read(utst, iostat=iexit) label_l, label_r, n2j, n2mj, n2kj, np, &
                           xover, xener, xpari, xprot, xneut, xamj2, &
                           xamsp, xamsn, xist2, xr2p, xr2n, xr2m
  xnucl = xprot + xneut
  if ( iexit /= 0 ) exit

  !!! Only considers reference states previously recorded   
  loclab = findloc(projme_label_read,label_l)
  if ( loclab(1) == 0 ) cycle 

  loclab = findloc(projme_label_read,label_r)
  if ( loclab(1) == 0 ) cycle 
 
  !!! Only considers the states with good quantum numbers
  if ( n2j /= block_2j ) cycle
  if ( np  /= block_p  ) cycle
 
  !!! Determines the presence/index of the left state 
  idx_l = 0 
  do i = 1, projme_bdim(0,block_p,1)
    if ( label_l == projme_label(i,0,block_p,1) ) then
      if ( n2mj == projme_2kj(i,0,block_p,1) ) then
        idx_l = i
        exit
      endif
    endif
  enddo

  if ( idx_l == 0 ) cycle

  !!! Determines the presence/index of the right state
  idx_r = 0 
  do i = 1, projme_bdim(0,block_p,1)
    if ( label_r == projme_label(i,0,block_p,1) ) then
      if ( n2kj == projme_2kj(i,0,block_p,1) ) then
        idx_r = i
        exit
      endif
    endif
  enddo

  if ( idx_r == 0 ) cycle

  !!! Computes the spin-orbit correction for the charge radii
  xr2so = (magmome_mup - 0.5d0) * xamsp + magmome_mun * xamsn

  !!! Stores the matrix elements (direct and symmetric)
  projme_over(idx_l,idx_r)  = xover
  projme_ener(idx_l,idx_r)  = xener
  projme_pari(idx_l,idx_r)  = xpari
  projme_prot(idx_l,idx_r)  = xprot
  projme_neut(idx_l,idx_r)  = xneut
  projme_nucl(idx_l,idx_r)  = xnucl
  projme_amj2(idx_l,idx_r)  = xamj2
  projme_ist2(idx_l,idx_r)  = xist2
  projme_r2p(idx_l,idx_r,0) = xr2p   
  projme_r2n(idx_l,idx_r,0) = xr2n
  projme_r2m(idx_l,idx_r,0) = xr2m
  projme_r2so(idx_l,idx_r)  = xr2so

  projme_over(idx_r,idx_l)  = xover
  projme_ener(idx_r,idx_l)  = xener
  projme_pari(idx_r,idx_l)  = xpari
  projme_prot(idx_r,idx_l)  = xprot
  projme_neut(idx_r,idx_l)  = xneut
  projme_nucl(idx_r,idx_l)  = xnucl
  projme_amj2(idx_r,idx_l)  = xamj2
  projme_ist2(idx_r,idx_l)  = xist2
  projme_r2p(idx_r,idx_l,0) = xr2p
  projme_r2n(idx_r,idx_l,0) = xr2n
  projme_r2m(idx_r,idx_l,0) = xr2m
  projme_r2so(idx_r,idx_l)  = xr2so

  !!! Indicates that the matrix element has been read
  projme_read(idx_l,idx_r) = 1
  projme_read(idx_r,idx_l) = 1
enddo

!!! Stores the number of missing matrix elements for later printing
warnings_zero(block_2j,block_p) = projme_bdim(0,block_p,1)**2 - sum(projme_read) 

rewind(utst)

end subroutine read_projmatelem_states    

!------------------------------------------------------------------------------!
! subroutine print_projmatelem_states                                          !
!                                                                              !
! Prints the properties of the states read and included for this quantum numb- !
! er block. The quantities printed are normalized (division by the overlap).   !
!                                                                              !
! Input: block_p = values of p for this symmetry block                         !
!------------------------------------------------------------------------------!
subroutine print_projmatelem_states(block_p)

integer, intent(in) :: block_p
integer :: i, j, kj, hkj, bdim
integer(i64) :: label 
real(r64) :: xover, xener, xpari, xprot, xneut, xnucl, xamj2, xist2, xj, xt
character(len=*), parameter :: format1 = "(1i5,1x,1i19,1i5,1f12.8,1x,1f12.5, &
                                          &4f12.6,2f10.5)", &
                               format2 = "(1i5,3x,1i5)"

print '(/,"States selected",/,15("="),//, &
       &"Index",8x,"Label",9x,1a3,3x,"Overlap",7x,"Energy",9x,"P",11x,"Z", &
       &11x,"N",11x,"A",10x,"J",9x,"T",/,123("-"))',char_Kleg

do i = 1, projme_bdim(0,block_p,1)
  label = projme_label(i,0,block_p,1)
  kj = projme_2kj(i,0,block_p,1)
  
  xover = projme_over(i,i)
  xener = projme_ener(i,i) / xover
  xpari = projme_pari(i,i) / xover
  xprot = projme_prot(i,i) / xover
  xneut = projme_neut(i,i) / xover
  xnucl = projme_nucl(i,i) / xover
  xamj2 = projme_amj2(i,i) / xover
  xist2 = projme_ist2(i,i) / xover
  
  xj = 0.5d0 * (-1 + sqrt(1+4*abs(xamj2)))
  xt = 0.5d0 * (-1 + sqrt(1+4*abs(xist2)))
 
  hkj = kj 
  if ( (-1)**kj == 1 ) hkj = kj/2
  
  write(uto,format1) i, label, hkj, xover, xener, xpari, xprot, xneut, &
                     xnucl, xj, xt
enddo

!!! Checks if all matrix elements have been read. If no, prints a warning array
!!! indicating the number of missing element.
bdim = projme_bdim(0,block_p,1)

if ( sum(projme_read(1:bdim,1:bdim)) /= bdim**2 ) then 
  print '(/,"Warning: missing matrix elements",/,32("="),//, &
          &"Label",2x,"# missing",/,16("-"))'
  do i = 1, bdim
    j = sum(projme_read(1:bdim,i))
    if ( j /= bdim ) write(uto,format2) i, bdim - j
  enddo
endif

end subroutine print_projmatelem_states    

!------------------------------------------------------------------------------!
! subroutine normalize_projmatelem_states                                      !
!                                                                              !
! Normalizes the matrix elements of the states. This is optional but this may  !
! help to have a better-conditioned general eigenvalue problem.                !
!                                                                              !
! Input: block_p = values of p for this symmetry block                         !
!------------------------------------------------------------------------------!
subroutine normalize_projmatelem_states(block_p)

integer, intent(in) :: block_p
integer :: i, j
real(r64) :: factor

!!! Computes the normalization factor for each state 
do i = 1, projme_bdim(0,block_p,1)
  projme_norm(i,block_p,0) = 1 / sqrt(projme_over(i,i))
enddo

!!! Normalizes all the quantities
do j = 1, projme_bdim(0,block_p,1)
  do i = 1, projme_bdim(0,block_p,1)
    factor = projme_norm(i,block_p,0) * projme_norm(j,block_p,0)
    projme_over(i,j) = projme_over(i,j) * factor
    projme_ener(i,j) = projme_ener(i,j) * factor
    projme_pari(i,j) = projme_pari(i,j) * factor
    projme_prot(i,j) = projme_prot(i,j) * factor
    projme_neut(i,j) = projme_neut(i,j) * factor
    projme_nucl(i,j) = projme_nucl(i,j) * factor
    projme_amj2(i,j) = projme_amj2(i,j) * factor
    projme_ist2(i,j) = projme_ist2(i,j) * factor
    projme_r2so(i,j) = projme_r2so(i,j) * factor
    projme_r2p(i,j,0) = projme_r2p(i,j,0) * factor
    projme_r2n(i,j,0) = projme_r2n(i,j,0) * factor
    projme_r2m(i,j,0) = projme_r2m(i,j,0) * factor
  enddo
enddo

end subroutine normalize_projmatelem_states    

!------------------------------------------------------------------------------!
! subroutine remove_negev_projmatelem_states                                   !
!                                                                              !
! Removes the projected states that give large negative norm eigenvalues.      !
! Multiple norm diagonalizations are performed, incrementally adding the       !
! projected states.                                                            !
!                                                                              !
! Input: block_2j = value of 2j of the current symmetry block                  !
!        block_p  = value of p of the current symmetry block                   !
!------------------------------------------------------------------------------!
subroutine remove_negev_projmatelem_states(block_2j,block_p)

integer, intent(in) :: block_2j, block_p
integer :: i, j, k, bdim, zdim, rdim, ialloc, info
integer, dimension(:), allocatable :: list_rm, list_zdim
real(r64) :: cutoff
real(r64), dimension(:), allocatable :: eigen_norm, work
real(r64), dimension(:,:), allocatable :: overlap
character(len=*), parameter :: format1 = "(1i5,3x,1i5)"

bdim = projme_bdim(0,block_p,1)
cutoff = cutoff_block_negev(block_2j,block_p)

allocate( list_rm(bdim), list_zdim(bdim), overlap(bdim,bdim), & 
          eigen_norm(bdim), work(3*bdim-1), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation when removing the states'

!!! Performs the successive diagonalization
rdim = 0
list_rm = 0
list_zdim = 0

do i = 1, bdim
  zdim = 0
  overlap = zero            
  overlap(1:i,1:i) = projme_over(1:i,1:i)

  do j = 1, rdim
    k = list_rm(j)
    overlap(k,:) = zero
    overlap(:,k) = zero
  enddo

  call dsyev('V','U',bdim,overlap,bdim,eigen_norm,work,3*bdim-1,info)
 
  do j = 1, bdim
    if ( eigen_norm(j) < -cutoff ) zdim = zdim + 1
  enddo

  if ( zdim > 0 ) then
    rdim = rdim + 1
    list_rm(rdim) = i
    list_zdim(rdim) = zdim
  endif
enddo

!!! Sets the matrix elements to zero
! Remark: for sqrt algorithm, the F will be zero with that so it works.
! Remark: for QZ algorithm all oberversables (including transitions) have
! to be set to zero. Need to store the index of removed states in an array
! that can be read by read_projmatelem_transitions_elm (not implemented yet).
do j = 1, rdim
  k = list_rm(j)
  projme_over(:,k) = zero
  projme_ener(:,k) = zero
! projme_pari(:,k) = zero
! projme_prot(:,k) = zero
! projme_neut(:,k) = zero
! projme_nucl(:,k) = zero 
! projme_amj2(:,k) = zero
! projme_ist2(:,k) = zero
! projme_r2so(:,k) = zero 
! projme_r2p(:,k,0)= zero
! projme_r2n(:,k,0)= zero
! projme_r2m(:,k,0)= zero

  projme_over(k,:) = zero
  projme_ener(k,:) = zero
! projme_pari(k,:) = zero
! projme_prot(k,:) = zero
! projme_neut(k,:) = zero
! projme_nucl(k,:) = zero 
! projme_amj2(k,:) = zero
! projme_ist2(k,:) = zero
! projme_r2so(k,:) = zero 
! projme_r2p(k,:,0)= zero
! projme_r2n(k,:,0)= zero
! projme_r2m(k,:,0)= zero
enddo

!!! Prints the states removed
if ( rdim > 0 ) then 
  print '(/,"Removing problematic states",/,27("="),//, &
          &"Cutoff value: ",1es10.3,//, &
          &"Label",2x,"# ev < 0",/,15("-"))', -cutoff
  do i = 1, rdim
    write(uto,format1) list_rm(i), list_zdim(i)
  enddo
endif

deallocate( list_rm, list_zdim, overlap, eigen_norm, work )

end subroutine remove_negev_projmatelem_states    

!------------------------------------------------------------------------------!
! subroutine read_projmatelem_occnumb                                          !
!                                                                              !
! Reads the file containing the occupation numbers.                            !
!                                                                              !
! Remark: using pdim, we could parse much more rapidly the file, as we know    !
! how many blocks to jump to go to the next label!                             !
!                                                                              !
! Input: block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "    "     "       "                         !
!------------------------------------------------------------------------------!
subroutine read_projmatelem_occnumb(block_2j,block_p)

integer, intent(in) :: block_2j, block_p
integer :: i, k, n2j, n2mj, n2kj, np, idx_l, idx_r, pdim, iexit, ialloc=0
integer(i64) :: label_l, label_r
integer, dimension(1) :: loclab
real(r64), dimension(:,:), allocatable  :: xoccn
real(r64) :: factor

!!! Initialization
projme_occn = zero

allocate( xoccn(HOsh_dim,2), source=zero, stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation when reading occ. numb.'

!!! Reading            
do 
  read(utoc, iostat=iexit) pdim
  if ( iexit /= 0 ) exit
                          
  do k = 1, pdim
    read(utoc, iostat=iexit) label_l, label_r, n2j, n2mj, n2kj, np, &
                             xoccn

    !!! Only considers reference states previously recorded   
    loclab = findloc(projme_label_read,label_l)
    if ( loclab(1) == 0 ) cycle 
  
    loclab = findloc(projme_label_read,label_r)
    if ( loclab(1) == 0 ) cycle 
   
    !!! Only considers the states with good quantum numbers
    if ( n2j /= block_2j ) cycle
    if ( np  /= block_p  ) cycle
   
    !!! Determines the presence/index of the left state 
    idx_l = 0 
    do i = 1, projme_bdim(0,block_p,1)
      if ( label_l == projme_label(i,0,block_p,1) ) then
        if ( n2mj == projme_2kj(i,0,block_p,1) ) then
          idx_l = i
          exit
        endif
      endif
    enddo
  
    if ( idx_l == 0 ) cycle
  
    !!! Determines the presence/index of the right state
    idx_r = 0 
    do i = 1, projme_bdim(0,block_p,1)
      if ( label_r == projme_label(i,0,block_p,1) ) then
        if ( n2kj == projme_2kj(i,0,block_p,1) ) then
          idx_r = i
          exit
        endif
      endif
    enddo
  
    if ( idx_r == 0 ) cycle
 
    !!! Stores the matrix elements (direct and symmetric)
    if ( hwg_norm == 1 ) then 
      factor = projme_norm(idx_l,block_p,0) * projme_norm(idx_r,block_p,0)
    else
      factor = one
    endif
  
    projme_occn(idx_l,idx_r,:,:) = factor * xoccn(:,:)
    projme_occn(idx_r,idx_l,:,:) = factor * xoccn(:,:)
  enddo
enddo

deallocate( xoccn )

rewind(utoc)

end subroutine read_projmatelem_occnumb   

!------------------------------------------------------------------------------!
! subroutine read_projmatelem_transitions_elm                                  !
!                                                                              !
! Reads the file containing the transition matrix elements for various electro-!
! magnetic transitions.                                                        !
! The files need to have been detected during the initialization.              !
!                                                                              !
! Input: block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "    "     "       "                         !
!------------------------------------------------------------------------------!
subroutine read_projmatelem_transitions_elm(block_2j,block_p)

integer, intent(in) :: block_2j, block_p
integer :: i, jf, kjf, tf, ktf, pf, ji, kji, ti, kti, pi, &
           m, j1, kj1, t1, kt1, p1, j2, kj2, t2, kt2, p2, &
           lt, idx_i, idx_f, dj, utbl, iexit
integer(i64) :: label_f, label_i, label_1, label_2
integer, dimension(1) :: loclab
real(r64) :: factor, factor2, xred, xredp, xredn, xreda, echn
logical :: do_Tl

!!! B(Tl) transitions
do i = 1, 5

  select case (i)
    !!! E1
    case (1)
      if ( .not. do_E1 ) cycle
      do_Tl = do_E1
      utbl = ute1
      lt = -1
      projme_E1 = zero
 
    !!! E2
    case (2)
      if ( .not. do_E2 ) cycle
      do_Tl = do_E2
      utbl = ute2
      lt = -2
      projme_E2 = zero
 
    !!! E3
    case (3)
      if ( .not. do_E3 ) cycle
      do_Tl = do_E3
      utbl = ute3
      lt = -3
      projme_E3 = zero
 
    !!! M1
    case (4)
      if ( .not. do_M1 ) cycle
      do_Tl = do_M1
      utbl = utm1
      lt = -1
      projme_M1 = zero
 
    !!! M2
    case (5)
      if ( .not. do_M2 ) cycle
      do_Tl = do_M2
      utbl = utm2
      lt = -2
      projme_M2 = zero
  end select

  !!! Reads the transitions
  if ( do_Tl ) then 
    transition = zero
    do 
      if ( i <= 3 ) then 
        read(utbl,iostat=iexit) label_1, label_2, j1, kj1, p1, j2, kj2, p2, &
                                xredp, xredn
    
        ! Effective charges: general case
        if ( i /= 1 ) then
          xred = hwg_echp * xredp + hwg_echn * xredn
         
        ! E1 caser: special attention to com contamination,
        ! see appendix B in Ring and Schuck
        elseif ( i == 1 ) then
     
          if ( hwg_Ac == 0 ) then 
            echn = hwg_echp
          else 
            echn = hwg_echn
          endif

          xred =  (hwg_Nv * hwg_echp / hwg_Av) * xredp &
                - (hwg_Zv *     echn / hwg_Av) * xredn
        endif

      else                 
        read(utbl,iostat=iexit) label_1, label_2, j1, kj1, p1, j2, kj2, p2, &
                                xreda
        xred = xreda
      endif                
      if ( iexit /= 0 ) exit

      !!! Only considers reference states previously recorded   
      loclab = findloc(projme_label_read,label_1)
      if ( loclab(1) == 0 ) cycle 
    
      loclab = findloc(projme_label_read,label_2)
      if ( loclab(1) == 0 ) cycle 

      !!! Only considers the states with good quantum numbers compatible with
      !!! the order of the loop on j,p
      if ( (j1 < hwg_2jmin) .or. (j1 > hwg_2jmax) ) cycle 
      if ( (j2 < hwg_2jmin) .or. (j2 > hwg_2jmax) ) cycle 
      if ( (j2 >= j1) .and. ((j2 /= block_2j) .or. (p2 /= block_p)) ) cycle 
      if ( (j1 >= j2) .and. ((j1 /= block_2j) .or. (p1 /= block_p)) ) cycle 
      if ( (j1 == j2) .and. (p1*p2 == -1) .and. (block_p == 1) ) cycle

      dj = abs(j1 - j2)/2

      !!! Reconstructs the correct matrix element
      if ( j2 >= j1 ) then 
        label_f = label_1
        label_i = label_2
        jf = j1
        ji = j2
        kjf = kj1
        kji = kj2
        tf = t1
        ti = t2
        ktf = kt1
        kti = kt2
        pf = p1
        pi = p2
        factor = one
      else
        label_f = label_2
        label_i = label_1
        jf = j2
        ji = j1
        kjf = kj2
        kji = kj1
        tf = t2
        ti = t1
        ktf = kt2
        kti = kt1
        pf = p2
        pi = p1
        factor = (-one)**( (jf-ji)/2 )
      endif

      !!! Determines the presence/index of the final state 
      idx_f = 0 
      do m = 1, projme_bdim(-dj,pf,1)
        if ( label_f == projme_label(m,-dj,pf,1) ) then
          if ( kjf == projme_2kj(m,-dj,pf,1) ) then
            idx_f = m
            exit
          endif
        endif
      enddo
      
      if ( idx_f == 0 ) cycle
      
      !!! Determines the presence/index of the initial state
      idx_i = 0 
      do m = 1, projme_bdim(0,pi,1)
        if ( label_i == projme_label(m,0,pi,1) ) then
          if ( kji == projme_2kj(m,0,pi,1) ) then
            idx_i = m
            exit
          endif 
        endif
      enddo
      
      if ( idx_i == 0 ) cycle

      !!! Stores the matrix element 
      if ( hwg_norm == 1 ) then 
        factor2 = projme_norm(idx_f,pf,-dj) * projme_norm(idx_i,pi,0)
      else
        factor2 = one
      endif

      transition(idx_f,idx_i,-dj) = factor * factor2 * xred   

      if ( dj == 0 ) then 
        transition(idx_i,idx_f,-dj) = factor * factor2 * xred   
      endif
    enddo
   
    !!! Copies the transitions in the correct array 
    select case (i)
      case (1)
        projme_E1(:,:,lt:0) = transition(:,:,lt:0)
      case (2)
        projme_E2(:,:,lt:0) = transition(:,:,lt:0)
      case (3)
        projme_E3(:,:,lt:0) = transition(:,:,lt:0)
      case (4)
        projme_M1(:,:,lt:0) = transition(:,:,lt:0)
      case (5)
        projme_M2(:,:,lt:0) = transition(:,:,lt:0)
    end select
 
    rewind(utbl)
  endif
enddo

end subroutine read_projmatelem_transitions_elm

END MODULE ProjMatElem       
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
