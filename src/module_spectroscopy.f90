!==============================================================================!
! MODULE Spectroscopy                                                          !
!                                                                              !
! This module contains the variables and routines related to the energy        !
! spectrum obtained after solving the HWG equations, i.e. PGCM states and      !
! gamma spectroscopy.                                                          !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine set_spectroscopy                                                !
! - subroutine read_weights_initial                                            !
! - subroutine read_weights_blockJP                                            !
! - subroutine solve_hwg_sqrt                                                  !
! - subroutine solve_hwg_qz                                                    !
! - subroutine calculate_properties_spectrum                                   !
! - subroutine calculate_occnumb_spectrum                                      !
! - subroutine calculate_transitions_elm                                       !
! - subroutine store_solution                                                  !
! - subroutine print_warnings                                                  !
! - subroutine print_complementary_files                                       !
! - subroutine print_spectrum                                                  !
! - subroutine print_transitions_elm                                           !
!==============================================================================!
MODULE Spectroscopy          

use MathMethods
use ProjMatElem   

implicit none
public

!!! Weights solution of the HWG equation
real(r64), dimension(:,:,:,:,:), allocatable :: weights_f 

!!! Properties of the GCM states 
real(r64), dimension(:,:,:), allocatable :: states_ener, & ! energy
                                            states_pari, & ! parity
                                            states_prot, & ! proton  number
                                            states_neut, & ! neutron   "
                                            states_nucl, & ! nucleon   "
                                            states_amj2, & ! ang. mom. J^2
                                            states_ist2, & ! isospin T^2
                                            states_r2so, & ! spin-orbit correc.
                                            states_occn    ! occupation numbers

!!! Reduced matrix elements for electromagnetic transitions
real(r64), dimension(:,:,:,:,:), allocatable :: transi_E1, & ! E1 transition ELM
                                                transi_E2, & ! E2      "      "
                                                transi_E3, & ! E3      "      "
                                                transi_M1, & ! M1      "      "
                                                transi_M2, & ! M2      "      "
                                                transi_r2p, & ! radius protons
                                                transi_r2n, & !   "    neutrons
                                                transi_r2m    !   "    matter  

!!! Various informations
integer :: states_odim, states_imax
integer, dimension(:,:), allocatable :: states_tabijp, &
                                        states_cdim, &
                                        warnings_norm, &
                                        warnings_cplx
integer, dimension(:,:,:), allocatable :: states_notzero

CONTAINS

!------------------------------------------------------------------------------!
! subroutine set_spectrum                                                      !
!                                                                              !
! Sets the arrays containing the energies and weights solution of the HWG equ- !
! ation as well as the expectation values of correlated states.                !
!------------------------------------------------------------------------------!
subroutine set_spectroscopy

integer :: ndim, ialloc=0, ialloctot=0

ndim = projme_tdim

!!! Sets the properties of the GCM states
allocate ( states_ener(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_pari(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_prot(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_neut(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_nucl(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_amj2(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_ist2(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_r2so(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_occn(ndim,HOsh_dim,2), &
           states_notzero(ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_cdim(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           warnings_norm(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           warnings_cplx(hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           weights_f(ndim,ndim,-3:0,hwg_pmin:hwg_pmax,1), &
           stat=ialloc )

ialloctot = ialloctot + ialloc

states_pari = zero
states_prot = zero
states_neut = zero
states_nucl = zero
states_amj2 = zero
states_ist2 = zero
states_r2so = zero
states_occn = zero
states_notzero = 0      
states_cdim = 0      
warnings_norm = 0
warnings_cplx = 0
weights_f = zero

!!! Sets the electromanetic transitions
allocate ( transi_r2p(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-0:0), &
           transi_r2n(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-0:0), &
           transi_r2m(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-0:0), &
           stat=ialloc )
ialloctot = ialloctot + ialloc
transi_r2p = zero
transi_r2n = zero
transi_r2m = zero
  
if ( do_E1 ) then
  allocate ( transi_E1(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-1:0), &
             stat=ialloc )
  ialloctot = ialloctot + ialloc
  transi_E1 = zero
endif

if ( do_E2 ) then
  allocate ( transi_E2(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-2:0), &
             stat=ialloc )
  ialloctot = ialloctot + ialloc
  transi_E2 = zero
endif

if ( do_E3 ) then
  allocate ( transi_E3(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-3:0), &
             stat=ialloc )
  ialloctot = ialloctot + ialloc
  transi_E3 = zero
endif

if ( do_M1 ) then
  allocate ( transi_M1(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-1:0), &
             stat=ialloc )
  ialloctot = ialloctot + ialloc
  transi_M1 = zero
endif

if ( do_M2 ) then
  allocate ( transi_M2(ndim,ndim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax,-2:0), &
             stat=ialloc )
  ialloctot = ialloctot + ialloc
  transi_M2 = zero
endif

if ( ialloctot /= 0 ) stop 'Error during allocation of arrays for the &
                           &spectroscopy' 
  
end subroutine set_spectroscopy

!------------------------------------------------------------------------------!
! subroutine read_weights_initial                                              !
!                                                                              !
! Reads the information from the files containing the weights ("X_weights.txt" !
! with X = left or right) to intialize the problem.                            !
!                                                                              !
! Remark: this routine is not used in the present version of the code.         !
!------------------------------------------------------------------------------!
subroutine read_weights_initial

integer :: ltdim, l2jmin, l2jmax, lpmin, lpmax, lZ, lN, ierr=0, &
           rtdim, r2jmin, r2jmax, rpmin, rpmax, rZ, rN, ialloc=0  
character(len=*), parameter :: format1 = '(7i5)', &
                               format2 = '(1a16,4x,1i3)'

!!! Reads the right weights
read(utrw,format1) rtdim, r2jmin, r2jmax, rpmin, rpmax, rZ, rN
rewind(utrw)

projme_tdim = rtdim

hwg_2jmin = max(hwg_2jmin,r2jmin)
hwg_2jmax = min(hwg_2jmax,r2jmax)

hwg_pmin = max(hwg_pmin,rpmin)
hwg_pmax = min(hwg_pmax,rpmax)

!!! Reads the left weights
if ( hwg_phys == 0 ) then
  ltdim = rtdim
  l2jmin = r2jmin
  l2jmax = r2jmax
  lpmin = rpmin
  lpmax = rpmax
  lZ = rZ   
  lN = rN
endif

projme_tdim = max(projme_tdim,ltdim)

!!! Prints the information
print '(/,60("%"),/,25x,"WEIGHTS F",/,60("%"),//, &
       &"Right weights",5x,"Value",/,23("-"))'
print format2, 'No. of protons  Z', rZ    
print format2, 'No. of neutrons N', rN
print format2, 'Minimum 2*J      ', r2jmin
print format2, 'Maximum 2*J      ', r2jmax
print format2, 'Minimum P        ', rpmin
print format2, 'Minimum P        ', rpmax

print '(/,"Left weights",6x,"Value",/,23("-"))'
print format2, 'No. of protons  Z', lZ    
print format2, 'No. of neutrons N', lN
print format2, 'Minimum 2*J      ', l2jmin
print format2, 'Maximum 2*J      ', l2jmax
print format2, 'Minimum P        ', lpmin
print format2, 'Minimum P        ', lpmax

!!! Checks for possible errors
if ( (r2jmin > hwg_2jmax) .or. (r2jmax < hwg_2jmin) ) ierr = ierr + 1
if ( (l2jmin > hwg_2jmax) .or. (l2jmax < hwg_2jmin) ) ierr = ierr + 1
if ( (rpmin > hwg_pmax) .or. (rpmax < hwg_pmin)  ) ierr = ierr + 1
if ( (lpmin > hwg_pmax) .or. (lpmax < hwg_pmin)  ) ierr = ierr + 1
if ( (hwg_phys == 0) .and. ((rZ /= hwg_Zv) .or. (rN /= hwg_Nv)) ) ierr = ierr+1

if ( ierr == 0 ) then
  print '(/,a,/,a)','Warning: the code takes into account the combined min/max &
         &values of','the input parameters and the weights.'
 
  hwg_2jmin = max(hwg_2jmin,l2jmin,r2jmin)
  hwg_2jmax = min(hwg_2jmax,l2jmax,r2jmax)
 
  hwg_pmin = max(hwg_pmin,lpmin,rpmin)
  hwg_pmax = min(hwg_pmax,lpmax,rpmax)
 
  print '(/,"HWG parameters",4x,"Value",/,23("-"))'
  print format2, 'Minimum 2*J      ', hwg_2jmin
  print format2, 'Maximum 2*J      ', hwg_2jmax
  print format2, 'Minimum P        ', hwg_pmin
  print format2, 'Minimum P        ', hwg_pmax
else
  print '(/,a,/,a)','Error: the min/max values of the input parameters are not &   
         &compatible','with the ones of the weights.'
  stop
endif

!!! Allocations
allocate( projme_label(projme_tdim,-3:0,hwg_pmin:hwg_pmax,2), &
          projme_2kj(projme_tdim,-3:0,hwg_pmin:hwg_pmax,2),   &
          stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of arrays for projected &
                        &matrix elements (weights)'

projme_label = 0
projme_2kj = 0

end subroutine read_weights_initial

!------------------------------------------------------------------------------!
! subroutine read_weights_blockJP                                              !
!                                                                              !
! Reads the weights from the files "X_weights.txt" (X = left or right) and     !
! store them.                                                                  !
!                                                                              !
! Input: block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "   "      "       "                         !
!                                                                              !
! Remark: right now, this is only used in the calculation of the nuclear matrix!
! elements of the neutrinoless double beta decay (not included in this version)!
!------------------------------------------------------------------------------!
subroutine read_weights_blockJP(block_2j,block_p)  

integer, intent(in) :: block_2j, block_p
integer :: i, j, k, n2j, np, n2kj, ntdim, n2jmin, n2jmax, npmin, npmax, nZ, nN,&
           m2j, mp, bdim, ldim, utkw
integer(i64) :: label 
character(len=*), parameter :: format1 = '(7i5)', &
                               format2 = '(2i3,1i5)', &
                               format3 = '(1i19,1i5,*(d24.16))'

!!! Initialization
projme_bdim(0,block_p,:) = 0
projme_ldim(0,block_p,:) = 0
projme_label(:,0,block_p,:) = 0
projme_2kj(:,0,block_p,:) = 0

!!! Reading the right and left weights
do j = 1, 2
  
  if ( j == 1 ) then
    utkw = utrw
  else
    utkw = utlw
  endif 

  read(utkw,format1) ntdim, n2jmin, n2jmax, npmin, npmax, nZ, nN

  do n2j = n2jmin, n2jmax, 2
    do np = npmax, npmin, -2
 
      read(utkw,format2) m2j, mp, bdim, ldim
      
      if ( (m2j /= block_2j) .or. (mp /= block_p) ) then 
        do i = 1, bdim
          read(utkw,format3) 
        enddo
      else
        projme_bdim(0,block_p,j) = bdim
        projme_ldim(0,block_p,j) = ldim
        
        do i = 1, bdim
          read(utkw,format3) label, n2kj, (weights_f(i,k,0,block_p,j),k=1,bdim)
          projme_label(i,0,block_p,j) = label
          projme_2kj(i,0,block_p,j) = n2kj   
        enddo
      endif
  
    enddo
  enddo
  
  rewind(utkw)
 
  if ( hwg_phys == 0 ) exit

enddo 

end subroutine read_weights_blockJP

!------------------------------------------------------------------------------!
! subroutine solve_hwg_sqrt                                                    !
!                                                                              !
! Solves the HWG equation through the steps:                                   !
! 1) Diagonalizes the norm matrix                                              !
! 2) Discards the norm eigenvalues below a cutoff and forms the "natural basis"!
!    from the normalized norm eigenvalues.                                     !
! 3) Diagonalizes the Hamiltonian in this basis.                               !
!                                                                              !
! Additionally, the collective wave functions are computed.                    !
!                                                                              !
! Input: ndim = dimension of the matrices (number of states)                   !
!        block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "   "      "       "                         !
!------------------------------------------------------------------------------!
subroutine solve_hwg_sqrt(block_2j,block_p,ndim)  

integer, intent(in) :: ndim, block_2j, block_p
integer :: i, k, m, n, imin, zdim, cdim, info 
real(r64) :: xener, xpari, xprot, xneut, xnucl, xamj2, xist2, xj, xt,  g2sum, &
             cond2
real(r64), dimension(ndim) :: eigen_norm, eigen_ener
real(r64), dimension(3*ndim-1) :: work1   
real(r64), dimension(ndim,ndim) :: A, C, D, Dn, F, L, G
character(1) :: chP
character(5) :: chJ
character(50) :: namefile
character(len=*), parameter :: format1 = '(10es14.6)', &
                               format2 = '(1i5,1x,1es14.6,1x,1f13.6,4f12.6, &
                                          &2f10.5)', &
                               format3 = '(1i5,1x,1f13.6,4f12.6,2f10.5)', &
                               format4 = '(1i5,2x,1i5,2x,1i19,2f12.6)'

!!!
!!! Opens the data files to write the information
!!!

chP = adjustr(chfile_P(block_p))
chJ = adjustr(chfile_J(block_2j))

namefile = "convergence_norm_eigenvalues_" // trim(adjustl(chJ)) // &
           trim(chP) // ".txt"
open(utd1, file=namefile, status='replace', action='write', form='formatted')

write(utd1,'(4x,"i",6x,"Energy",9x,"P",11x,"Z",11x,"N",11x,"A",10x,"J", &
           &9x,"T",/,87("-"))')

namefile = "collective_wavefunction_" // trim(adjustl(chJ)) // &
           trim(chP) // ".txt"
open(utd2, file=namefile, status='replace', action='write', form='formatted')

write(utd2,'(4x,"i",2x,"index",9x,"label",12x,"g(cwf)",6x,"g**2", &
           &/,58("-"))')

!!!
!!! Diagonalizes the norm
!!!

D(:,:) = projme_over(1:ndim,1:ndim)

call dsyev('V','U',ndim,D,ndim,eigen_norm,work1,3*ndim-1,info)

print '(/,"Norm eigenvalues",/,16("="),/)'

do i = 1, 1+ndim/8
  n = 8 * (i-1)
  m = min(8,ndim,ndim-n)
  if ( m /= 0 ) write(uto,format1) (eigen_norm(ndim+1-(k+n)),k=1,m)
enddo

cond2 = maxval(abs(eigen_norm)) / minval(abs(eigen_norm))

print '(/,"Condition number (2-norm) of the matrix :",1es9.2)', cond2   

!!! Determines the eigenvalues to keep/discard according to the cutoff
cdim = 0
zdim = 0

do i = 1, ndim
  if ( eigen_norm(i) >= eigen_norm(ndim) * &
                        cutoff_block_algo(block_2j,block_p) ) exit
  if ( eigen_norm(i) < -cutoff_block_negev(block_2j,block_p) ) zdim = i
  cdim = i
enddo

states_cdim(block_2j,block_p) = cdim
warnings_norm(block_2j,block_p) = zdim

print '(/,"Number of norm eigenvalues < zero_eps   :",1i5)', zdim
print '(  "Number of norm eigenvalues before cutoff:",1i5)', ndim
print '(  "Number of norm eigenvalues after  cutoff:",1i5)', ndim-cdim

if ( cdim == projme_bdim(0,block_p,1) ) then 
  print '(/,"No states left after cutoff on the norm eigenvalues.")'  
  return
endif

!!! Builds the collective subspace with proper normalization             
do i = 1, ndim
  if ( i <= cdim ) then
    Dn(:,i) = zero
  else
    Dn(:,i) = D(:,i) / sqrt(eigen_norm(i))
  endif 
enddo

call calculate_properties_spectrum(Dn(1:ndim,1:ndim),block_2j,block_p,ndim)

print '(/,"Natural states",/,14("="),//, &
       &4x,"i",5x,"Norm ev",9x,"Energy", &
       &9x,"P",11x,"Z",11x,"N",11x,"A",10x,"J",9x,"T",/,102("-"))'

do i = ndim, cdim+1, -1
  xener = states_ener(i,block_2j,block_p)
  xpari = states_pari(i,block_2j,block_p)
  xprot = states_prot(i,block_2j,block_p)
  xneut = states_neut(i,block_2j,block_p)
  xnucl = states_nucl(i,block_2j,block_p)
  xamj2 = states_amj2(i,block_2j,block_p)
  xist2 = states_ist2(i,block_2j,block_p)

  xj = 0.5d0 * (-1 + sqrt(1+4*abs(xamj2)))
  xt = 0.5d0 * (-1 + sqrt(1+4*abs(xist2)))

  write(uto,format2) ndim+1-i, eigen_norm(i), xener, xpari, xprot, xneut, & 
                     xnucl, xj, xt
enddo 

!!!
!!! Diagonalizes the Hamiltonian and computes the weigths f
!!!

print '(/,"Energy eigenvectors",/,19("="),//, &
       &4x,"i",6x,"Energy", &
       &9x,"P",11x,"Z",11x,"N",11x,"A",10x,"J",9x,"T",/,87("-"))'

!!! If the convergence analysis option is on, the diagonalization is performed
!!! several times, increasing the number of norm eigenvalues included in the
!!! calculation.
if ( hwg_conv == 0 ) then
  imin = ndim-cdim
else
  imin = 1
endif

do i = imin, ndim-cdim
  L = zero
  L(:,ndim+1-i:ndim) = Dn(:,ndim+1-i:ndim)

  C(:,:) = projme_ener(1:ndim,1:ndim)

  call dgemm('t','n',ndim,ndim,ndim,one,L,ndim,C,ndim,zero,A,ndim)
  call dgemm('n','n',ndim,ndim,ndim,one,A,ndim,L,ndim,zero,C,ndim)

  !!! Sets the energy of the zero blocks to a very high value (to be ignored)
  do k = 1, ndim
    if ( sum(abs(C(k,:))) < zero_eps ) C(k,k) = 100000.d0
  enddo

  !!! Diagonalization
  call dsyev('V','U',ndim,C,ndim,eigen_ener,work1,3*ndim-1,info)

  !!! Weights F
  call dgemm('n','n',ndim,ndim,ndim,one,L,ndim,C,ndim,zero,F,ndim)

  !!! Properties of the states
  call calculate_properties_spectrum(F(1:ndim,1:ndim),block_2j,block_p,ndim)
  states_ener(1:ndim,block_2j,block_p) = eigen_ener(:)

  !!! Tag the states that are non-vanishing     
  if ( i == ndim-cdim ) then 
    states_notzero(:,block_2j,block_p) = 0 
    
    do k = 1, ndim
      if (     states_ener(k,block_2j,block_p)  < 99999.0d0 .and. &
           abs(states_pari(k,block_2j,block_p)) > zero_eps  .and. &
           abs(states_nucl(k,block_2j,block_p)) > zero_eps ) then   
        states_notzero(k,block_2j,block_p) = 1
      endif 
    enddo
  endif

  !!! Prints/writes the properties
  do k = 1, min(i,ndim-cdim) 
    xener = states_ener(k,block_2j,block_p)
    xpari = states_pari(k,block_2j,block_p)
    xprot = states_prot(k,block_2j,block_p)
    xneut = states_neut(k,block_2j,block_p)
    xnucl = states_nucl(k,block_2j,block_p)
    xamj2 = states_amj2(k,block_2j,block_p)
    xist2 = states_ist2(k,block_2j,block_p)
 
    xj = 0.5d0 * (-1 + sqrt(1+4*abs(xamj2)))
    xt = 0.5d0 * (-1 + sqrt(1+4*abs(xist2)))
 
    if ( i == ndim-cdim ) then 
      write(uto,format3) k, xener, xpari, xprot, xneut, xnucl, xj, xt
    endif
    write(utd1,format3) k, xener, xpari, xprot, xneut, xnucl, xj, xt
  enddo
  write(utd1,'(" ")') 

  if ( i == ndim-cdim ) then 
    weights_f(:,:,0,block_p,1) = zero 
    weights_f(1:ndim,1:ndim,0,block_p,1) = F(1:ndim,1:ndim)
  endif
enddo

!!!
!!! Computes the collective wave functions
!!!

call dgemm('n','n',ndim,ndim,ndim,one,D,ndim,C,ndim,zero,G,ndim)

do i = 1, ndim-cdim
  g2sum = zero
  do k = 1, projme_bdim(0,block_p,1)
    write(utd2,format4) i, k, projme_label(k,0,block_p,1), G(k,i), G(k,i)**2
    g2sum = g2sum +  G(k,i)**2
  enddo
  write(utd2,'(49x,8("-"),/,45x,1f12.6,/)') g2sum
enddo

close(utd1, status='keep')
close(utd2, status='keep')

end subroutine solve_hwg_sqrt     

!------------------------------------------------------------------------------!
! subroutine solve_hwg_qz                                                      !
!                                                                              !
! Solves the HWG equation using the QZ algorithm implemented in LAPACK. This   !
! option is experimental and should be avoided if possible.                    !
!                                                                              !
! Additionally, the collective wave functions are computed.                    !
!                                                                              !
! Input: ndim = dimension of the matrices (number of states)                   !
!        block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "   "      "       "                         !
!------------------------------------------------------------------------------!
subroutine solve_hwg_qz(block_2j,block_p,ndim)  

integer, intent(in) :: ndim, block_2j, block_p
integer :: i, j, k, m, n, zdim, cdim, sdim, info 
integer, dimension(1) :: ntab
integer, dimension(ndim) :: eord
real(r64) :: xener, xpari, xprot, xneut, xnucl, xamj2, xist2, xj, xt,  g2sum, &
             cond2, maxbeta 
real(r64), dimension(ndim) :: eigen_norm, eigen_ener, alphar, alphai, beta
real(r64), dimension(3*ndim-1) :: work1   
real(r64), dimension(8*ndim) :: work2   
real(r64), dimension(ndim,ndim) :: A, C, D, F, G, VL, VR
character(1) :: chP
character(5) :: chJ
character(50) :: namefile
character(len=*), parameter :: format1 = '(10es14.6)', &
                               format2 = '(1i5,1x,1es14.6,1x,1f13.6,4f12.6, &
                                          &2f10.5)', &
                               format3 = '(1i5,1x,1f13.6,4f12.6,2f10.5)', &
                               format4 = '(1i5,2x,1i5,2x,1i19,2f12.6)', &
                               format5 = '(1i5,2x,4es14.5)'

!!!
!!! Opens the data files to write the information
!!!

chP = adjustr(chfile_P(block_p))
chJ = adjustr(chfile_J(block_2j))

namefile = "collective_wavefunction_" // trim(adjustl(chJ)) // &
           trim(chP) // ".txt"
open(utd2, file=namefile, status='replace', action='write', form='formatted')

write(utd2,'(4x,"i",2x,"index",9x,"label",12x,"g(cwf)",6x,"g**2", &
           &/,58("-"))')

!!!
!!! Diagonalizes the norm
!!!

D(:,:) = projme_over(1:ndim,1:ndim)

call dsyev('V','U',ndim,D,ndim,eigen_norm,work1,3*ndim-1,info)

print '(/,"Norm eigenvalues",/,16("="),/)'

do i = 1, 1+ndim/8
  n = 8 * (i-1)
  m = min(8,ndim,ndim-n)
  if ( m /= 0 ) write(uto,format1) (eigen_norm(ndim+1-(k+n)),k=1,m)
enddo

cond2 = maxval(abs(eigen_norm)) / minval(abs(eigen_norm))

print '(/,"Condition number (2-norm) of the matrix :",1es9.2)', cond2   

!!! Builds the square root of the norm discarding negative eigenvalues 
sdim = 0

do i = 1, ndim
  if ( eigen_norm(i) < zero_eps ) sdim = i
enddo

warnings_norm(block_2j,block_p) = sdim

print '(/,"Number of norm eigenvalues < zero_eps  :",1i5)', sdim

A = zero
do i = sdim+1, ndim
  A(i,i) = sqrt(eigen_norm(i))
enddo

call dgemm('n','n',ndim,ndim,ndim,one,D,ndim,A,ndim,zero,C,ndim)
call dgemm('n','t',ndim,ndim,ndim,one,C,ndim,D,ndim,zero,A,ndim)

!!!
!!! Diagonalizes the Hamiltonian and computes the weigths f
!!!

!!! QZ algorithm
D(:,:) = projme_over(1:ndim,1:ndim)
C(:,:) = projme_ener(1:ndim,1:ndim)

call dggev('N','V',ndim,C,ndim,D,ndim,alphar,alphai,beta,VL,ndim,VR,ndim, &
           work2,8*ndim,info)

print '(/,"Generalzed eigenvalues",/,22("="),//, &
       &"Solutions of beta H f = alpha N f",//, &
       &4x,"i",5x,"real(alpha)",3x,"imag(alpha)",6x,"beta",7x, &
       &"r(alp)/beta",/,63("-"))'

do i = 1, ndim
  write(uto,format5) i, alphar(i), alphai(i), beta(i), alphar(i) / beta(i)
enddo

!!! Orders and selects the eigenvalues
zdim = 0
cdim = 0

eigen_ener = zero
do i = 1, ndim
  if ( abs(alphai(i)) > 1.d-16 ) then 
    zdim = zdim + 1
    cdim = cdim + 1
  elseif ( abs(beta(i)) < cutoff_block_algo(block_2j,block_p) ) then
    cdim = cdim + 1
  else 
    eigen_ener(i) = alphar(i) / beta(i)
  endif
enddo

warnings_cplx(block_2j,block_p) = zdim

print '(/,"Number of complex eigenvalues      :",1i5)', zdim
print '(  "Number of eigenvalues before cutoff:",1i5)', ndim
print '(  "Number of eigenvalues after  cutoff:",1i5)', ndim-cdim

VL = VR
beta = eigen_ener
ntab = 0
eord = 0
maxbeta = 666.d0

do i = 1, ndim                                                                   
  ntab = minloc(beta,mask=beta<maxbeta)
  eord(i) = ntab(1)      
  beta(ntab(1)) = maxbeta + i
enddo

beta = eigen_ener

do i = 1, ndim
  eigen_ener(i) = beta(eord(i))
  VR(:,i) = VL(:,eord(i))
enddo

!!! Normalization of the eigenvectors
D(:,:) = projme_over(1:ndim,1:ndim)
call dgemm('n','n',ndim,ndim,ndim,one,D,ndim,VR,ndim,zero,C,ndim)
call dgemm('t','n',ndim,ndim,ndim,one,VR,ndim,C,ndim,zero,D,ndim)

F = zero
do i = 1, ndim
  do j = 1, ndim
    F(i,j) = VR(i,j) / sqrt(D(j,j))
  enddo
enddo

!!! Prints the properties of the eigenstates
print '(/,"Energy eigenvectors",/,19("="),//, &
       &4x,"i",6x,"Energy", &
       &9x,"P",11x,"Z",11x,"N",11x,"A",10x,"J",9x,"T",/,87("-"))'

call calculate_properties_spectrum(F(1:ndim,1:ndim),block_2j,block_p,ndim)
states_ener(1:ndim,block_2j,block_p) = eigen_ener(:)

do i = 1, ndim-cdim
  xener = states_ener(i,block_2j,block_p)
  xpari = states_pari(i,block_2j,block_p)
  xprot = states_prot(i,block_2j,block_p)
  xneut = states_neut(i,block_2j,block_p)
  xnucl = states_nucl(i,block_2j,block_p)
  xamj2 = states_amj2(i,block_2j,block_p)
  xist2 = states_ist2(i,block_2j,block_p)

  xj = 0.5d0 * (-1 + sqrt(1+4*abs(xamj2)))
  xt = 0.5d0 * (-1 + sqrt(1+4*abs(xist2)))

  write(uto,format3) i, xener, xpari, xprot, xneut, xnucl, xj, xt
enddo

weights_f(:,:,0,block_p,1) = zero 
weights_f(1:ndim,1:ndim,0,block_p,1) = F(1:ndim,1:ndim)

!!!
!!! Computes the collective wave functions
!!!

call dgemm('n','n',ndim,ndim,ndim,one,A,ndim,F,ndim,zero,G,ndim)

do i = 1, ndim-cdim
  g2sum = zero
  do j = 1, projme_bdim(0,block_p,1)
    write(utd2,format4) i, j, projme_label(j,0,block_p,1), G(j,i), G(j,i)**2
    g2sum = g2sum +  G(j,i)**2
  enddo
  write(utd2,'(49x,8("-"),/,45x,1f12.6,/)') g2sum
enddo

close(utd2, status='keep')

end subroutine solve_hwg_qz       

!------------------------------------------------------------------------------!
! subroutine calculate_properties_spectrum                                      !
!                                                                              !
! Computes the expectation values of the various operators for the eigenstates !
! of N or H.                                                                   !
!                                                                              !
! Input: ndim = dimension of the matrices (number of states)                   !
!        block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "   "      "       "                         !
!        F = weights solutions of HWG                                          !
!------------------------------------------------------------------------------!
subroutine calculate_properties_spectrum(F,block_2j,block_p,ndim)  

integer, intent(in) :: ndim, block_2j, block_p
real(r64), dimension(ndim,ndim), intent(in) :: F

!!! Energy
call multiply_weights_diag(F,projme_ener(1:ndim,1:ndim), &
                           states_ener(1:ndim,block_2j,block_p),ndim)

!!! Parity 
call multiply_weights_diag(F,projme_pari(1:ndim,1:ndim), &
                           states_pari(1:ndim,block_2j,block_p),ndim)

!!! Number of nucleons
call multiply_weights_diag(F,projme_prot(1:ndim,1:ndim), &
                           states_prot(1:ndim,block_2j,block_p),ndim)
call multiply_weights_diag(F,projme_neut(1:ndim,1:ndim), &
                           states_neut(1:ndim,block_2j,block_p),ndim)
call multiply_weights_diag(F,projme_nucl(1:ndim,1:ndim), &
                           states_nucl(1:ndim,block_2j,block_p),ndim)

!!! Angular momentum and isospin
call multiply_weights_diag(F,projme_amj2(1:ndim,1:ndim), &
                           states_amj2(1:ndim,block_2j,block_p),ndim)
call multiply_weights_diag(F,projme_ist2(1:ndim,1:ndim), &
                           states_ist2(1:ndim,block_2j,block_p),ndim)

!!! Spin-orbit correction for the charge radius
call multiply_weights_diag(F,projme_r2so(1:ndim,1:ndim), &
                           states_r2so(1:ndim,block_2j,block_p),ndim)

end subroutine calculate_properties_spectrum

!------------------------------------------------------------------------------!
! subroutine calculate_occnumb_spectrum                                        !
!                                                                              !
! Computes the occupation numbers of the final eigenstates of H.               !
!                                                                              !
! Input: ndim = dimension of the matrices (max. number of states)              !
!        block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "   "      "       "                         !
!------------------------------------------------------------------------------!
subroutine calculate_occnumb_spectrum(block_2j,block_p,ndim)  

integer, intent(in) :: ndim, block_2j, block_p
integer :: i, k
real(r64), dimension(ndim,ndim) :: F
character(1) :: chP
character(5) :: chJ
character(50) :: namefile
character(len=*), parameter :: format1 = "(5i6,1i9,2f13.6)", &
                               format2 = "(9x,1a3,27x,2f13.6)"

!!! Opens the data files to write the information
chP = adjustr(chfile_P(block_p))
chJ = adjustr(chfile_J(block_2j))

namefile = "occupation_numbers_" // trim(adjustl(chJ)) // &
           trim(chP) // ".txt"
open(utd2, file=namefile, status='replace', action='write', form='formatted')

write(utd2,'(5x,"i",5x,"#",5x,"n",5x,"l",3x,"2*j",4x,"label",5x,"protons",6x, &
           &"neutrons",/,65("-"))')

!!! Computes the occupation numbers using the weights of the current states
F(:,:) = weights_f(:,:,0,block_p,1)

do i = 1, 2
  do k = 1, HOsh_dim
    call multiply_weights_diag(F,projme_occn(1:ndim,1:ndim,k,i), &
                               states_occn(1:ndim,k,i),ndim)
  enddo
enddo

!!! Writes the results in a file
do i = 1, ndim
  if ( i > 1 ) write(utd2,'(/)')
  do k = 1, HOsh_dim
    write(utd2,format1) i, k, HOsh_n(k), HOsh_l(k), HOsh_2j(k), & 
                        HOsh_na(k), states_occn(i,k,1), states_occn(i,k,2)
  enddo
  write(utd2,format2) "sum", sum(states_occn(i,:,1)), sum(states_occn(i,:,2))
enddo

close(utd2, status='keep')                                                       

end subroutine calculate_occnumb_spectrum

!------------------------------------------------------------------------------!
! subroutine calculate_transitions_elm                                         !
!                                                                              !
! Computes transition matrix elements for the various operators between the    !
! eigenstates of H.                                                            !
!                                                                              !
! Input: ndim = dimension of the matrices (max. number of states)              !
!        block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "   "      "       "                         !
!------------------------------------------------------------------------------!
subroutine calculate_transitions_elm(block_2j,block_p,ndim)  

integer, intent(in) :: ndim, block_2j, block_p
integer :: k, dj, dp, dpmin, djmin, ji, jf, pi, pf, lt, pt, ltmin, ltmax
real(r64), dimension(ndim,ndim) :: Fi, Ff

!!! Weights of the current states
Fi(:,:) = weights_f(:,:,0,block_p,1)

!!! Loop over all possible transitions
djmin = max(-6,-(block_2j-hwg_2jmin))
dpmin = hwg_pmin * hwg_pmax

do dj = djmin, 0, 2
  do dp = dpmin, 1, 2

    ji = block_2j
    jf = block_2j + dj 
    pi = block_p
    pf = block_p * dp

    ltmin = min(abs(ji-jf),ji+jf)
    ltmax = max(abs(ji-jf),ji+jf)

    if ( projme_bdim(dj/2,pf,1) == 0 ) cycle 
   
    !!! Weights of the final states 
    Ff(:,:) = weights_f(:,:,dj/2,pf,1)

    do k = 1, 5
      select case (k)
        !!! E0
        case (0)
          lt = 0
          pt = +1  
          if ( dp /= pt ) cycle
          if ( (lt < ltmin) .or. (lt > ltmax) ) cycle
          call multiply_weights_full(Ff,projme_r2p(:,:,dj/2),Fi, &
                                     transi_r2p(:,:,ji,pi,dj/2),ndim)
          call multiply_weights_full(Ff,projme_r2n(:,:,dj/2),Fi, &
                                     transi_r2n(:,:,ji,pi,dj/2),ndim)
          call multiply_weights_full(Ff,projme_r2m(:,:,dj/2),Fi, &
                                     transi_r2m(:,:,ji,pi,dj/2),ndim)
        !!! E1
        case (1)
          if ( .not. do_E1 ) cycle
          lt = 2
          pt = -1  
          if ( dp /= pt ) cycle
          if ( (lt < ltmin) .or. (lt > ltmax) ) cycle
          call multiply_weights_full(Ff,projme_E1(:,:,dj/2),Fi, &
                                     transi_E1(:,:,ji,pi,dj/2),ndim)
        !!! E2
        case (2)
          if ( .not. do_E2 ) cycle
          lt = 4
          pt = +1  
          if ( dp /= pt ) cycle
          if ( (lt < ltmin) .or. (lt > ltmax) ) cycle
          call multiply_weights_full(Ff,projme_E2(:,:,dj/2),Fi, &
                                     transi_E2(:,:,ji,pi,dj/2),ndim)
        !!! E3
        case (3)
          if ( .not. do_E3 ) cycle
          lt = 6
          pt = -1  
          if ( dp /= pt ) cycle
          if ( (lt < ltmin) .or. (lt > ltmax) ) cycle
          call multiply_weights_full(Ff,projme_E3(:,:,dj/2),Fi, &
                                     transi_E3(:,:,ji,pi,dj/2),ndim)
        !!! M1
        case (4)
          if ( .not. do_M1 ) cycle
          lt = 2
          pt = +1  
          if ( dp /= pt ) cycle
          if ( (lt < ltmin) .or. (lt > ltmax) ) cycle
          call multiply_weights_full(Ff,projme_M1(:,:,dj/2),Fi, &
                                     transi_M1(:,:,ji,pi,dj/2),ndim)
        !!! M2
        case (5)
          if ( .not. do_M2 ) cycle
          lt = 4
          pt = -1  
          if ( dp /= pt ) cycle
          if ( (lt < ltmin) .or. (lt > ltmax) ) cycle
          call multiply_weights_full(Ff,projme_M2(:,:,dj/2),Fi, &
                                     transi_M2(:,:,ji,pi,dj/2),ndim)
      end select

    enddo
  enddo
enddo

end subroutine calculate_transitions_elm

!------------------------------------------------------------------------------!
! subroutine store_solution                                                    !
!                                                                              !
! Stores the quantities (e.g. the weights) that will be needed in the next     !
! quantum number blocks.                                                       !
!                                                                              !
! Remark: we only store the blocks that are close in 2J/P to be used in        !
! subsequent calculations of transitions. But we could probably store all of   !
! them.                                                                        !
!                                                                              !
! Input: block_2j = value of 2*j for this symmetry block                       !
!        block_p  =   "   "   p   "   "      "       "                         !
!------------------------------------------------------------------------------!
subroutine store_solution(block_2j,block_p)

integer, intent(in) :: block_2j, block_p
integer :: i, k
character(len=*), parameter :: format1 = '(7i5)', &
                               format2 = '(2i3,1i5)', &
                               format3 = '(1i19,1i5,*(d24.16))'

!!! Writes the weights to tape and stores them for the next block
if ( (block_2j == hwg_2jmin) .and. (block_p == hwg_pmin) ) then 
  open(utrw, file="weights_f.txt", status='replace', action='write', &
             form='formatted')

  write(utrw,format1) projme_tdim, hwg_2jmin, hwg_2jmax, hwg_pmin, hwg_pmax, &
                      hwg_Zv, hwg_Nv
endif

write(utrw,format2) block_2j, block_p, projme_bdim(0,block_p,1),  & 
                    projme_ldim(0,block_p,1)

do i = 1, projme_bdim(0,block_p,1)
  write(utrw,format3) projme_label(i,0,block_p,1), projme_2kj(i,0,block_p,1), &
                      (weights_f(i,k,0,block_p,1),k=1,projme_bdim(0,block_p,1))  
enddo

!!! Stores the information on the states
if ( block_p == hwg_pmin ) then 
  do i = -2, 0 
    projme_bdim(i-1,:,1) = projme_bdim(i,:,1)
    projme_ldim(i-1,:,1) = projme_ldim(i,:,1)
    projme_label(:,i-1,:,1) = projme_label(:,i,:,1)
    projme_2kj(:,i-1,:,1) = projme_2kj(:,i,:,1)
    weights_f(:,:,i-1,:,1) = weights_f(:,:,i,:,1)
    if ( hwg_norm == 1 ) then
      projme_norm(:,:,i-1) = projme_norm(:,:,i)
    endif
  enddo
endif

end subroutine store_solution     

!------------------------------------------------------------------------------!
! subroutine print_warnings                                                    !
!                                                                              !
! Prints some warnings concerning the solving of the HWG equations. For now,   !
! only the number of significative negative norm eigenvalues.                  !
!------------------------------------------------------------------------------!
subroutine print_warnings

integer :: n2j, np
character(1) :: chP
character(5) :: chJ
character(len=*), parameter :: format1 = '(1a5,2x,1a1,1x,1i4)'

print '(/,60("%"),/,27x,"WARNINGS",/,60("%"),//, &
        &"Matrix elements missing",/,23("="),//, &
        &2x,1a3,2x,"P",2x,"Number",/,16("-"))',char_Jleg

do n2j = hwg_2jmin, hwg_2jmax, 2
  do np = hwg_pmax, hwg_pmin, -2
    chP = adjustr(char_P(np))
    chJ = adjustr(char_J(n2j))
    write(uto,format1) chJ, chP, warnings_zero(n2j,np)
  enddo
enddo

print '(/,"Significative norm ev < 0",/,25("="),//, &
        &2x,1a3,2x,"P",2x,"Number",/,16("-"))',char_Jleg

do n2j = hwg_2jmin, hwg_2jmax, 2
  do np = hwg_pmax, hwg_pmin, -2
    chP = adjustr(char_P(np))
    chJ = adjustr(char_J(n2j))
    write(uto,format1) chJ, chP, warnings_norm(n2j,np)
  enddo
enddo

if ( hwg_algo == 1 ) then
  print '(/,"Complex general eigenvalues",/,27("="),//, &
          &2x,1a3,2x,"P",2x,"Number",/,16("-"))',char_Jleg
 
  do n2j = hwg_2jmin, hwg_2jmax, 2
    do np = hwg_pmax, hwg_pmin, -2
      chP = adjustr(char_P(np))
      chJ = adjustr(char_J(n2j))
      write(uto,format1) chJ, chP, warnings_cplx(n2j,np)
    enddo
  enddo
endif

end subroutine print_warnings     

!------------------------------------------------------------------------------!
! subroutine print_complementary_files                                         !
!                                                                              !
! Prints the name of the complementary files that were generated during the    !
! execution of the code.                                                       !
!------------------------------------------------------------------------------!
subroutine print_complementary_files

integer :: k, n2j, np
character(1) :: chP
character(5) :: chJ
character(50) :: namefile
logical :: is_exist            

print '(/,60("%"),/,20x,"COMPLEMENTARY FILES",21x,/,60("%"),//, &
       &7x,"Description",22x,"File",/,62("-"))'

!!! Weights
inquire(file="weights_f.txt", exist=is_exist)
if ( is_exist ) print '(a)',"Weights f (solution HWG) : weights_f.txt"

!!! Spectrum
inquire(file="spectrum_sorted_energy.txt", exist=is_exist)
if ( is_exist ) print '(a)', "Spectrum (sorted E)      : & 
                             &spectrum_sorted_energy.txt"

inquire(file="spectrum_sorted_spinparity.txt", exist=is_exist)
if ( is_exist ) print '(a)',"Spectrum (sorted J^pi)   : & 
                       &spectrum_sorted_spinparity.txt"

!!! Transitions
inquire(file="transitions_sorted_energy.txt", exist=is_exist)
if ( is_exist ) print '(a)',"Transitions (sorted E)   : &
                       &transitions_sorted_energy.txt"

inquire(file="transitions_sorted_spinparity.txt", exist=is_exist)
if ( is_exist ) print '(a)',"Transitions (sorted J^pi): & 
                       &transitions_sorted_spinparity.txt"

!!! Convergence as a function of the number of norm eigenvalues
if ( hwg_conv == 0 ) then
  k = -1 
else
  k = hwg_2jmax
endif

do n2j = hwg_2jmin, k, 2
  do np = hwg_pmax, hwg_pmin, -2
    chP = adjustr(chfile_P(np))
    chJ = adjustr(chfile_J(n2j))

    namefile = "convergence_norm_eigenvalues_" // trim(adjustl(chJ)) // &
               trim(chP) // ".txt"
    inquire(file=namefile, exist=is_exist)
    if ( is_exist ) print '(a,a)',"Convergence (# norm ev)  : ", &
                                   adjustl(namefile)
  enddo
enddo

!!! Collective wave functions
do n2j = hwg_2jmin, hwg_2jmax, 2
  do np = hwg_pmax, hwg_pmin, -2
    chP = adjustr(chfile_P(np))
    chJ = adjustr(chfile_J(n2j))

    namefile = "collective_wavefunction_" // trim(adjustl(chJ)) // &
               trim(chP) // ".txt"
    inquire(file=namefile, exist=is_exist)
    if ( is_exist ) print '(a,a)',"Collective wave function : ", & 
                                  adjustl(namefile)
  enddo
enddo

!!! Occupation numbers
do n2j = hwg_2jmin, hwg_2jmax, 2
  do np = hwg_pmax, hwg_pmin, -2
    chP = adjustr(chfile_P(np))
    chJ = adjustr(chfile_J(n2j))

    namefile = "occupation_numbers_" // trim(adjustl(chJ)) // &
               trim(chP) // ".txt"
    inquire(file=namefile, exist=is_exist)
    if ( is_exist ) print '(a,a)',"Occupation numbers       : ", & 
                                  adjustl(namefile)
  enddo
enddo

end subroutine print_complementary_files

!------------------------------------------------------------------------------!
! subroutine print_spectrum                                                    !
!                                                                              !
! Prints the energy spectrum up to the desired Edis excitation energy. The full!
! energy spectrum is also written into two files (sorted according to J^Pi and !
! energy).                                                                     !
!------------------------------------------------------------------------------!
subroutine print_spectrum

integer :: i, ni, nj, np, n1, n2, n3, ialloc=0
integer, dimension(3) :: ntab
real(r64) :: xener, xexci, xener_gs, xpari, xprot, xneut, xnucl, xamj2, xist2, &
             xj, xt, xqs, xmu, xr2p=0.d0, xr2n=0.d0, xr2m=0.d0, xr2so, & 
             xr2ch=0.d0, cg
real(r64), dimension(:,:,:), allocatable :: xener_tab
character(1) :: chP
character(5) :: chJ
character(len=*), parameter :: format1 = '(1a5,2x,1a1,1i4,1f11.3,3f9.3,4f9.4, &
                                           &4f12.6,2f10.5)'

print '(/,60("%"),/,23x,"ENERGY SPECTRUM",/,60("%"))'

!!!
!!! Opens the data files to write the information
!!!
open(utd1, file='spectrum_sorted_energy.txt', status='replace', &
           action='write', form='formatted')

open(utd2, file='spectrum_sorted_spinparity.txt', status='replace', &
           action='write', form='formatted')

!!!
!!! Sorts the spectrum by excitation energy 
!!!
states_odim = projme_tdim * ((hwg_2jmax-hwg_2jmin+2)/2) &
                          * ((hwg_pmax-hwg_pmin+2)/2)   &
              - sum(states_cdim)

allocate ( xener_tab(projme_tdim,hwg_2jmin:hwg_2jmax,hwg_pmin:hwg_pmax), &
           states_tabijp(states_odim,3), stat=ialloc )
if ( ialloc /= 0 ) stop 'Error during allocation of arrays for the spectrum'

where ( states_notzero == 1 ) 
  xener_tab = states_ener
elsewhere
  xener_tab = 100000.0d0
end where

ntab = 0
states_imax = 0
states_tabijp = 0
  
do i = 1, states_odim
  ntab = minloc(xener_tab)                               
  if ( ntab(3) == 0 ) cycle                                                   
  n1 = ntab(1)
  n2 = ntab(2) + hwg_2jmin - 1
  n3 = ntab(3) - 2*kdelta(hwg_pmin,-1)
  xener = xener_tab(n1,n2,n3) 
  if ( xener > 99999.0d0 ) exit 
  states_imax = i
  states_tabijp(i,1) = n1
  states_tabijp(i,2) = n2
  states_tabijp(i,3) = n3
  xener_tab(n1,n2,n3) = 100001.0d0
enddo
 
xener_gs = minval(states_ener)

!!!
!!! Spectrum sorted by energy
!!!
print '(/,"Displayed up to E_exc = ",1i3," MeV",/,31("="),//, &
       &2x,1a3,2x,"P",3x,"n",4x,"Energy",5x,"E_exc",5x,"Qs",7x,"mu", & 
       &7x,"r_p",6x,"r_n",6x,"r_m",5x,"r_ch",8x,"P",11x,"Z",11x,"N",11x, & 
       &"A",10x,"J",9x,"T",/,154("-"))', hwg_Edis, char_Jleg

write(utd1,'(2x,1a3,2x,"P",3x,"n",4x,"Energy",5x,"E_exc",5x,"Qs",7x,"mu", & 
         &7x,"r_p",6x,"r_n",6x,"r_m",5x,"r_ch",8x,"P",11x,"Z",11x,"N",11x, & 
         &"A",10x,"J",9x,"T",/,154("-"))') char_Jleg

do i = 1, states_imax
  ni = states_tabijp(i,1)
  nj = states_tabijp(i,2)
  np = states_tabijp(i,3)

  xener = states_ener(ni,nj,np)
  xexci = xener - xener_gs

  if ( xener > 99999.0 ) cycle
  if ( ni > 999 ) cycle 

  xpari = states_pari(ni,nj,np)
  xprot = states_prot(ni,nj,np)
  xneut = states_neut(ni,nj,np)
  xnucl = states_nucl(ni,nj,np)
  xamj2 = states_amj2(ni,nj,np)
  xist2 = states_ist2(ni,nj,np)
  xr2so = states_r2so(ni,nj,np)

  xj = 0.5d0 * (-1 + sqrt(1+4*abs(xamj2)))
  xt = 0.5d0 * (-1 + sqrt(1+4*abs(xist2)))

  !!! Spectroscopic moment
  xqs = zero
  if ( do_E2 ) then 
    call ClebschGordan(nj,4,nj,nj,0,nj,cg)                                 
    xqs = sqrt(16*pi/5) * cg * transi_E2(ni,ni,nj,np,0) / sqrt(nj + one)
  endif

  !!! Magnetic moment       
  xmu = zero
  if ( do_M1 ) then 
    call ClebschGordan(nj,2,nj,nj,0,nj,cg)                                 
    xmu = sqrt(4*pi/3) * cg * transi_M1(ni,ni,nj,np,0) / sqrt(nj + one)
  endif

  !!! Charge radii (correction taken from Cipollone.2015.PhysRevC.92.014306 and
  !!! Reinhard.2021.PhysRevC.103.054310)
  xr2p  = sqrt( transi_r2p(ni,ni,nj,np,0) ) 
  xr2n  = sqrt( transi_r2n(ni,ni,nj,np,0) ) 
  xr2m  = sqrt( transi_r2m(ni,ni,nj,np,0) )
  if ( xr2p > zero_eps ) then
    xr2ch = sqrt( transi_r2p(ni,ni,nj,np,0) &
            + radius_rp2 + (hwg_Nn * one / hwg_Zn) * radius_rn2 &
            + 0.75d0 * (hbarc / mass_ma)**2 &
            + (one / hwg_Zn) * ((hbarc / mass_ma)**2) * xr2so )
  endif


  chP = adjustr(char_P(np))
  chJ = adjustr(char_J(nj))

  !!! Printing
  if ( xexci <= one*hwg_Edis ) then 
    write(uto,format1) chJ, chP, ni, xener, xexci, xqs, xmu, xr2p, xr2n, &
                       xr2m, xr2ch, xpari, xprot, xneut, xnucl, xj, xt
  endif 
  write(utd1,format1) chJ, chP, ni, xener, xexci, xqs, xmu, xr2p, xr2n, & 
                      xr2m, xr2ch, xpari, xprot, xneut, xnucl, xj, xt
enddo

!!!
!!! Spectrum sorted by spin-parity
!!!
write(utd2,'(2x,1a3,2x,"P",3x,"n",4x,"Energy",5x,"E_exc",5x,"Qs",7x,"mu",  &
         &7x,"r_p",6x,"r_n",6x,"r_m",5x,"r_ch",8x,"P",11x,"Z",11x,"N",11x, &
         &"A",10x,"J",9x,"T",/,154("-"))') char_Jleg

do nj = hwg_2jmin, hwg_2jmax, 2 
  do np = hwg_pmax, hwg_pmin, -2
    do ni = 1, projme_tdim

      if ( states_notzero(ni,nj,np) == 0 ) cycle

      xener = states_ener(ni,nj,np)
      xexci = xener - xener_gs

      xpari = states_pari(ni,nj,np)
      xprot = states_prot(ni,nj,np)
      xneut = states_neut(ni,nj,np)
      xnucl = states_nucl(ni,nj,np)
      xamj2 = states_amj2(ni,nj,np)
      xist2 = states_ist2(ni,nj,np)
      xr2so = states_r2so(ni,nj,np)

      xj = 0.5d0 * (-1 + sqrt(1+4*abs(xamj2)))
      xt = 0.5d0 * (-1 + sqrt(1+4*abs(xist2)))

      !!! Spectroscopic moment
      xqs = zero
      if ( do_E2 ) then 
        call ClebschGordan(nj,4,nj,nj,0,nj,cg)                                 
        xqs = sqrt(16*pi/5) * cg * transi_E2(ni,ni,nj,np,0) / sqrt(nj + one)
      endif
  
      !!! Magnetic moment       
      xmu = zero
      if ( do_M1 ) then 
        call ClebschGordan(nj,2,nj,nj,0,nj,cg)                                 
        xmu = sqrt(4*pi/3) * cg * transi_M1(ni,ni,nj,np,0) / sqrt(nj + one)
      endif

      !!! Charge radii (correction taken from Cipollone.2015.PhysRevC.92.014306)
      xr2p  = sqrt( transi_r2p(ni,ni,nj,np,0) ) 
      xr2n  = sqrt( transi_r2n(ni,ni,nj,np,0) ) 
      xr2m  = sqrt( transi_r2m(ni,ni,nj,np,0) ) 
      if ( xr2p > zero_eps ) then
        xr2ch = sqrt( transi_r2p(ni,ni,nj,np,0) &
                + radius_rp2 + (hwg_Nn * one / hwg_Zn) * radius_rn2 &
                + 0.75d0 * (hbarc / mass_ma)**2 &
                + (one / hwg_Zn) * ((hbarc / mass_ma)**2) * xr2so )
      endif

      chP = adjustr(char_P(np))
      chJ = adjustr(char_J(nj))

      !!! Printing
      write(utd2,format1) chJ, chP, ni, xener, xexci, xqs, xmu, xr2p, xr2n, & 
                          xr2m, xr2ch, xpari, xprot, xneut, xnucl, xj, xt
    enddo
  enddo
enddo

close(utd1, status='keep')
close(utd2, status='keep')

deallocate( xener_tab )

end subroutine print_spectrum     

!------------------------------------------------------------------------------!
! subroutine print_transitions_elm                                             !
!                                                                              !
! Prints the electromagnetic transitions up to the desired Edis excitation en- !
! ergy. All the transitions are also written into two files (sorted according  !
! to J^Pi and energy).                                                         !
!------------------------------------------------------------------------------!
subroutine print_transitions_elm

integer :: i, f, k, ni, nf, ji, jf, dj, pi, pf, dp, lt, pt, ltmin, ltmax, & 
           n1, n2, j1, j2, dj12, p1, p2, mdj, jfmin, jfmax, tprem
real(r64) :: ei, ef, de, exci, excf, ener_gs, proba_fm, proba_wu, reduced_me, &
             factor_conv
character(1) :: chPi, chPf
character(2) :: namet
character(4) :: chJileg, chJfleg
character(5) :: chJi, chJf
character(len=*), parameter :: format1 = '(1a5,2x,1a1,1i4,1f8.3,1a5,2x,1a1,1i4,&
                                          &2f8.3,5x,1a2,2f16.8)', &
                               format2 = '(53x,1a2,2f16.8)'
 
print '(/,60("%"),/,16x,"ELECTROMAGNETIC TRANSITIONS",/,60("%"))'

!!!
!!! Opens the data files to write the information
!!!
open(utd1, file='transitions_sorted_energy.txt', status='replace', &
           action='write', form='formatted')

open(utd2, file='transitions_sorted_spinparity.txt', status='replace', &
           action='write', form='formatted')

!!!
!!! Transition sorted by energy
!!!

chJileg = char_Jleg // 'i'
chJfleg = char_Jleg // 'f'

print '(/,"Displayed up to E_exc = ",1i3," MeV",/,31("="),/, &
       &69x,"B(Tl:i->f)",/, &
       &1x,1a4,1x,"Pi",2x,"ni",3x,"Ei_ex",1x,1a4,1x,"Pf",2x,"nf",3x,"Ef_ex",&
       &3x,"Ei-Ef",4x,"Type",7x,"(e fm)",10x,"(wu)",/,87("-"))', hwg_Edis, &
       chJileg, chJfleg

write(utd1,'(69x,"B(Tl:i->f)",/, &
        &1x,1a4,1x,"Pi",2x,"ni",3x,"Ei_ex",1x,1a4,1x,"Pf",2x,"nf",3x,"Ef_ex",&
        &3x,"Ei-Ef",4x,"Type",7x,"(e fm)",10x,"(wu)",/,87("-"))') chJileg, & 
        chJfleg

ener_gs = minval(states_ener)

do i = 2, states_imax
  ni = states_tabijp(i,1)
  ji = states_tabijp(i,2)
  pi = states_tabijp(i,3)
  if ( states_notzero(ni,ji,pi) == 0 ) cycle

  ei = states_ener(ni,ji,pi)
  exci = ei - ener_gs

  do f = 1, i-1
    nf = states_tabijp(f,1)
    jf = states_tabijp(f,2)
    pf = states_tabijp(f,3)
    if ( states_notzero(nf,jf,pf) == 0 ) cycle

    ef = states_ener(nf,jf,pf)
    excf = ef - ener_gs

    dj = (jf - ji)/2
    dp = pi * pf
    de = ei - ef 
   
    if ( abs(dj) > 3 ) cycle
    if ( (abs(dj) > 2) .and. (dp == 1) ) cycle

    !!! Inverse the indices if necessary (order of loop on j,p in the code)
    if ( (dj < 0) .or. ((dj == 0) .and. (pi == +1) .and. (pf == -1)) ) then 
      n1 = ni 
      n2 = nf 
      j1 = ji
      j2 = jf
      p1 = pi
      p2 = pf
    else 
      dj = -dj
      n1 = nf 
      n2 = ni
      j1 = jf
      j2 = ji
      p1 = pf
      p2 = pi
    endif

    ltmin = min(abs(ji-jf),ji+jf)
    ltmax = max(abs(ji-jf),ji+jf)
    tprem = 0 

    do k = 1, 5

      reduced_me = zero
      proba_fm = zero
      proba_wu = zero

      !!! Select the type of transition
      select case (k)
        ! ATTENTION: never benchmarked. Need to check definition & factor!
        case (0)
          namet = 'E0'
          lt = 0
          pt = +1  
          mdj = max(-lt/2,dj)
          reduced_me = transi_r2p(n2,n1,j1,p1,mdj) * sqrt(j1 + one)
          factor_conv = (radius_r0 * hwg_An)**4

        case (1)
          if ( .not. do_E1 ) cycle
          namet = 'E1'
          lt = 2
          pt = -1  
          mdj = max(-lt/2,dj)
          reduced_me = transi_E1(n2,n1,j1,p1,mdj)
          factor_conv = 0.06446d0 * (hwg_An**(2.d0/3.d0))

        case (2)
          if ( .not. do_E2 ) cycle
          namet = 'E2'
          lt = 4
          pt = +1  
          mdj = max(-lt/2,dj)
          reduced_me = transi_E2(n2,n1,j1,p1,mdj)
          factor_conv = 0.05940d0 * (hwg_An**(4.d0/3.d0))

        case (3)
          if ( .not. do_E3 ) cycle
          namet = 'E3'
          lt = 6
          pt = -1  
          mdj = max(-lt/2,dj)
          reduced_me = transi_E3(n2,n1,j1,p1,mdj)
          factor_conv = 0.05940d0 * (hwg_An**(6.d0/3.d0))
 
        case (4)
          if ( .not. do_M1 ) cycle
          namet = 'M1'
          lt = 2
          pt = +1  
          mdj = max(-lt/2,dj)
          reduced_me = transi_M1(n2,n1,j1,p1,mdj)
          factor_conv = 1.790d0

        case (5)
          if ( .not. do_M2 ) cycle
          namet = 'M2'
          lt = 4
          pt = -1  
          mdj = max(-lt/2,dj)
          reduced_me = transi_M2(n2,n1,j1,p1,mdj)
          factor_conv = 1.650d0 * (hwg_An**(2.d0/3.d0))
      end select

      if ( dp /= pt ) cycle
      if ( (lt < ltmin) .or. (lt > ltmax) ) cycle

      proba_fm = reduced_me**2 / (j1 + 1)  
      proba_wu = proba_fm / factor_conv

      if ( proba_fm < 5.d-8 ) cycle 

      chPi = adjustr(char_P(pi))
      chJi = adjustr(char_J(ji))
      chPf = adjustr(char_P(pf))
      chJf = adjustr(char_J(jf))

      !!! Printing
      ! Remark: possibly adapt the printing to the value (B(M1) are small)
      tprem = tprem + 1
        
      if ( tprem == 1 ) then 
        if ( exci <= one*hwg_Edis ) then 
          write(uto,format1) chJi, chPi, ni, exci, chJf, chPf, nf, excf, de, & 
                             namet, proba_fm, proba_wu
        endif
        write(utd1,format1) chJi, chPi, ni, exci, chJf, chPf, nf, excf, de, & 
                            namet, proba_fm, proba_wu
      else
        if ( exci <= one*hwg_Edis ) then 
          write(uto,format2) namet, proba_fm, proba_wu 
        endif
        write(utd1,format2) namet, proba_fm, proba_wu 
      endif

    enddo
  enddo 
enddo

!!!
!!! Transition sorted by spin-parity
!!!
write(utd2,'(69x,"B(Tl:i->f)",/, &
        &3x,"Ji",1x,"Pi",2x,"ni",3x,"Ei_ex",3x,"Jf",1x,"Pf",2x,"nf",3x,"Ef_ex",&
        &3x,"Ei-Ef",4x,"Type",7x,"(e fm)",10x,"(wu)",/,87("-"))')

do ji = hwg_2jmin, hwg_2jmax, 2 
  do pi = hwg_pmax, hwg_pmin, -2

    jfmin = max(ji-6,hwg_2jmin)
    jfmax = min(ji+6,hwg_2jmax)

    do jf = jfmin, jfmax, 2 
      do pf = hwg_pmax, hwg_pmin, -2

        dj = (jf - ji)/2
        dp = pi * pf
        if ( (abs(dj) > 2) .and. (dp == 1) ) cycle

        do ni = 1, projme_tdim
          if ( states_notzero(ni,ji,pi) == 0 ) cycle

          ei = states_ener(ni,ji,pi)
          exci = ei - ener_gs

          do nf = 1, projme_tdim
            if ( states_notzero(nf,jf,pf) == 0 ) cycle

            ef = states_ener(nf,jf,pf) 
            excf = ef - ener_gs
            de = ei - ef
            if ( de <= zero ) cycle

            !!! Inverse the indices if necessary (order of loops on j,p)
            if ( (dj < 0) .or. ((dj == 0) .and. (pi == +1) .and. &
                                (pf == -1)) ) then 
              dj12 = +dj
              n1 = ni 
              n2 = nf 
              j1 = ji
              j2 = jf
              p1 = pi
              p2 = pf
            else 
              dj12 = -dj
              n1 = nf 
              n2 = ni
              j1 = jf
              j2 = ji
              p1 = pf
              p2 = pi
            endif
        
            ltmin = min(abs(ji-jf),ji+jf)
            ltmax = max(abs(ji-jf),ji+jf)
            tprem = 0 
        
            do k = 1, 5
        
              reduced_me = zero
              proba_fm = zero
              proba_wu = zero
        
              !!! Select the type of transition
              select case (k)
                ! Never benchmarked. Need to check definition & factor!
                case (0)
                  namet = 'E0'
                  lt = 0
                  pt = +1  
                  mdj = max(-lt/2,dj12)
                  reduced_me = transi_r2p(n2,n1,j1,p1,mdj) * sqrt(j1 + one)
                  factor_conv = (radius_r0 * hwg_An)**4
        
                case (1)
                  if ( .not. do_E1 ) cycle
                  namet = 'E1'
                  lt = 2
                  pt = -1  
                  mdj = max(-lt/2,dj12)
                  reduced_me = transi_E1(n2,n1,j1,p1,mdj)
                  factor_conv = 0.06446d0 * (hwg_An**(2.d0/3.d0))
        
                case (2)
                  if ( .not. do_E2 ) cycle
                  namet = 'E2'
                  lt = 4
                  pt = +1  
                  mdj = max(-lt/2,dj12)
                  reduced_me = transi_E2(n2,n1,j1,p1,mdj)
                  factor_conv = 0.05940d0 * (hwg_An**(4.d0/3.d0))
        
                case (3)
                  if ( .not. do_E3 ) cycle
                  namet = 'E3'
                  lt = 6
                  pt = -1  
                  mdj = max(-lt/2,dj12)
                  reduced_me = transi_E3(n2,n1,j1,p1,mdj)
                  factor_conv = 0.05940d0 * (hwg_An**(6.d0/3.d0))
         
                case (4)
                  if ( .not. do_M1 ) cycle
                  namet = 'M1'
                  lt = 2
                  pt = +1  
                  mdj = max(-lt/2,dj12)
                  reduced_me = transi_M1(n2,n1,j1,p1,mdj)
                  factor_conv = 1.790d0
        
                case (5)
                  if ( .not. do_M2 ) cycle
                  namet = 'M2'
                  lt = 4
                  pt = -1  
                  mdj = max(-lt/2,dj12)
                  reduced_me = transi_M2(n2,n1,j1,p1,mdj)
                  factor_conv = 1.650d0 * (hwg_An**(2.d0/3.d0))
              end select

              if ( dp /= pt ) cycle
              if ( (lt < ltmin) .or. (lt > ltmax) ) cycle
        
              proba_fm = reduced_me**2 / (j1 + 1)  
              proba_wu = proba_fm / factor_conv
        
              if ( proba_fm < 5.d-8 ) cycle 

              chPi = adjustr(char_P(pi))
              chJi = adjustr(char_J(ji))
              chPf = adjustr(char_P(pf))
              chJf = adjustr(char_J(jf))
        
              !!! Printing
              tprem = tprem + 1
        
              if ( tprem == 1 ) then 
                write(utd2,format1) chJi, chPi, ni, exci, chJf, chPf, nf, excf,& 
                                    de, namet, proba_fm, proba_wu
              else
                write(utd2,format2) namet, proba_fm, proba_wu 
              endif
            enddo

          enddo
        enddo
      enddo
    enddo
  enddo
enddo

close(utd1, status='keep')
close(utd2, status='keep')

end subroutine print_transitions_elm

END MODULE Spectroscopy      
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
