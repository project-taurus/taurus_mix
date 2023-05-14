!==============================================================================!
! MODULE Initialization                                                        !
!                                                                              !
! This module contains the variables and routines related to the reading of    !
! the input parameters and files.                                              !
!                                                                              !
! List of routines and functions:                                              !
! - subroutine print_version                                                   !
! - subroutine read_input                                                      !
! - subroutine print_input                                                     !
! - subroutine check_input                                                     !
! - subroutine open_files_projmatelem                                          !
!==============================================================================!
MODULE Initialization

use Cutoffs   

implicit none
private 

!!! Cosmetic input names
character(30), dimension(31) :: input_names ! general inputs
character(30), dimension(:), allocatable :: input_spec ! special cutoffs

!!! Public routines
public :: print_version, read_input

CONTAINS 

!------------------------------------------------------------------------------!
! subroutine print_version                                                     !
!                                                                              !
! Prints the panel with the version, the description and the logo of the code. !
! This is the first action of the code.                                        !
!------------------------------------------------------------------------------!
subroutine print_version

print '(" __________________________________________________________",/, &
      & "|                                                          |",/, &
      & "|  (______)                                                |",/, &
      & "|  <(0  0)>   TAURUS_mix, version 2023.05.14               |",/, &
      & "|    (°°)                                                  |",/, &
      & "|                                                          |",/, &
      & "| This code performs the configuration mixing of symmetry  |",/, &
      & "| -projected states by solving the Hill-Wheeler-Griffin    |",/, &
      & "| equation.                                                |",/, &
      & "|                                                          |",/, &
      & "| Licence: GNU General Public License version 3 or later   |",/, &
      & "| DOI: 10.5281/zenodo.XXXXXXX                              |",/, &
      & "| Git: https://github.com/project-taurus/taurus_mix.git    |",/, &
      & "|                                                          |",/, &
      & "| Contributors: B. Bally, T. Rodríguez                     |",/, &
      & "|__________________________________________________________|",/)' 

end subroutine print_version

!------------------------------------------------------------------------------!
! subroutine read_input                                                        !
!                                                                              !
! Reads the input parameters from the input file and perform somes checks to   !
! see if they have approriate values. Also, opens the files containing the     !
! matrix elements.                                                             !
!------------------------------------------------------------------------------!
subroutine read_input

integer :: i, ialloc=0
character(len=*), parameter :: format1 = '(1a)', &
                               format2 = '(1a30,1i1)', &
                               format3 = '(1a30,1i3)', &
                               format4 = '(1a30,1i5)', &
                               format5 = '(1a30,1f5.2)', &
                               format6 = '(1a30,1es10.3)', &
                               format7 = '(1a30,1a1,2i3,1x,1es10.3)', &
                               format8 = '(1a30,1a1,2i3,1i19,1i3)'

!!! Reads the input parameters
read(uti,format1) input_names(1)
read(uti,format1) input_names(2)
read(uti,format2) input_names(3),  hwg_phys
read(uti,format2) input_names(4),  hwg_algo
read(uti,format2) input_names(5),  hwg_norm
read(uti,format2) input_names(6),  hwg_rmev
read(uti,format2) input_names(7),  hwg_conv
read(uti,format3) input_names(8),  hwg_Edis
read(uti,format1) input_names(9) 
read(uti,format1) input_names(10)
read(uti,format1) input_names(11)
read(uti,format3) input_names(12), hwg_Z
read(uti,format3) input_names(13), hwg_N
read(uti,format3) input_names(14), hwg_Zc
read(uti,format3) input_names(15), hwg_Nc
read(uti,format3) input_names(16), hwg_2jmin
read(uti,format3) input_names(17), hwg_2jmax
read(uti,format3) input_names(18), hwg_pmin
read(uti,format3) input_names(19), hwg_pmax
read(uti,format5) input_names(20), hwg_echp
read(uti,format5) input_names(21), hwg_echn
read(uti,format1) input_names(22) 
read(uti,format1) input_names(23)
read(uti,format1) input_names(24)
read(uti,format4) input_names(25), cutoff_ldim  
read(uti,format6) input_names(26), cutoff_over
read(uti,format6) input_names(27), cutoff_algo
read(uti,format6) input_names(28), cutoff_negev
read(uti,format6) input_names(29), cutoff_J   
read(uti,format6) input_names(30), cutoff_A
read(uti,format4) input_names(31), cutoff_spec_dim
                                
allocate( input_spec(cutoff_spec_dim),        &
          cutoff_spec_2j(cutoff_spec_dim),    &
          cutoff_spec_p(cutoff_spec_dim),     &
          cutoff_spec_type(cutoff_spec_dim),  &
          cutoff_spec_value(cutoff_spec_dim), &
          cutoff_spec_label(cutoff_spec_dim), &
          cutoff_spec_lab2k(cutoff_spec_dim), &
          stat=ialloc ) 
if ( ialloc /= 0 ) stop 'Error during allocation of arrays for inputs'   

cutoff_spec_2j = 0
cutoff_spec_p = 0
cutoff_spec_type = 'X'
cutoff_spec_value = zero
cutoff_spec_label = 0
cutoff_spec_lab2k = 0
cutoff_label_dim = 0

do i = 1, cutoff_spec_dim
  read(uti,format7) input_spec(i), cutoff_spec_type(i)
  backspace uti
  if ( cutoff_spec_type(i) /= 'L' ) then
    read(uti,format7) input_spec(i), cutoff_spec_type(i), cutoff_spec_2j(i), &
                      cutoff_spec_p(i), cutoff_spec_value(i)
  else
    cutoff_label_dim = cutoff_label_dim + 1
    read(uti,format8) input_spec(i), cutoff_spec_type(i), cutoff_spec_2j(i), &
                      cutoff_spec_p(i), cutoff_spec_label(i), &
                      cutoff_spec_lab2k(i)
  endif
enddo

!!! Special values 
call set_characters_JP

hwg_A  = hwg_Z  + hwg_N
hwg_Ac = hwg_Zc + hwg_Nc

if ( hwg_Edis == 0 ) hwg_Edis=20

!!! Performs some tests on the value of the inputs and link the units for
!!! the files containing the matrix elements
call print_input
call check_input
call open_files_projmatelem   

end subroutine read_input

!------------------------------------------------------------------------------!
! subroutine print_input                                                       !
!                                                                              !
! Prints the input parameters at the beginning of the calculation in the same  !
! format such that it can be copied and reused in an input file.               !
!------------------------------------------------------------------------------!
subroutine print_input

integer :: i
character(3) :: hwg_Z_ch, hwg_N_ch, hwg_Zc_ch, hwg_Nc_ch, hwg_2jmin_ch, & 
                hwg_2jmax_ch, hwg_pmin_ch, hwg_pmax_ch, hwg_Edis_ch
character(5) :: hwg_echp_ch, hwg_echn_ch, cutoff_ldim_ch, cutoff_spec_dim_ch
character(10) :: cutoff_over_ch, cutoff_negev_ch, cutoff_algo_ch, cutoff_J_ch, &
                 cutoff_A_ch

character(len=*), parameter :: format1 = '(1a)', &
                               format2 = '(1a30,1i1)', &
                               format3 = '(1a30,1a3)', &
                               format4 = '(1a30,1a5)', &
                               format5 = '(1a30,1a10)', &
                               format6 = '(1a30,1a1,2i3,1x,1es10.3)', &
                               format7 = '(1a30,1a1,2i3,1i19,1i3)'

!!! Formats the variable to eliminate the unpleasant blank spaces
write(hwg_Z_ch,'(1i3)') hwg_Z  
write(hwg_N_ch,'(1i3)') hwg_N
write(hwg_Zc_ch,'(1i3)') hwg_Zc
write(hwg_Nc_ch,'(1i3)') hwg_Nc
write(hwg_2jmin_ch,'(1i3)') hwg_2jmin  
write(hwg_2jmax_ch,'(1i3)') hwg_2jmax  
write(hwg_pmin_ch,'(1i3)') hwg_pmin  
write(hwg_pmax_ch,'(1i3)') hwg_pmax  
write(hwg_echp_ch,'(1f5.2)') hwg_echp  
write(hwg_echn_ch,'(1f5.2)') hwg_echn  
write(hwg_Edis_ch,'(1i3)') hwg_Edis  
hwg_Z_ch = adjustl(hwg_Z_ch)
hwg_N_ch = adjustl(hwg_N_ch)
hwg_Zc_ch = adjustl(hwg_Zc_ch)
hwg_Nc_ch = adjustl(hwg_Nc_ch)
hwg_2jmin_ch = adjustl(hwg_2jmin_ch)
hwg_2jmax_ch = adjustl(hwg_2jmax_ch)
hwg_pmin_ch = adjustl(hwg_pmin_ch)
hwg_pmax_ch = adjustl(hwg_pmax_ch)
hwg_echp_ch = adjustl(hwg_echp_ch)
hwg_echn_ch = adjustl(hwg_echn_ch)
hwg_Edis_ch = adjustl(hwg_Edis_ch)

write(cutoff_ldim_ch,'(1i5)') cutoff_ldim
write(cutoff_spec_dim_ch,'(1i5)') cutoff_spec_dim
cutoff_ldim_ch = adjustl(cutoff_ldim_ch)
cutoff_spec_dim_ch = adjustl(cutoff_spec_dim_ch)

write(cutoff_over_ch,'(1es10.3)') cutoff_over
write(cutoff_negev_ch,'(1es10.3)') cutoff_negev
write(cutoff_algo_ch,'(1es10.3)') cutoff_algo
write(cutoff_J_ch,'(1es10.3)') cutoff_J
write(cutoff_A_ch,'(1es10.3)') cutoff_A
cutoff_over_ch = adjustl(cutoff_over_ch)
cutoff_negev_ch= adjustl(cutoff_negev_ch)
cutoff_algo_ch = adjustl(cutoff_algo_ch)
cutoff_J_ch = adjustl(cutoff_J_ch)
cutoff_A_ch = adjustl(cutoff_A_ch)

!!! Prints the input parameters
print '(60("%"),/,22x,"INPUT PARAMETERS",22x,/,60("%"),/)'
write(uto,format1) input_names(1)
write(uto,format1) input_names(2)
write(uto,format2) input_names(3),  hwg_phys
write(uto,format2) input_names(4),  hwg_algo
write(uto,format2) input_names(5),  hwg_norm
write(uto,format2) input_names(6),  hwg_rmev
write(uto,format2) input_names(7),  hwg_conv
write(uto,format3) input_names(8),  hwg_Edis_ch
write(uto,format1) input_names(9)  
write(uto,format1) input_names(10)
write(uto,format1) input_names(11) 
write(uto,format3) input_names(12), hwg_Z_ch     
write(uto,format3) input_names(13), hwg_N_ch     
write(uto,format3) input_names(14), hwg_Zc_ch     
write(uto,format3) input_names(15), hwg_Nc_ch     
write(uto,format3) input_names(16), hwg_2jmin_ch
write(uto,format3) input_names(17), hwg_2jmax_ch
write(uto,format3) input_names(18), hwg_pmin_ch
write(uto,format3) input_names(19), hwg_pmax_ch
write(uto,format4) input_names(20), hwg_echp_ch
write(uto,format4) input_names(21), hwg_echn_ch
write(uto,format1) input_names(22) 
write(uto,format1) input_names(23)
write(uto,format1) input_names(24)
write(uto,format4) input_names(25), cutoff_ldim_ch
write(uto,format5) input_names(26), cutoff_over_ch
write(uto,format5) input_names(27), cutoff_algo_ch
write(uto,format5) input_names(28), cutoff_negev_ch
write(uto,format5) input_names(29), cutoff_J_ch
write(uto,format5) input_names(30), cutoff_A_ch
write(uto,format4) input_names(31), cutoff_spec_dim_ch
do i = 1, cutoff_spec_dim
  if ( cutoff_spec_type(i) /= 'L' ) then
    write(uto,format6) input_spec(i), cutoff_spec_type(i), cutoff_spec_2j(i), &
                       cutoff_spec_p(i), cutoff_spec_value(i)
  else
    write(uto,format7) input_spec(i), cutoff_spec_type(i), cutoff_spec_2j(i), & 
                       cutoff_spec_p(i), cutoff_spec_label(i), &
                       cutoff_spec_lab2k(i)
  endif
enddo

print*,' '
 
end subroutine print_input

!------------------------------------------------------------------------------!
! subroutine check_input                                                       !
!                                                                              !
! Checks the values of the input parameters to see if they are appropriate.    !
! If not, the code will stop.                                                  !
!------------------------------------------------------------------------------!
subroutine check_input

integer :: i, ierror

!!! Counter for the number of errors (should be 0 at the end)
ierror = 0

!!!
!!! General parameters
!!!

if ( (hwg_phys < 0) .or. (hwg_phys > 0) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The physics case studied (hwg_phys) = ", & 
         hwg_phys," should be 0."
endif 

if ( (hwg_conv < 0) .or. (hwg_conv > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to perform a convegence analysis & 
         &(hwg_conv) = ",hwg_conv," should be 0 or 1."
endif 

if ( (hwg_norm < 0) .or. (hwg_norm > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to perform a normalization of the matrices & 
         &(hwg_norm) = ",hwg_norm," should be 0 or 1."
endif 

if ( (hwg_rmev < 0) .or. (hwg_rmev > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The option to remove the projected states given negatve & 
         &eigenvalues (hwg_rmev) = ",hwg_rmev," should be 0 or 1."
endif 

if ( (hwg_algo < 0) .or. (hwg_algo > 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The type of algorithm to solve the HWG equation & 
         &(hwg_algo) = ",hwg_algo," should be 0 or 1."
endif 

if ( (hwg_rmev == 1) .and. (hwg_algo == 1) ) then
  ierror = ierror + 1
  print "(a,1i1,a)","The removal of projected states has not been implemented & 
         &yet for the QZ algorithm."
endif 

if ( hwg_Edis < 0 ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The value of the maximum energy displayed in the tables & 
         &(hwg_Edis) = ",hwg_Edis," should be positive."
endif 

!!!
!!! Quantum numbers    
!!!

if ( hwg_Z < 0 ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The number of active protons (hwg_Z) = ",hwg_Z, &
        " should be positive or null."
endif 

if ( hwg_N < 0 ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The number of active neutrons (hwg_N) = ",hwg_N, &
        " should be positive or null."
endif 

if ( hwg_Zc < 0 ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The number of core protons (hwg_Zc) = ",hwg_Zc, &
        " should be positive or null."
endif 

if ( hwg_Nc < 0 ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The number of core neutrons (hwg_Nc) = ",hwg_Nc, &
        " should be positive or null."
endif 

if ( hwg_2jmin < 0 ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The minimum value for the angular momentum (hwg_2jmin) &
        &= ",hwg_2jmin," should be positive or null."
endif 

if ( hwg_2jmax < 0 ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The maximum value for the angular momentum (hwg_2jmax) &
        &= ",hwg_2jmax," should be positive or null."
endif 

if ( hwg_2jmin > hwg_2jmax ) then
  ierror = ierror + 1
  print "(a,1i3,a,1i3,a)","The maximum value for the angular momentum &
        &(hwg_2jmin) = ",hwg_2jmin," should be lesser than or equal to the & 
        &minimum value for the angular momentum (hwg_2jmax) = ",hwg_2jmax,"."
endif 

if ( (-1)**hwg_2jmin /= (-1)**hwg_2jmax ) then
  ierror = ierror + 1
  print "(a,1i2,a,1i2,a)","The number parity (-1**hwg_2jmin) = ", &
        (-1)**hwg_2jmin," should be equal to the number parity (-1**hwg_2jmax) &
        &= ",(-1)**hwg_2jmax,"."
endif 

if ( (-1)**(hwg_N+hwg_Z) /= (-1)**hwg_2jmin ) then
  ierror = ierror + 1
  print "(a,1i2,a,1i2,a)","The number parity (-1**(hwg_N+hwg_Z)) = ", &
        (-1)**(hwg_N+hwg_Z)," should be equal to the number parity &
        &(-1**hwg_2jmin) = ",(-1)**hwg_2jmin,"."
endif 

if ( (hwg_pmin /= -1) .and. (hwg_pmin /= 1) ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The minimum value for the parity (hwg_pmin) = ", & 
         hwg_pmin," should be -1 or 1."
endif 

if ( (hwg_pmax /= -1) .and. (hwg_pmax /= 1) ) then
  ierror = ierror + 1
  print "(a,1i3,a)","The maximum value for the parity (hwg_pmax) = ", & 
         hwg_pmax," should be -1 or 1."
endif 

if ( hwg_pmax < hwg_pmin ) then
  ierror = ierror + 1
  print "(a,1i3,a,1i3,a)","The maximum value for the parity (hwg_pmax) &
        &= ",hwg_pmax," should be greater than or equal to the minimum value & 
        &for the parity (hwg_pmin) = ",hwg_pmin,"."
endif 

!!!
!!! Cutoffs       
!!!

if ( cutoff_ldim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The cutoff for the maximum number of states &               
       &(cutoff_ldim) = ",cutoff_ldim," should be positive or null."
endif 

if ( cutoff_algo < zero ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The cutoff for the algorithm eigenvalues. &
        &(cutoff_algo) = ",cutoff_algo," should be positive or null."
endif 

if ( cutoff_J < zero ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The cutoff for the expectation values of the angular &     
        &momentum (cutoff_J) = ",cutoff_J," should be positive or null."
endif 

if ( cutoff_A < zero ) then
  ierror = ierror + 1
  print "(a,1es10.3,a)","The cutoff for the expectation values of the particle &     
        &numbers (cutoff_A) = ",cutoff_A," should be positive or null."
endif 

if ( cutoff_spec_dim < 0 ) then
  ierror = ierror + 1
  print "(a,1i5,a)","The number of specific cutoffs (cutoff_spec_dim) = ", &               
        cutoff_spec_dim," should be positive or null."
endif 

!!! Specific cutoffs

do i = 1, cutoff_spec_dim

  if ( cutoff_spec_2j(i) < 0 ) then
    ierror = ierror + 1
    print "(a,1i32,a,1i5,a)","The angular momentum of the specific cutoff &
          &(cutoff_spec_2j(i)) = ",cutoff_spec_2j(i)," for i = ",i, & 
          " should be positive or null."
  endif 

  if ( (cutoff_spec_p(i) /= -1) .and. (cutoff_spec_p(i) /= +1) ) then
    ierror = ierror + 1
    print "(a,1i3,a,1i5,a)","The parity of the specific cutoff &
          &(cutoff_spec_p(i)) = ",cutoff_spec_p(i)," for i = ",i, & 
          " should be +1 or -1."
  endif 

  if ( (cutoff_spec_type(i) /= 'O') .and. (cutoff_spec_type(i) /= 'S') .and. &
       (cutoff_spec_type(i) /= 'N') .and. (cutoff_spec_type(i) /= 'J') .and. &
       (cutoff_spec_type(i) /= 'A') .and. (cutoff_spec_type(i) /= 'L') ) then
    ierror = ierror + 1
    print "(a,1a1,a,1i5,a)","The type of the specific cutoff &
          &(cutoff_spec_type(i)) = ",cutoff_spec_type(i)," for i = ",i, & 
          " should be O, S, N, J, A, or L."
  endif 

  if ( cutoff_spec_value(i) < zero ) then
    ierror = ierror + 1
    print "(a,1es10.3,a,1i5,a)","The value of the specific cutoff &
          &(cutoff_spec_value(i)) = ",cutoff_spec_value(i)," for i = ",i, & 
          " should be positive or null."
  endif 

enddo

!!!
!!! Stops the code if an error has been found in the input file
!!!

if ( ierror /= 0 ) then
  print "(a,1i2,a)", "The code has dectected ",ierror," problem(s) with the &
        &input parameters and will stop. Please check the manual."
  stop 
endif

end subroutine check_input

!------------------------------------------------------------------------------!
! subroutine open_files_projmatelem                                            !
!                                                                              !
! Opens the matrix elements filmes and links the units if they exist. If not,  !
! the code will continue but will not carry the calculations related to the    !
! missing files (e.g. no B(M1)).                                               !
!------------------------------------------------------------------------------!
subroutine open_files_projmatelem   
 
integer :: ierror, iwarn
logical :: is_exist            
 
!!! Counter for the number of errors (should be 0 at the end)
ierror = 0
iwarn  = 0

!!!
!!! Files with the matrix elements
!!!

if ( hwg_phys == 0 ) then 
 
  !!! States    
  inquire (file="projmatelem_states.bin", exist=is_exist)
 
  if ( is_exist ) then
    open(utst, file="projmatelem_states.bin", status='old', action='read', &
         form='unformatted')
  else
    ierror = ierror + 1
    print "(a)","The file with the projected states (projmatelem_states.bin) &
          &can not be found."
  endif 
 
  !!! B(E1)                                
  if ( hwg_pmin /= hwg_pmax ) then
    inquire (file="projmatelem_E1.bin", exist=is_exist)
    
    if ( is_exist ) then
      do_E1 = .true.
      open(ute1, file="projmatelem_E1.bin", status='old', action='read', &
           form='unformatted')
    else
      iwarn = iwarn + 1
      print "(a)","The file with the B(E1) probabilities (projmatelem_E1.bin) &
            &can not be found."
    endif 
  endif 
  
  !!! B(E2)                                
  inquire (file="projmatelem_E2.bin", exist=is_exist)

  if ( is_exist ) then
    do_E2 = .true.
    open(ute2, file="projmatelem_E2.bin", status='old', action='read', &
         form='unformatted')
  else
    iwarn = iwarn + 1
    print "(a)","The file with the B(E2) probabilities (projmatelem_E2.bin) & 
          &can not be found."
  endif 
  
  !!! B(E3)                                
  if ( hwg_pmin /= hwg_pmax ) then
    inquire (file="projmatelem_E3.bin", exist=is_exist)
    
    if ( is_exist ) then
      do_E3 = .true.
      open(ute3, file="projmatelem_E3.bin", status='old', action='read', & 
           form='unformatted')
    else
      iwarn = iwarn + 1
      print "(a)","The file with the B(E3) probabilities (projmatelem_E3.bin) &
            &can not be found."
    endif 
  endif 
  
  !!! B(M1)                                
  inquire (file="projmatelem_M1.bin", exist=is_exist)

  if ( is_exist ) then
    do_M1 = .true.
    open(utm1, file="projmatelem_M1.bin", status='old', action='read', &
         form='unformatted')
  else
    iwarn = iwarn + 1
    print "(a)","The file with the B(M1) probabilities (projmatelem_M1.bin) &
          &can not be found."
  endif 
  
  !!! B(M2)                                
  if ( hwg_pmin /= hwg_pmax ) then
    inquire (file="projmatelem_M2.bin", exist=is_exist)
    
    if ( is_exist ) then
      do_M2 = .true.
      open(utm2, file="projmatelem_M2.bin", status='old', action='read', &
           form='unformatted')
    else
      iwarn = iwarn + 1
      print "(a)","The file with the B(M2) probabilities (projmatelem_M2.bin) &
            &can not be found."
    endif 
  endif 
 
  !!! Any BTl
  if ( do_E1 .or. do_E2 .or. do_E3 .or. do_M1 .or. do_M2 ) then
    do_Tl = .true.
  endif
 
endif 

!!!
!!! Stops the code if an error has been found with the files
!!!

if ( ierror /= 0 ) then
  print "(a,1i1,a)", "The code has dectected ",ierror," problem(s) with the &
        &matrix elements files and will stop. Please check the files."
  stop 
endif

end subroutine open_files_projmatelem   

END MODULE Initialization
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
