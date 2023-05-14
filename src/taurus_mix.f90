!==============================================================================!
! PROGRAM TAURUS_mix                                                           !
!                                                                              !
! This program performs the configuration mixing of symmetry-projected states  !
! by solving the Hill-Wheeler-Griffin equation.                                !
!                                                                              !
! Licence: GNU General Public License version 3 or later                       !
! DOI: 10.5281/zenodo.XXXXXXX                                                  !
! Repository: github.com/project-taurus/taurus_mix                             !
!                                                                              !
! Article describing the code:                                                 !
!==============================================================================!
PROGRAM TAURUS_mix

use Spectroscopy       
use Initialization

implicit none

!!! Quantum number blocks
integer :: block_2j, block_p 

!!!
!!! INITIALIZATION
!!!

call print_version        
call read_input           

if ( hwg_phys == 0 ) then 
  call set_cutoffs
  call read_projmatelem_initial  
  call set_spectroscopy
endif 

!!!
!!! SOLVING THE HWG EQUATION
!!!

do block_2j = hwg_2jmin, hwg_2jmax, 2
  do block_p = hwg_pmax, hwg_pmin, -2
  
    if ( hwg_phys == 0 ) then
      print '(/,60("%"),/,24x,"J^pi = ",1a7,/,60("%"))', &
            char_JP(block_2j,block_p)
      call print_cutoffs_blockJP(block_2j,block_p)
      call read_projmatelem_blockJP(block_2j,block_p)
    endif
  
    if ( projme_bdim(0,block_p,1) /= 0 ) then 
      call read_projmatelem_states(block_2j,block_p)
      call print_projmatelem_states(block_p)
      if ( hwg_norm == 1 ) call normalize_projmatelem_states(block_p)
      if ( hwg_rmev == 1 ) call remove_negev_projmatelem_states(block_2j, &
                                                                block_p)
      if ( do_Tl ) call read_projmatelem_transitions_elm(block_2j,block_p)
  
      if ( hwg_algo == 0 ) then 
        call solve_hwg_sqrt(block_2j,block_p,projme_bdim(0,block_p,1))
      else 
        call solve_hwg_qz(block_2j,block_p,projme_bdim(0,block_p,1))
      endif
   
      call calculate_transitions_elm(block_2j,block_p,projme_tdim)
    else
      print '(/,"No projected states have been found or selected.")'
    endif

    call store_solution(block_2j,block_p)
  enddo
enddo

!!!
!!! PRINTING THE RESULTS
!!!

if ( hwg_phys == 0 ) then 
  call print_warnings
  call print_spectrum
  if ( do_Tl ) call print_transitions_elm
  call print_complementary_files
endif 

print '(/,"This is the end, my only friend, the end.")'  

END PROGRAM TAURUS_mix
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
