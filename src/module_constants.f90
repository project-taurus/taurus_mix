!==============================================================================!
! MODULE Constants                                                             !
!                                                                              !
! This module contains the variables related to the physical, mathematical and !
! numerical constants.                                                         !
!==============================================================================!
MODULE Constants

use Iso_fortran_env

implicit none
public

!!! Definition of input/output units (for portability)
integer, parameter :: uti = input_unit,  &
                      uto = output_unit, &
                      utlw = uti + uto +  1, & ! left_weights.bin
                      utrw = uti + uto +  2, & ! right_weights.bin
                      utst = uti + uto +  3, & ! projmatelm_states.bin
                      ute1 = uti + uto +  4, & ! projmatelm_E1.bin        
                      ute2 = uti + uto +  5, & ! projmatelm_E2.bin         
                      ute3 = uti + uto +  6, & ! projmatelm_E3.bin          
                      utm1 = uti + uto +  7, & ! projmatelm_M1.bin         
                      utm2 = uti + uto +  8, & ! projmatelm_M2.bin        
                      ut0n = uti + uto +  9, & ! 0nu2beta.me   
                      utd1 = uti + uto + 15, & ! generic data file
                      utd2 = uti + uto + 16    ! generic data file

!!! Definition of kind parameters (for portability)
integer, parameter :: i8  = int8,   & ! integer  8 bits
                      i16 = int16,  & !    "    16  "   
                      i32 = int32,  & !    "    32  "   
                      i64 = int64,  & !    "    64  "   
                      r32 = real32, & ! real 32 bits (single precision)
                      r64 = real64    !  "   64  "   (double     "    )

!!! Definition of simple names for numerical values
real(r64), parameter :: one  = 1.0d0, &
                        zero = 0.0d0, &
                        zero_eps = 1.0d-16
complex(r64), parameter :: zone  = (1.0d0,0.0d0), &
                           zzero = (0.0d0,0.0d0), &
                           zimag = (0.0d0,1.0d0)

!!! Definition of physical constants
real(r64), parameter :: pi = 4.0d0 * atan(1.0d0),  & ! 3.14159265...
                        hbarc = 197.3269788d0,     & ! hbar*c in Mev.fm
                        radius_r0 = 1.2d0,         & ! radius factor 
                        radius_rp2 = +0.8414d0**2, & ! ms charge radius proton 
                        radius_rn2 = -0.1161d0,    & ! ms charge radius neutron
                        mass_mp = 938.27208816d0,  & ! proton mass
                        mass_mn = 939.56542052d0,  & ! neutron mass
                        mass_ma = (mass_mp + mass_mn)/2, & ! nucleon mass
                        magmome_mup = +2.79284734d0, & ! magnetic moment proton
                        magmome_mun = -1.91304273d0, & !    "       "    neutron
                        gyro_glp = 1.0d0, & ! gyromagnetic factor L proton 
                        gyro_gln = 0.0d0, & !        "       "    " neutron 
                        gyro_gsp = +5.58569469d0, & !"       "    S proton  
                        gyro_gsn = -3.82608545d0    !"       "    " neutron 

END MODULE Constants
!==============================================================================!
! End of file                                                                  !
!==============================================================================!
