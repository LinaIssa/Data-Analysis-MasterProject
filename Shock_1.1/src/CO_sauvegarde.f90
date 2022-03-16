! Copyright (C) 2019 - Observatoire de Paris
!                      CNRS (France)
!                      Contacts : antoine.gusdorf@lra.ens.fr
!                                 sylvie.cabrit@obspm.fr
!                                 benjamin.godard@obspm.fr
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

MODULE MODULE_CO
  !*****************************************************************************
  !** The module 'MODULE_CO' contains variables and subroutines related to    **
  !** the CO molecule                                                         **
  !**     * levels                                                            **
  !**     * Aij, quadrupole lines                                             **
  !**     * escape probabilities, optical depths                              **
  !*****************************************************************************

  USE MODULE_PHYS_VAR, ONLY : NCO_lev, choice_ep

  IMPLICIT NONE
  INCLUDE "precision.f90"

  !-----------------------------------------
  ! useful variables
  !-----------------------------------------
  REAL(KIND=DP)            :: jj
  INTEGER(KIND=LONG)       :: j

  !-----------------------------------------
  ! parameters relative to CO
  !-----------------------------------------
  REAL(KIND=DP), PARAMETER :: Brot_CO   = 57635.96828_DP  ! rotational constant (MHz) NIST value
  REAL(KIND=DP), PARAMETER :: Moment_CO = 0.11011_DP      ! dipolar moment (debye ; 1 debye = 1.e-18 CGS) NIST value
  
  !-----------------------------------------
  ! energy levels
  ! variables defined in READ_CO_LEVELS
  !-----------------------------------------
  INTEGER(KIND=LONG)            :: Jmin_CO , Jmax_CO ! min. and max. values for J (rotation)

  !-----------------------------------------
  ! other useful quantities and constants
  !-----------------------------------------
  REAL(KIND=DP)                 :: gradv                              ! velocity gradient
  REAL(KIND=DP)                 :: Tb = 2.7_DP                        ! cosmic background temperature
!  REAL(KIND=DP)                 :: TAU_LIM = 1.D-20
!  REAL(KIND=DP)                 :: TAU_LIM = 0._DP
  REAL(KIND=DP)                 :: population_limit = 1.D-20
  REAL(KIND=DP)                 :: PF, bolj                           ! partition function
  REAL(KIND=DP)                 :: CST1, CST2, CST3, CST4, CST5, CST6 ! useful constants through the code
  INTEGER                       :: err                                ! allocation  parameter

  !-----------------------------------------
  !  radiative and collisional source terms
  !-----------------------------------------
  REAL(KIND=DP), DIMENSION(:),   ALLOCATABLE :: A_CO
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Radi, Coll
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Coll_oH2, Coll_pH2, Coll_H, Coll_He
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Matrix
  REAL(KIND=DP)                              :: Dens_paraH2, Dens_orthoH2

  !-----------------------
  ! data type of one level
  !-----------------------
  TYPE TYPE_CO_LEVEL
     INTEGER(KIND=LONG) :: J              ! rotational quantum number
     REAL(KIND=DP)      :: Weight         ! statistical weight 
     REAL(KIND=DP)      :: Energy         ! energy (K)
     REAL(KIND=DP)      :: DREL           ! calculation of Boltzmann's populations (no unit)
     REAL(KIND=DP)      :: fracbol        ! Boltzmann's fractional populations (no unit)
     REAL(KIND=DP)      :: fracpop        ! fractional population (no unit)
  END TYPE TYPE_CO_LEVEL

  ! initialized in READ_CO_LEVELS
  TYPE(TYPE_CO_LEVEL),DIMENSION(:), ALLOCATABLE  :: CO_lev ! vector containing the CO levels

  !--------------------------
  ! data type of one CO line
  !--------------------------
  TYPE TYPE_CO_LINE
     CHARACTER(LEN=9)   :: name              ! name of the line
     INTEGER(KIND=LONG) :: NCO_up, NCO_low   ! index of the upper and the lower levels
     REAL(KIND=DP)      :: A_CO              ! Aij (s-1)
     REAL(KIND=DP)      :: BOLB, BOLE        ! useful terms for radiative matrix calculation (no unit)
     REAL(KIND=DP)      :: DeltaE            ! energy of the line : Eup-Elow (K)
     REAL(KIND=DP)      :: Tdv, Tdv2         ! integrated intensity of the line (K.km/s)
     REAL(KIND=DP)      :: Tdv_old, Tdv2_old ! integrated intensity of the line (K.km/s) at the last call to DVODE
     REAL(KIND=DP)      :: intensity         ! intensity (Goldreich & Kwan style, K)
     REAL(KIND=DP)      :: intensity2        ! intensity (Surdej style, K)
     REAL(KIND=DP)      :: Texc              ! excitation temperature (K)
     REAL(KIND=DP)      :: TAU               ! optical depth (no unit)
     REAL(KIND=DP)      :: EP1               ! escape probability, statistical equilibrium (no unit)
     REAL(KIND=DP)      :: EP2               ! escape probability, calculation of intensity (no unit)
  END TYPE TYPE_CO_LINE

  !--- vector containing the CO lines ---
  TYPE (TYPE_CO_LINE), DIMENSION(:), SAVE, ALLOCATABLE :: CO_lines ! lines

  !---------------------------------------------
  ! collision rates
  ! read in READ_CO_RATES, used in EVOLUTION_CO
  !---------------------------------------------

  !--- CO-oH2 collision rates (Flower, 2001) ---
  ! max value of J for which collision rates are available
  INTEGER(KIND=LONG), PARAMETER                                        :: Jmin_oH2 = 0, Jmax_oH2 = 40
  ! number of available temperatures for tabulated collision rates
  INTEGER(KIND=LONG), PARAMETER                                        :: nt_oH2   = 14  
  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
  REAL(KIND=DP), DIMENSION(Jmin_oH2:Jmax_oH2,Jmin_oH2:Jmax_oH2,nt_oH2) :: Rate_oH2
  ! temperatures which defines Rate_CO(i,j,k)
  REAL(KIND=DP), DIMENSION(nt_oH2)                                     :: Temp_oH2 
  ! new rates, calculated at Tn (useful collision rate array)
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                           :: CR_oH2
  ! name of the file
  CHARACTER(LEN=75)                                                    :: filen1_oH2$ = 'input/coeff_12CO_oH2_leiden.in'
  ! format of temperatures in the file
  CHARACTER(LEN=*), PARAMETER                                          :: format_oH2_temp = '(8X,14F16.1)'  

  !--- CO-pH2 collision rates (Flower, 2001) ---
  ! max value of J for which collision rates are available
  INTEGER(KIND=LONG), PARAMETER                                        :: Jmin_pH2 = 0, Jmax_pH2 = 40 
  ! number of available temperatures for tabulated collision rates
  INTEGER(KIND=LONG), PARAMETER                                        :: nt_pH2   = 14  
  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
  REAL(KIND=DP), DIMENSION(Jmin_pH2:Jmax_pH2,Jmin_pH2:Jmax_pH2,nt_pH2) :: Rate_pH2
  ! temperatures which defines Rate_CO(i,j,k)
  REAL(KIND=DP), DIMENSION(nt_pH2)                                     :: temp_pH2 
  ! new rates, calculated at Tn (useful collision rate array)
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                           :: CR_pH2
  ! name of the file
  CHARACTER(LEN=75)                                                    :: filen1_pH2$ = 'input/coeff_12CO_pH2_leiden.in'
  ! format of temperatures in the file
  CHARACTER(LEN=*), PARAMETER                                          :: format_pH2_temp = '(8X,14F16.1)'  

  !--- CO-H collision rates (5 K < T < 3000 K ; Balakrishnan et al., ApJ, 2002) ---
  ! max value of J for which collision rates are available
  INTEGER(KIND=LONG), PARAMETER                              :: Jmin_H = 0, Jmax_H = 16
  ! number of available temperatures for tabulated collision rates
  INTEGER(KIND=LONG), PARAMETER                              :: nt_H   = 20  
  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
  REAL(KIND=DP), DIMENSION(Jmin_H:Jmax_H,Jmin_H:Jmax_H,nt_H) :: Rate_H
  ! temperatures which defines Rate_CO(i,j,k)
  REAL(KIND=DP), DIMENSION(nt_H)                             :: temp_H 
  ! new rates, calculated at Tn (useful collision rate array)
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                 :: CR_H
  ! name of the file
  CHARACTER(LEN=75)                                          :: filen1_H$ = 'input/coeff_12CO_H_inc.in'
  ! format of temperatures in the file
  CHARACTER(LEN=*), PARAMETER                                :: format_H_temp = '(5X,20F10.1)'  

  !--- CO-He collision rates (Cecchi-Pestellini et al., ApJ, 2002) ---
  ! max value of J for which collision rates are available
  INTEGER(KIND=LONG), PARAMETER                                   :: Jmin_He = 0, Jmax_He = 14 
  ! number of available temperatures for tabulated collision rates
  INTEGER(KIND=LONG), PARAMETER                                   :: nt_He   = 10  
  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
  REAL(KIND=DP), DIMENSION(Jmin_He:Jmax_He,Jmin_He:Jmax_He,nt_He) :: Rate_He
  ! temperatures which defines Rate_CO(i,j,k)
  REAL(KIND=DP), DIMENSION(nt_He)                                 :: temp_He
  ! new rates, calculated at Tn (useful collision rate array)
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                      :: CR_He
  ! name of the file
  CHARACTER(LEN=75)                                               :: filen1_He$ = 'input/coeff_12CO_He_inc.in'
  ! format of temperatures in the file
  CHARACTER(LEN=*), PARAMETER                                     :: format_He_temp = '(5X,10F9.1)'  

  !--------------------------------------------------
  ! source terms for CO level population (cm-3.s-1)
  !--------------------------------------------------
  ! change in CO levels number density (cm-3.s-1)
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YN_rot_CO         

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


CONTAINS

  SUBROUTINE READ_CO_LEVELS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads CO levels : J, weight, energy.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    CO_lev
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero, h
    USE MODULE_PHYS_VAR, ONLY: NCO_lev

    IMPLICIT NONE

    CHARACTER(len=1) :: charact
    INTEGER(KIND=LONG) :: i, ii, J

    IF (.not. ALLOCATED(CO_lev)) THEN
       ALLOCATE (CO_lev(0:NCO_lev-1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of CO_lev'
       END IF
    END IF

    ! initialization
    CO_lev(:)%J        = Zero
    CO_lev(:)%weight   = Zero
    CO_lev(:)%Energy   = Zero
    CO_lev(:)%fracpop  = Zero
    CO_lev(:)%fracbol  = Zero

    DO i = 0, NCO_lev-1
       ii = DBLE(i)
       CO_lev(i)%J      = i
       CO_lev(i)%weight = 2._DP * ii + 1._DP 
       CO_lev(i)%Energy = h * Brot_CO * 1.D6 * ii * (ii + 1._DP)
    END DO

 
  END SUBROUTINE READ_CO_LEVELS

  SUBROUTINE INITIALIZE_ROT_CO
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read collision rates for CO-oH2, CO-pH2, CO-H, and CO-He.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    rate_CO, temp_CO
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : Tn, NCO_lev
    USE MODULE_CONSTANTS, ONLY: pi, Zero, kB, h, c
    USE MODULE_DEBUG_JLB

    IMPLICIT NONE
    
    INTEGER(KIND=LONG)       :: j
    REAL(KIND=DP)            :: jj
    
    CST1 = h * Brot_CO * 1D6 / kB
    
    ! calculation of the partition function
    PF = 0._DP
    DO j = 0, NCO_lev - 1
       jj             = DBLE(j)
       bolj           = CST1 * (jj * (jj + 1._DP)) / Tn
       CO_lev(j)%DREL = CO_lev(j)%Weight * EXP(- bolj)
    END DO
    DO j = 0, NCO_lev - 1
       PF = PF + CO_lev(j)%DREL
    END DO
    
    ! calculation of Bolzmann's fractional populations
    DO j = 0, NCO_lev - 1
       CO_lev(j)%fracbol = CO_lev(j)%DREL / PF
       CO_lev(j)%fracpop = CO_lev(j)%fracbol
    END DO
    
    WHERE (CO_lev(:)%fracpop < population_limit)
       CO_lev(:)%fracpop = population_limit
    END WHERE


  END SUBROUTINE INITIALIZE_ROT_CO
  
  

  SUBROUTINE READ_CO_LINES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads CO lines (levels, energies, Aij).
    ! subroutine/function needed :
    !    NAME_CO_LINE
    ! input variables :
    ! ouput variables :
    ! results :
    !    Aij_CO, sum_Aij_CO
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero, h, kB
    USE MODULE_PHYS_VAR, ONLY : NCO_lev

    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i, ii, Jup, Jlow, NCO_up, NCO_low

    IF (.not. ALLOCATED(CO_lines)) THEN
       ALLOCATE (CO_lines(0:NCO_lev-2),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of CO_lines'
       END IF
    END IF
    IF (.not. ALLOCATED(A_CO)) THEN
       ALLOCATE (A_CO(0:39),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of A_CO'
       END IF
    END IF

    CST1 = h * Brot_CO * 1D6 / kB

    ! allocation and initialization of CO_lines
    CO_lines(:)%name       = ''
    CO_lines(:)%NCO_up     = Zero
    CO_lines(:)%NCO_low    = Zero
    CO_lines(:)%DeltaE     = Zero
    CO_lines(:)%BOLB       = Zero
    CO_lines(:)%BOLE       = Zero
    CO_lines(:)%Tdv        = Zero
    CO_lines(:)%Tdv2       = Zero
    CO_lines(:)%intensity  = Zero
    CO_lines(:)%intensity2 = Zero
    CO_lines(:)%Texc       = Zero
    CO_lines(:)%TAU        = Zero
    CO_lines(:)%EP1        = Zero
    CO_lines(:)%EP2        = Zero
    CO_lines(:)%A_CO       = Zero

    ! Leiden list of Einstein coefficients, s^-1
    A_CO(0)  = 7.203D-08 
    A_CO(1)  = 6.910D-07
    A_CO(2)  = 2.497D-06
    A_CO(3)  = 6.126D-06
    A_CO(4)  = 1.221D-05
    A_CO(5)  = 2.137D-05
    A_CO(6)  = 3.422D-05
    A_CO(7)  = 5.134D-05
    A_CO(8)  = 7.330D-05
    A_CO(9)  = 1.006D-04
    A_CO(10) = 1.339D-04
    A_CO(11) = 1.735D-04
    A_CO(12) = 2.200D-04
    A_CO(13) = 2.739D-04
    A_CO(14) = 3.354D-04
    A_CO(15) = 4.050D-04
    A_CO(16) = 4.829D-04
    A_CO(17) = 5.695D-04
    A_CO(18) = 6.650D-04
    A_CO(19) = 7.695D-04
    A_CO(20) = 8.833D-04
    A_CO(21) = 1.006D-03
    A_CO(22) = 1.139D-03
    A_CO(23) = 1.281D-03
    A_CO(24) = 1.432D-03
    A_CO(25) = 1.592D-03
    A_CO(26) = 1.761D-03
    A_CO(27) = 1.940D-03
    A_CO(28) = 2.126D-03
    A_CO(29) = 2.321D-03
    A_CO(30) = 2.524D-03
    A_CO(31) = 2.735D-03
    A_CO(32) = 2.952D-03
    A_CO(33) = 3.175D-03
    A_CO(34) = 3.404D-03
    A_CO(35) = 3.638D-03
    A_CO(36) = 3.878D-03
    A_CO(37) = 4.120D-03
    A_CO(38) = 4.365D-03
    A_CO(39) = 4.613D-03

    DO i = 0,NCO_lev - 2
       CO_lines(i)%A_CO    = A_CO(i)
       CO_lines(i)%NCO_up  = i+1
       CO_lines(i)%NCO_low = i
       CO_lines(i)%DeltaE  = CO_lev(i+1)%energy-CO_lev(i)%energy
       ii = DBLE(i)
       CO_lines(i)%BOLB    = 1._DP / (exp(2._DP * (ii + 1._DP) * CST1 / Tb) - 1._DP)
    END DO


  END SUBROUTINE READ_CO_LINES



  SUBROUTINE NAME_CO_LINES 
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_H2_LINES
    ! purpose :
    !    computes the name of one line of the form (1-0)
    ! subroutine/function needed :
    ! input variables :
    !    * Vup, Jup  -> (V,J) of the upper level  (two digits : < 100)
    !    * Jlow      -> J of th lower level       (two digits : < 100)
    !    * line_type -> 'S', 'O', 'Q'
    ! output variables :
    ! results :
    !    name -> name of the line, without blanks
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : NCO_lev

    IMPLICIT NONE
    INTEGER(KIND=LONG) :: Jlow, Jup
    CHARACTER(LEN=2) :: str_Jup, str_Jlow
    CHARACTER(LEN=9) :: name

    DO j = 0, NCO_lev - 2
       Jup  = j + 1
       Jlow = j
       WRITE(str_Jup,'(I2)')Jup
       WRITE(str_Jlow,'(I2)')Jlow
       name = '(' // TRIM(ADJUSTL(str_Jup)) // '-' // TRIM(ADJUSTL(str_Jlow)) // ')'
       CO_lines(j)%name = name
!       write(*,*) CO_lines(j)%name
    END DO

  END SUBROUTINE NAME_CO_LINES



  SUBROUTINE EVOLUTION_CO
    !---------------------------------------------------------------------------
    ! called by :
    !    COMPUTE_CO
    ! purpose :
    !    Calculate the source term (cm-3.s-1) for rotational levels of CO,
    !    using Aij_CO (radiation, s-1) read in READ_CO_LINE and Cij_CO
    !    (collisions, s-1) calculated in this subroutine.
    !    The result, YN_rot_CO, is used in COMPUTE_CO and in DIFFUN.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rot_CO
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY         : Tn, NCO_lev, distance, vsound, dVn, vn, op_H2
    USE MODULE_CHEMICAL_SPECIES, ONLY : Dens_H, Dens_H2, Dens_He, Dens_CO
    USE MODULE_H2, ONLY               : Dens_orthoH2, Dens_paraH2
    USE MODULE_CONSTANTS, ONLY        : pi, Zero, kB, h, c
    USE MODULE_DEBUG_JLB
    USE MODULE_INITIALIZE
!    USE MODULE_PHYS_VAR, ONLY         : NCO_lev

    IMPLICIT NONE

!    REAL(KIND=DP) :: distance, Tn, vsound, dVn, vn,op_H2
!    REAL(KIND=DP) :: Dens_H, Dens_He, Dens_CO, Dens_H2, Dens_orthoH2, Dens_paraH2

!    distance = 6.4768D14
!    Tn       = 8.1733_DP
!    vsound   = 21850._DP
!    dVn      = 1.4482D-14
!    vn       = 2.D6
!    op_H2    = 3._DP
!    Dens_H   = 3.0582_DP
!    Dens_He  = 9.9997D3
!    Dens_H2  = 4.9998D4
!    Dens_CO  = 8.2297
!    Dens_orthoH2 = 3.74985D4
!    Dens_paraH2  = 1.24995D4


    CST1 = h * Brot_CO * 1D6 / kB

    CALL NAME_CO_LINES

    ! calculation of the partition function
    PF = 0._DP
    DO j = 0, NCO_lev - 1
       jj             = DBLE(j)
       bolj           = CST1 * (jj * (jj + 1._DP)) / Tn
       CO_lev(j)%DREL = CO_lev(j)%Weight * EXP(- bolj)
    END DO
    DO j = 0, NCO_lev - 1
       PF = PF + CO_lev(j)%DREL
    END DO
    ! calculation of Bolzmann's fractional populations
    DO j = 0, NCO_lev - 1
       CO_lev(j)%fracbol = CO_lev(j)%DREL / PF
    END DO

    IF (.not. ALLOCATED(CR_oH2)) THEN
       ALLOCATE (CR_oH2(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of CR_oH2'
       END IF
    END IF
    IF (.not. ALLOCATED(CR_pH2)) THEN
       ALLOCATE (CR_pH2(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of CR_pH2'
       END IF
    END IF
    IF (.not. ALLOCATED(CR_H)) THEN
       ALLOCATE (CR_H(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of CR_H'
       END IF
    END IF
    IF (.not. ALLOCATED(CR_He)) THEN
       ALLOCATE (CR_He(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of CR_He'
       END IF
    END IF
    IF (.not. ALLOCATED(Coll_oH2)) THEN
       ALLOCATE (Coll_oH2(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of Coll_oH2'
       END IF
    END IF
    IF (.not. ALLOCATED(Coll_pH2)) THEN
       ALLOCATE (Coll_pH2(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of Coll_pH2'
       END IF
    END IF
    IF (.not. ALLOCATED(Coll_H)) THEN
       ALLOCATE (Coll_H(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of Coll_H'
       END IF
    END IF
    IF (.not. ALLOCATED(Coll_He)) THEN
       ALLOCATE (Coll_He(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of Coll_He'
       END IF
    END IF
    IF (.not. ALLOCATED(Coll)) THEN
       ALLOCATE (Coll(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of Coll'
       END IF
    END IF
    IF (.not. ALLOCATED(Radi)) THEN
       ALLOCATE (Radi(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of Radi'
       END IF
    END IF
    IF (.not. ALLOCATED(MATRIX)) THEN
       ALLOCATE (MATRIX(0:NCO_lev - 1,0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of MATRIX'
       END IF
    END IF
    IF (.not. ALLOCATED(YN_rot_CO)) THEN
       ALLOCATE (YN_rot_CO(0:NCO_lev - 1),stat=err)
       IF (err /= 0) THEN
          WRITE(*,*) 'Error in allocation of YN_rot_CO'
       END IF
    END IF

    CR_oH2 = Zero
    CR_pH2 = Zero
    CR_H   = Zero
    CR_He  = Zero

!    Tn = 300._DP

    ! calls the interpolation procedure for the rate coefficients
    CALL INTERPOLATION(Jmin_oH2, Jmax_oH2, Rate_oH2, temp_oH2, CR_oH2, nt_oH2)
    CALL INTERPOLATION(Jmin_pH2, Jmax_pH2, Rate_pH2, temp_pH2, CR_pH2, nt_pH2)
    CALL INTERPOLATION(Jmin_H  , Jmax_H  , Rate_H  , temp_H  , CR_H  , nt_H  )
    CALL INTERPOLATION(Jmin_He , Jmax_He , Rate_He , temp_He , CR_He , nt_He )

    write(*,*) 'coucou'

!    write(*,*) 'CR_oH2(3,5)=',CR_oH2(3,5)
!    write(*,*) 'CR_oH2(19,32)=',CR_oH2(19,32)

!    write(*,*) 'CR_pH2(5,3)=',CR_pH2(5,3)
!    write(*,*) 'CR_pH2(32,19)=',CR_pH2(32,19)

    ! minimum velocity gradient based on sound velocity
    gradv = sqrt(dVn ** 2 + vsound ** 2 / (distance + 1D13) ** 2)

    CST2 = 8._DP * ((pi) ** 3) * ((Moment_CO * 1D-18) ** 2) / (3._DP * h)
    CST3 = CST2 * Dens_CO / abs(gradv)

    CALL ESCAPE_PROBABILITY(moment_CO,gradv,Dens_CO,CO_lev%fracpop,CO_lines%TAU,CO_lines%EP1,CO_lines%EP2,NCO_lev - 2,choice_ep)
    CALL Rad_Matrix
    CALL Coll_Matrix

    Coll = Dens_orthoH2 * Coll_oH2  + &
           Dens_paraH2  * Coll_pH2  + &
           Dens_H       * Coll_H    + &
           Dens_He      * Coll_He

    ! fill the global matrix
    MATRIX = 0._DP
    MATRIX = Coll + Radi
    
!    write(*,*) 'Matrix(3,5)=',Matrix(3,5)
!    write(*,*) 'Matrix(5,3)=',Matrix(5,3)
!    write(*,*) 'Matrix(19,32)=',Matrix(19,32)
!    write(*,*) 'Matrix(32,19)=',Matrix(32,19)
!    write(*,*) 'Matrix(5,5)=',Matrix(5,5)
!    write(*,*) 'Matrix(25,25)=',Matrix(25,25)

!    write(*,*) 'Dens_orthoH2=', Dens_orthoH2
!    write(*,*) 'Dens_paraH2=', Dens_paraH2
!    write(*,*) 'Dens_H=', Dens_H
!    write(*,*) 'Dens_He=', Dens_He

!    write(*,*) 'Coll(5,5)=', Coll(5,5)
!    write(*,*) 'Coll_oH2(5,5)=', Coll_oH2(5,5)
!    write(*,*) 'Coll_pH2(5,5)=', Coll_pH2(5,5)
!    write(*,*) 'Coll_H(5,5)=', Coll_H(5,5)
!    write(*,*) 'Coll_He(5,5)=', Coll_He(5,5)
!    stop

!    MATRIX(NCO_lev - 1,:) = 1._DP

!    DO j = 0, 2
!       write(*,*) MATRIX(j,:)
!    END DO

!       write(*,*) 'distance =',     distance
!       write(*,*) 'vsound =',       vsound
!       write(*,*) 'Tn = ',          Tn
!       write(*,*) 'dVn=',           dVn
!       write(*,*) 'gradv =',        gradv
!       write(*,*) 'Tot_dens_CO =',  Dens_CO
!       write(*,*) 'op_H2 =',        op_H2
!       write(*,*) 'n_H2 =',         Dens_H2
!       write(*,*) 'n_He =',         Dens_He
!       write(*,*) 'n_H =',          Dens_H
!       write(*,*) 'n_pH2 =',        Dens_paraH2
!       write(*,*) 'n_oH2 =',        Dens_orthoH2

    ! calculate evolution terms for populations in s-1
    !    = matrix multiplication of population density and
    !    (de)excitation probabilities
    YN_rot_CO = MATMUL(Matrix,CO_lev%fracpop)

  END SUBROUTINE EVOLUTION_CO


  SUBROUTINE READ_CO_RATES(filen1$,format_CO_temp,rate_CO,temp_CO,nt,Jmin_rate,Jmax_rate)
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read collision rates for CO-oH2, CO-pH2, CO-H, and CO-He.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    rate_CO, temp_CO
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE

    INTEGER(KIND=LONG)                                                    :: nt, Jmin_rate, Jmax_rate, Jmax_rate_2
    INTEGER(KIND=LONG)                                                    :: i, Ji, Jf, error, Nrates, nch1
    CHARACTER(LEN=75)                                                     :: name_file_CO
    CHARACTER(LEN=75)                                                     :: format_CO_temp
    CHARACTER(LEN=75)                                                     :: filen1$
    INTEGER                                                               :: file_CO
    REAL(KIND=DP), DIMENSION(Jmin_rate:Jmax_rate,Jmin_rate:Jmax_rate, nt) :: Rate_CO
    REAL(KIND=DP), DIMENSION(nt)                                          :: Temp_CO

    DO i = 1,LEN(filen1$)
       IF (filen1$(i:i) /= '') THEN
          nch1 = i
       END IF
    END DO
    filen1$ = filen1$(1:nch1)

    name_file_CO = filen1$(1:nch1)

    ! initialization
    Rate_CO = 0._DP
    Temp_CO = 0._DP

    ! file opening
    file_CO = GET_FILE_NUMBER()
    OPEN(file_CO,file=name_file_CO,STATUS='OLD',  &
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! comments
    DO i = 1, 8
       READ(file_CO,*)
    END DO

    ! temperatures
    READ(file_CO,format_CO_temp) temp_CO

    ! comments
    DO i = 1,3
       READ(file_CO,*)
    END DO

    ! collision rate coefficients
    Nrates = SIZE(Rate_CO,dim=1)*SIZE(Rate_CO,dim=2)
    error = 0
    DO i = 1, Nrates
       READ(file_CO,*,iostat=error)Ji,Jf,Rate_CO(Ji-1,Jf-1,:)
       IF (error>0) STOP "*** WARNING, error in READ_CO_RATES"
    END DO

    CLOSE(file_CO)

  END SUBROUTINE READ_CO_RATES



  SUBROUTINE INTERPOLATION(Jmin_rate,Jmax_rate,Rate_CO,Temp_CO,CR,nt)
    !---------------------------------------------------------------------------
    ! purpose :
    !    interpolates the value of Tkin among the ones for which the collision
    !    rates are available
    !    method = simple linear interpolation
    ! subroutine/function needed :
    ! input variables :
    !   Rate_CO, temp_CO, tkin 
    ! ouput variables :
    ! results :
    !   NRate_CO, CR
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : h, kB
    USE MODULE_PHYS_VAR, ONLY : Tn

    IMPLICIT NONE

    INTEGER(KIND=LONG)                                                    :: ji,jf,ilo,ihi,it,nt
    INTEGER(KIND=LONG)                                                    :: Jmin_rate, Jmax_rate, Jmax_rate_2
    REAL(KIND=DP)                                                         :: xp, CST1
    REAL(KIND=DP)                                                         :: DE, R, gl, gu, jji, jjf
    REAL(KIND=DP), DIMENSION(Jmin_rate:Jmax_rate,Jmin_rate:Jmax_rate,nt)  :: Rate_CO
    REAL(KIND=DP), DIMENSION(nt)                                          :: Temp_CO
    REAL(KIND=DP), DIMENSION(0:NCO_lev - 1,0:NCO_lev - 1)                 :: CR
    
    CST1 = h * Brot_CO * 1.D6 / kB

    ! initializations
    ilo = 1
    ihi = nt
    CR  = 0._DP

    ! determination of the nearest temperatures (higher and lower) in the array of
    ! available collision rates that brackett the kinetic temperature
    ! indexes : ilo, ihi
    it = nt+1
    DO WHILE (ilo /= it)
       it = it-1
       IF (temp_CO(it) <= Tn) THEN
          ilo = it
       END IF
    END DO
    it = 0
    DO WHILE (ihi /= it)
       it = it+1
       IF (temp_CO(it) >= Tn) THEN
          ihi = it
       END IF
    END DO
    
    Jmax_rate_2 = MIN(Jmax_rate,NCO_lev - 1)

    ! first case : the kinetic temperature belongs to the list of temperatures for which 
    ! collision rates are available
    ! second case : the interpolation is necessary
    ! just use the deexcitation rates and detailed balance
    IF (ihi == ilo) THEN
       DO ji = 0, Jmax_rate_2
          DO jf = 0, Jmax_rate_2
             IF (jf < ji) THEN
                jji = DBLE(ji)
                jjf = DBLE(jf)
                CR(jf,ji) = Rate_CO(ji,jf,ilo)
                DE = CST1 * (ji * (ji + 1) - jf * (jf + 1)) / Tn
                gl = 2._DP * jjf + 1._DP
                gu = 2._DP * jji + 1._DP
                R = gl /gu
                CR(ji,jf) = CR(jf,ji) * exp(- DE) / R
             END IF
          END DO
       END DO
    ELSE
       DO ji = 0, Jmax_rate_2
          DO jf = 0, Jmax_rate_2
             IF (jf < ji) THEN
                jji = DBLE(ji)
                jjf = DBLE(jf)
                CR(jf,ji) = Rate_CO(ji,jf,ilo) + (Tn - temp_CO(ilo)) * (Rate_CO(ji,jf,ihi) - &
                            Rate_CO(ji,jf,ilo))/ (temp_CO(ihi) - temp_CO(ilo))
                DE = CST1 * (ji * (ji + 1) - jf * (jf + 1)) / Tn
                gl = 2._DP * jjf + 1._DP
                gu = 2._DP * jji + 1._DP
                R = gl /gu
                CR(ji,jf) = CR(jf,ji) * exp(- DE) / R
             END IF
          END DO
       END DO
    END IF
    

  END SUBROUTINE INTERPOLATION

  SUBROUTINE ESCAPE_PROBABILITY(moment,grad,dens,pop,TAU,EP1,EP2,jmax,expr)
    !---------------------------------------------------------------------------
    ! purpose :
    !    calculate the optical depth and escape probability of a line
    ! subroutine/function needed :
    !    EINSTEIN_COEFF_SiO
    ! input variables :
    !   vgrad, dens_SiO
    ! ouput variables :
    ! results :
    !   TAU, EP1, EP2
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : h, kB, pi
    USE MODULE_PHYS_VAR, ONLY  : NCO_lev
    
    IMPLICIT NONE

    INTEGER(KIND=LONG) :: jmax,expr
    REAL(KIND=DP)      :: grad,moment,dens
    REAL(KIND=DP), DIMENSION(0:NCO_lev - 1), INTENT(IN)    :: pop
    REAL(KIND=DP), DIMENSION(0:NCO_lev - 2), INTENT(OUT) :: TAU
    REAL(KIND=DP), DIMENSION(0:NCO_lev - 2), INTENT(OUT) :: EP1, EP2

    CST2 = 8._DP * ((pi) ** 3) * ((Moment_CO * 1D-18) ** 2) / (3._DP * h)
    CST3 = CST2 * dens / abs(grad)

    DO j = 0, jmax
       jj = DBLE(j) ! conversion to real
       TAU(J) = CST3 * (jj + 1._DP) * (pop(j)/(2._DP * jj + 1._DP) - pop(j+1)/(2._DP * jj + 3._DP))
!       write(*,*) 'TAU(',j,')=',TAU(j)

       IF (TAU(J) <= 1.D-20) TAU(J) = 1.D-20

       IF (expr == 1) THEN

!          IF (TAU(J)<=0.0001) THEN
!             EP1(J) = 1._DP - 0.5_DP * TAU(J)
!             EP1(J) = 1._DP - 3._DP * TAU(J)
!          ELSE IF (TAU(J) > 0.00001 .and. TAU(J) <= 10000) THEN
!             EP1(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)      ! same direction as the v_grad, Surdej 1977
!             EP1(J) = 1._DP / (1._DP + 3._DP * TAU(J))
!          ELSE
!             EP1(J) = 1._DP / TAU(J)
!             EP1(J) = 1._DP
!          END IF
!          IF (TAU(J) > 100 .and. TAU(J) < -100) THEN
!             EP1(J) = 1._DP / TAU(J)
!          ELSE IF (TAU(J) >= 1D-4 .and. TAU(J) <= 100) THEN
!             EP1(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)
!          ELSE IF (TAU(J) <= -1D-4 .and. TAU(J) >= -100) THEN
!             EP1(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)
!          ELSE IF (TAU(J) >= 0 .and. TAU(J) < 1D-4) THEN
!             EP1(J) = 1 - TAU(J) / 2._DP + TAU(J) ** 2 / 6._DP & 
!                        - TAU(J) ** 3 / 24._DP + TAU(J) ** 4 / 120._DP
!         ELSE
!             EP1(J) = 1 - TAU(J) / 2._DP - TAU(J) ** 2 / 6._DP & 
!                        - TAU(J) ** 3 / 24._DP - TAU(J) ** 4 / 120._DP
!         END IF



!          EP2(J) = EP1(J)

!          EP1(J) = 1._DP / (1._DP + 3._DP * TAU(J))             ! plan parallel expr, Neufeld & Kaufman 1993
!          IF (TAU(J) <= 0.0001) THEN
!             EP1(j) = 1._DP - 3._DP * TAU(J)
!             EP1(j) = 1._DP - 0.5_DP * TAU(j)
!          ELSE IF (TAU(J) > 0.0001 .and. TAU(J) <= 100) THEN
!             EP1(J) = 1._DP / (1._DP + 3._DP * TAU(J))
!             EP1(j) = (1._DP - exp(-abs(TAU(j)))) / TAU(j)
!          ELSE
!             EP1(j) = 0._DP
!             EP1(j) = 1._DP / TAU(j)
!          END IF 

          EP1(J) = 1._DP
          EP2(J) = 1._DP

!          EP2(J) = EP1(J)

!          IF (TAU(J)<=0.0001) THEN
!             EP2(j) = 1._DP - 0.5_DP * TAU(j)
!          ELSE IF (TAU(J) > 0.0001 .and. TAU(J) <= 100) THEN
!             EP2(j) = (1._DP - exp(-abs(TAU(j)))) / TAU(j)                 ! classical expression, Surdej 1977
!          ELSE
!             EP2(j) = 1._DP / TAU(j)
!          END IF 


!         IF (TAU(J)<=0.0001) THEN
!            EP2(J) = 1._DP - 0.5_DP * TAU(J)
!         ELSE IF (TAU(J) > 0.0001 .and. TAU(J) <= 100) THEN
!            EP2(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)      ! same direction as the v_grad, Surdej 1977
!         ELSE
!            EP2(J) = 1._DP / TAU(J)
!         END IF


!         IF (TAU(J) > 100 .and. TAU(J) < -100) THEN
!            EP2(J) = 1._DP / TAU(J)
!         ELSE IF (TAU(J) >= 1D-4 .and. TAU(J) <= 100) THEN
!            EP2(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)
!         ELSE IF (TAU(J) <= -1D-4 .and. TAU(J) >= -100) THEN
!            EP2(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)
!         ELSE IF (TAU(J) >= 0 .and. TAU(J) < 1D-4) THEN
!            EP2(J) = 1 - TAU(J) / 2._DP + TAU(J) ** 2 / 6._DP & 
!                       - TAU(J) ** 3 / 24._DP + TAU(J) ** 4 / 120._DP
!         ELSE
!            EP2(J) = 1 - TAU(J) / 2._DP - TAU(J) ** 2 / 6._DP & 
!                       - TAU(J) ** 3 / 24._DP - TAU(J) ** 4 / 120._DP
!         END IF
!       ELSE IF (expr == 2) THEN                                                         ! Surdej expr
!          IF (TAU(J) > 100 .and. TAU(J) < -100) THEN
!             EP1(J) = 1._DP / TAU(J)
!             EP2(J) = EP1(J)
!          ELSE IF (TAU(J) >= 1D-4 .and. TAU(J) <= 100) THEN
!             EP1(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)
!             EP2(J) = EP1(J)
!          ELSE IF (TAU(J) <= -1D-4 .and. TAU(J) >= -100) THEN
!             EP1(J) = (1._DP - exp(-abs(TAU(J)))) / TAU(J)
!             EP2(J) = EP1(J)
!          ELSE IF (TAU(J) >= 0 .and. TAU(J) < 1D-4) THEN
!             EP1(J) = 1 - TAU(J) / 2._DP + TAU(J) ** 2 / 6._DP & 
!                        - TAU(J) ** 3 / 24._DP + TAU(J) ** 4 / 120._DP
!             EP2(J) = EP1(J)
!          ELSE
!             EP1(J) = 1 - TAU(J) / 2._DP - TAU(J) ** 2 / 6._DP & 
!                        - TAU(J) ** 3 / 24._DP - TAU(J) ** 4 / 120._DP
!             EP2(J) = EP1(J)
!          END IF
!       ELSE IF (expr == 3) THEN
!          IF (TAU(J)<=0.0001) THEN
!             EP1(J) = 1._DP - 1.5_DP * TAU(J)
!             EP2(J) = EP1(J)
!          ELSE IF (TAU(J) > 0.0001 .and. TAU(J) <= 100) THEN
!             EP1(J) = (1._DP - exp(-abs(3._DP * TAU(J)))) / (3._DP * TAU(J))   ! classical 'slab' expr
!             EP2(J) = EP1(J)
!          ELSE
!             EP1(J) = 1._DP / (3._DP * TAU(J))
!             EP2(J) = EP1(J)
!          END IF
!!       ELSE IF (expr == 4) THEN                                                         ! Radex expression
!!          IF (TAU(J) == 0._DP) THEN
!!             EP1(J) = 1._DP
!!             EP2(J) = 1._DP
!          IF (ABS(TAU(J)) <= 1.D-20) THEN
!             EP1(J) = 1._DP
!             EP2(J) = 1._DP
!!          ELSE IF (TAU(J) >= -1.D-05 .and. TAU(J) < 0) THEN
!          ELSE IF (TAU(J) >= -1.D-04 .and. TAU(J) < - 1.D-20) THEN
!             EP1(J) = 2._DP + 6._DP / TAU(J) ** 2 &
!                            + 6._DP / TAU(J) &
!                            + 5._DP * TAU(J) / 6._DP & 
!                            + 3._DP * TAU(J) ** 2 / 20._DP  ! &
!!                            + 7._DP * TAU(J) ** 3 / 240._DP & 
!!                            + TAU(J) ** 4 / 210._DP &
!!                            + 3._DP * TAU(J) ** 5 / 4480._DP & 
!!                            + TAU(J) ** 6 / 12096._DP & 
!!                            + 11._DP * TAU(J) ** 7 / 1209600._DP & 
!!                            + TAU(J) ** 8 / 1108800._DP & 
!!                            + 13._DP * TAU(J) ** 9 / 159667200._DP & 
!!                            + TAU(J) ** 10 / 148262400._DP &
!!                            + TAU(J) ** 11 / 1937295360._DP & 
!!                            + TAU(J) ** 12 / 27243216000._DP &
!!                            + TAU(J) ** 13 / 435891456000._DP 
!             EP2(J) = EP1(J)
!!          ELSE  IF (TAU(J) <= 1.D-05 .and. TAU(J) > 0) THEN
!          ELSE  IF (TAU(J) <= 1.D-04 .and. TAU(J) > 1.D-20) THEN
!             EP1(J) = 1._DP - 3._DP * TAU(J) / 8._DP & 
!                            + TAU(J) ** 2 / 10._DP  ! & 
!!                            - TAU(J) ** 3 / 48._DP & 
!!                            + TAU(J) ** 4 / 280._DP &
!!                            - TAU(J) ** 5 / 1920._DP & 
!!                            + TAU(J) ** 6 / 15120._DP & 
!!                            - TAU(J) ** 7 / 134400._DP &
!!                            + TAU(J) ** 8 / 1330560._DP & 
!!                            - TAU(J) ** 9 / 14515200._DP & 
!!                            + TAU(J) ** 10 / 172972800._DP &
!!                            - TAU(J) ** 11 / 2235340800._DP & 
!!                            + TAU(J) ** 12 / 31135104000._DP &
!!                            - TAU(J) ** 13 / 435891456000._DP 
 !            EP2(J) = EP1(J)
 !         ELSE IF (TAU(J) > 1.D-04 .and. TAU(J) <= 100000) THEN
 !            EP1(J) = (1.5_DP / TAU(J)) * (1._DP - 2._DP / TAU(J) ** 2 & 
 !                              + (2._DP / TAU(J) + 2._DP / TAU(J) ** 2) * exp(-abs(TAU(J))))
 !            EP2(J) = EP1(J)
 !         ELSE IF (TAU(J) < -1.D-04 .and. TAU(J) > -100000) THEN
 !            EP1(J) = (1.5_DP / TAU(J)) * (1._DP - 2._DP / TAU(J) ** 2 & 
 !                              + (2._DP / TAU(J) + 2._DP / TAU(J) ** 2) * exp(-abs(TAU(J))))
 !            EP2(J) = EP1(J)
 !         ELSE
 !            EP1(J) = 3._DP / (2._DP * TAU(J))
 !            EP2(J) = EP1(J)
 !         END IF
       ELSE
          WRITE(*,*) 'erreur dans le choix de votre expression de EP1'
          WRITE(*,*) 'votre choix ne peut etre que 1, 2, 3, ou 4'
      END IF
!      write(*,*) 'EP1(',j,')=',EP1(j)
!      write(*,*) 'EP2(',j,')=',EP2(j)
    END DO

!    stop

  END SUBROUTINE ESCAPE_PROBABILITY  

  SUBROUTINE RAD_MATRIX
    !---------------------------------------------------------------------------
    ! purpose :
    !    construct the radiative part of the statistical equilibrium matrix
    ! subroutine/function needed : 
    !    EINSTEIN_COEFF_SiO, B_nu
    ! input variables :
    ! ouput variables :
    ! results :
    !   Radi
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : NCO_lev

    IMPLICIT NONE
    
    INTEGER :: i, j

    ! initialization
    Radi = 0._DP

    ! particular cases : J = 0 and J = NCO_lev - 1
    Radi(0,0)                     = - CO_lines(0)%EP1 * CO_lines(0)%A_CO * 3._DP * CO_lines(0)%BOLB
    Radi(NCO_lev - 1,NCO_lev - 1) = - CO_lines(NCO_lev - 2)%EP1 * CO_lines(NCO_lev - 2)%A_CO * (1._DP + CO_lines(NCO_lev - 2)%BOLB)

!    write(*,*) 'Radi(0,0)=', Radi(0,0)
!    write(*,*) 'Radi(NCO_lev - 1,NCO_lev - 1)=',Radi(NCO_lev - 1,NCO_lev - 1)
!    write(*,*) 'NCO_lev - 1 =', NCO_lev - 1

!    stop

    ! diagonal
    DO j = 1,NCO_lev - 2
       jj = DBLE(j)
       Radi(j,j) = - CO_lines(j)%EP1 * CO_lines(j)%A_CO * ((2._DP * jj + 3._DP) /  (2._DP * jj + 1._DP)) &
                   * CO_lines(j)%BOLB - CO_lines(j - 1)%EP1 * CO_lines(j - 1)%A_CO * (1._DP + CO_lines(j - 1)%BOLB) 
!       write(*,*) 'Radi(',j,',',j,')=',Radi(j,j)
    END DO

!    stop

    ! above the diagonal
    DO j = 0,NCO_lev - 2
       Radi(j,j+1) = CO_lines(j)%EP1 * CO_lines(j)%A_CO * (1._DP + CO_lines(j)%BOLB)
!       write(*,*) 'Radi(',j,',',j+1,')=',Radi(j,j+1)
    END DO

!    stop

    ! below the diagonal
    DO j = 1,NCO_lev - 1 
       jj = DBLE(j)
       Radi(j,j-1) = CO_lines(j - 1)%EP1 * CO_lines(j - 1)%A_CO * CO_lines(j - 1)%BOLB &
                                         * (2._DP * jj + 1._DP) / (2._DP * jj - 1._DP)
!       write(*,*) 'Radi(',j,',',j-1,')=',Radi(j,j-1)
    END DO

!    stop

  END SUBROUTINE RAD_MATRIX


  SUBROUTINE COLL_MATRIX
    !---------------------------------------------------------------------------
    ! purpose :
    !    construct the collisionnal part of the statistical equilibrium matrix
    ! subroutine/function needed : 
    !    EINSTEIN_COEFF_SiO, B_nu
    ! input variables :  
    ! ouput variables :
    ! results :
    !   Coll
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : NCO_lev

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND = DP) :: CR_limit

    CR_limit = 1.D-30

    ! initialization
    Coll_oH2 = 0._DP
    Coll_pH2 = 0._DP
    Coll_H   = 0._DP
    Coll_He  = 0._DP

    ! simplest cases
    DO i = 0, NCO_lev - 1
       DO j = 0, NCO_lev - 1
          Coll_oH2(i,j) = CR_oH2(i,j)
          Coll_pH2(i,j) = CR_pH2(i,j)
          Coll_H(i,j)   = CR_H(i,j)
          Coll_He(i,j)  = CR_He(i,j)
       END DO
    END DO

!    write(*,*) 'Coll_oH2(5,4)',Coll_oH2(5,4)
!    write(*,*) 'Coll_oH2(4,5)',Coll_oH2(4,5)
!    write(*,*) 'Coll_oH2(32,19)',Coll_oH2(32,19)
!    write(*,*) 'Coll_oH2(19,32)',Coll_oH2(19,32)
!    write(*,*) 'Coll_pH2(5,4)',Coll_pH2(5,4)
!    write(*,*) 'Coll_pH2(4,5)',Coll_pH2(4,5)
!    write(*,*) 'Coll_pH2(32,19)',Coll_pH2(32,19)
!    write(*,*) 'Coll_pH2(19,32)',Coll_pH2(19,32)
!    write(*,*) 'Coll_H(5,4)',Coll_H(5,4)
!    write(*,*) 'Coll_H(4,5)',Coll_H(4,5)
!    write(*,*) 'Coll_H(32,19)',Coll_H(32,19)
!    write(*,*) 'Coll_H(19,32)',Coll_H(19,32)
!    write(*,*) 'Coll_He(5,4)',Coll_He(5,4)
!    write(*,*) 'Coll_He(4,5)',Coll_He(4,5)
!    write(*,*) 'Coll_He(32,19)',Coll_He(32,19)
!    write(*,*) 'Coll_He(19,32)',Coll_He(19,32)

!    stop

    ! diagonal
    DO j = 0,NCO_lev - 1
       Coll_oH2(j,j) = - SUM(CR_oH2(:,j))
       Coll_pH2(j,j) = - SUM(CR_pH2(:,j))
       Coll_H(j,j)   = - SUM(CR_H(:,j))
       Coll_He(j,j)  = - SUM(CR_He(:,j))
    END DO

!    write(*,*) 'Coll_He(0,0)=',Coll_He(0,0)
!    write(*,*) 'Coll_He(5,5)=',Coll_He(5,5)
!    write(*,*) 'Coll_He(10,10)=',Coll_He(10,10)
!    write(*,*) 'Coll_He(15,15)=',Coll_He(15,15)
!    write(*,*) 'Coll_He(20,20)=',Coll_He(20,20)
!    write(*,*) 'Coll_He(25,25)=',Coll_He(25,25)
!    write(*,*) 'Coll_He(30,30)=',Coll_He(30,30)
!    write(*,*) 'Coll_He(35,35)=',Coll_He(35,35)

!    stop

    ! above the diagonal(1)
    DO j = 0,NCO_lev - 2
       Coll_oH2(j,j+1) = CR_oH2(j,j+1)
       Coll_pH2(j,j+1) = CR_pH2(j,j+1)
       Coll_H(j,j+1)   = CR_H(j,j+1)
       Coll_He(j,j+1)  = CR_He(j,j+1)
    END DO

!    write(*,*) 'Coll_pH2(0,1)=',Coll_pH2(0,1)
!    write(*,*) 'Coll_pH2(5,6)=',Coll_pH2(5,6)
!    write(*,*) 'Coll_pH2(10,11)=',Coll_pH2(10,11)
!    write(*,*) 'Coll_pH2(15,16)=',Coll_pH2(15,16)
!    write(*,*) 'Coll_pH2(20,21)=',Coll_pH2(20,21)
!    write(*,*) 'Coll_pH2(25,26)=',Coll_pH2(25,26)
!    write(*,*) 'Coll_pH2(30,31)=',Coll_pH2(30,31)
!    write(*,*) 'Coll_pH2(35,36)=',Coll_pH2(35,36)

!    stop

    ! below the diagonal(1)
    DO j = 1,NCO_lev - 1
       Coll_oH2(j,j-1) = CR_oH2(j,j-1)
       Coll_pH2(j,j-1) = CR_pH2(j,j-1)
       Coll_H(j,j-1)   = CR_H(j,j-1)
       Coll_He(j,j-1)   = CR_He(j,j-1)
    END DO

!    write(*,*) 'Coll_He(1,0)=',Coll_He(1,0)
!    write(*,*) 'Coll_He(6,5)=',Coll_He(6,5)
!    write(*,*) 'Coll_He(11,10)=',Coll_He(11,10)
!    write(*,*) 'Coll_He(16,15)=',Coll_He(16,15)
!    write(*,*) 'Coll_He(21,20)=',Coll_He(21,20)
!    write(*,*) 'Coll_He(26,25)=',Coll_He(26,25)
!    write(*,*) 'Coll_He(31,30)=',Coll_He(31,30)
!    write(*,*) 'Coll_He(36,35)=',Coll_He(36,35)

!    write(*,*) 'Coll_oH2(5,5)=', Coll_oH2(5,5)

!    stop

    WHERE (abs(Coll_oH2(:,:)) < CR_limit)
       Coll_oH2(:,:) = CR_limit
    END WHERE
    WHERE (abs(Coll_pH2(:,:)) < CR_limit)
       Coll_pH2(:,:) = CR_limit
    END WHERE

  END SUBROUTINE COLL_MATRIX  


  SUBROUTINE COMPUTE_CO
    !---------------------------------------------------------------------------
    ! called by :
    !    SOURCE
    ! purpose :
    !    Calculate the source terms for CO levels (cm-3.s-1)
    ! subroutine/function needed :
    !    EVOLUTION_CO
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rot_CO
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY        : kB,h,c,pi
    USE MODULE_CHEMICAL_SPECIES, ONLY : Dens_CO
    USE MODULE_PHYS_VAR, ONLY         : NCO_lev,dVn
    USE MODULE_DEBUG_JLB
    IMPLICIT NONE

    !----------------------------------------------------------
    ! In EVOLUTION_CO, the source terms YN_rot_CO (cm-3.s-1)
    ! is calculated, with treatment of radiation and collision
    !----------------------------------------------------------
    CALL EVOLUTION_CO

    !----------------------------------------------------------
    ! Calculate the intensity (K) and excitation temperature (K) 
    ! of each CO rotational line.
    ! result in CO_lines%intensity, CO_lines%intensity2
    !----------------------------------------------------------

    CST1 = h * Brot_CO * 1D6 / kB
    CST4 = h * (c ** 3) / (8._DP * kB * 1D12 * pi)

!    WHERE (CO_lev(:)%fracpop < population_limit)
!       CO_lev(:)%fracpop = population_limit
!    END WHERE

    CALL ESCAPE_PROBABILITY(moment_CO,gradv,Dens_CO,CO_lev%fracpop,CO_lines%TAU,CO_lines%EP1,CO_lines%EP2,NCO_lev - 2,choice_ep)

    ! calculation of line intensities, substracting the cosmical nackground contribution
    DO j = 0, NCO_lev - 2
       jj = DBLE(j)
!       CO_lines(j)%Texc = (2._DP * CST1 * (jj + 1._DP)) /  &
!                          log((CO_lev(j+1)%weight * abs(CO_lev(j)%fracpop  )) / &
!                              (CO_lev(j)%weight   * abs(CO_lev(j+1)%fracpop)))
       CO_lines(j)%Texc = (2._DP * CST1 * (jj + 1._DP)) / (  &
                          log(CO_lev(j+1)%weight / CO_lev(j)%weight) + log(CO_lev(j)%fracpop / CO_lev(j+1)%fracpop))
       CO_lines(j)%BOLE = 1._DP / ((CO_lev(j)%fracpop   * CO_lev(j+1)%weight) / &
                                   (CO_lev(j+1)%fracpop * CO_lev(j)%weight  ) - 1._DP)
       CO_lines(j)%intensity  = 2._DP * CST1 * (jj + 1._DP) * (1._DP - EXP(- CO_lines(j)%TAU)) * &
                              (CO_lines(j)%BOLE - CO_lines(j)%BOLB)
       CO_lines(j)%intensity2 = CST4 / abs(gradv) * (CO_lines(j)%A_CO * CO_lev(j)%fracpop * Dens_CO * &
                              CO_lines(j)%EP2 / (2._DP * (jj + 1._DP) * Brot_CO) ** 2) * &
                              (1._DP - CO_lines(j)%BOLB/CO_lines(j)%BOLE)
    END DO


  END SUBROUTINE COMPUTE_CO

  
END MODULE MODULE_CO
