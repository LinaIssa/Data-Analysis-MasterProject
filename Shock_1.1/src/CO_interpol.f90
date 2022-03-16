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

MODULE MODULE_CO_INTERPOL
  !*****************************************************************************
  !** The module 'INTERPOL_CO' contains variables and subroutines related to  **
  !** the interpolation of CO rate coefficients with Tn                       **
  !*****************************************************************************

  USE MODULE_PHYS_VAR, ONLY         : Tn, NCO_lev
  USE MODULE_CONSTANTS, ONLY        : pi, Zero, kB, h, c
  USE MODULE_CO
  USE MODULE_DEF
  USE MODULE_DEBUG_JLB
  USE MODULE_INITIALIZE

  IMPLICIT NONE
  INCLUDE "precision.f90"

!  !-----------------------------------------
!  ! parameters relative to CO
!  !-----------------------------------------
!  REAL(KIND=DP), PARAMETER :: Brot_CO   = 57635.96828_DP  ! rotational constant (MHz) NIST value
!  REAL(KIND=DP), PARAMETER :: Moment_CO = 0.11011_DP      ! dipolar moment (debye ; 1 debye = 1.e-18 CGS) NIST value
!  INTEGER                  :: err                         ! allocation  parameter
!
!  !---------------------------------------------
!  ! collision rates
!  ! read in READ_CO_RATES, used in EVOLUTION_CO
!  !---------------------------------------------
!
!  !--- CO-oH2 collision rates (Flower, 2001) ---
!  ! max value of J for which collision rates are available
!  INTEGER(KIND=LONG), PARAMETER                                        :: Jmin_oH2 = 0, Jmax_oH2 = 40
!  ! number of available temperatures for tabulated collision rates
!  INTEGER(KIND=LONG), PARAMETER                                        :: nt_oH2   = 14  
!  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
!  REAL(KIND=DP), DIMENSION(Jmin_oH2:Jmax_oH2,Jmin_oH2:Jmax_oH2,nt_oH2) :: Rate_oH2
!  ! temperatures which defines Rate_CO(i,j,k)
!  REAL(KIND=DP), DIMENSION(nt_oH2)                                     :: Temp_oH2 
!  ! new rates, calculated at Tn (useful collision rate array)
!  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                           :: CR_oH2
!
!  !--- CO-pH2 collision rates (Flower, 2001) ---
!  ! max value of J for which collision rates are available
!  INTEGER(KIND=LONG), PARAMETER                                        :: Jmin_pH2 = 0, Jmax_pH2 = 40 
!  ! number of available temperatures for tabulated collision rates
!  INTEGER(KIND=LONG), PARAMETER                                        :: nt_pH2   = 14  
!  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
!  REAL(KIND=DP), DIMENSION(Jmin_pH2:Jmax_pH2,Jmin_pH2:Jmax_pH2,nt_pH2) :: Rate_pH2
!  ! temperatures which defines Rate_CO(i,j,k)
!  REAL(KIND=DP), DIMENSION(nt_pH2)                                     :: temp_pH2 
!  ! new rates, calculated at Tn (useful collision rate array)
!  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                           :: CR_pH2
!
!  !--- CO-H collision rates (5 K < T < 3000 K ; Balakrishnan et al., ApJ, 2002) ---
!  ! max value of J for which collision rates are available
!  INTEGER(KIND=LONG), PARAMETER                              :: Jmin_H = 0, Jmax_H = 16
!  ! number of available temperatures for tabulated collision rates
!  INTEGER(KIND=LONG), PARAMETER                              :: nt_H   = 20  
!  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
!  REAL(KIND=DP), DIMENSION(Jmin_H:Jmax_H,Jmin_H:Jmax_H,nt_H) :: Rate_H
!  ! temperatures which defines Rate_CO(i,j,k)
!  REAL(KIND=DP), DIMENSION(nt_H)                             :: temp_H 
!  ! new rates, calculated at Tn (useful collision rate array)
!  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                 :: CR_H
!
!  !--- CO-He collision rates (Cecchi-Pestellini et al., ApJ, 2002) ---
!  ! max value of J for which collision rates are available
!  INTEGER(KIND=LONG), PARAMETER                                   :: Jmin_He = 0, Jmax_He = 14 
!  ! number of available temperatures for tabulated collision rates
!  INTEGER(KIND=LONG), PARAMETER                                   :: nt_He   = 10  
!  ! rates (cm3/s) : Rate_CO(i,j,k) -> i = initial state, j = final state, k = temperature
!  REAL(KIND=DP), DIMENSION(Jmin_He:Jmax_He,Jmin_He:Jmax_He,nt_He) :: Rate_He
!  ! temperatures which defines Rate_CO(i,j,k)
!  REAL(KIND=DP), DIMENSION(nt_He)                                 :: temp_He
!  ! new rates, calculated at Tn (useful collision rate array)
!  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE                      :: CR_He

CONTAINS

  SUBROUTINE INTERPOL

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
  
  CR_oH2 = Zero
  CR_pH2 = Zero
  CR_H   = Zero
  CR_He  = Zero

  ! calls the interpolation procedure for the rate coefficients
!  CALL INTERPOLATION(Jmin_oH2, Jmax_oH2, Rate_oH2, temp_oH2, CR_oH2, nt_oH2)
!  CALL INTERPOLATION(Jmin_pH2, Jmax_pH2, Rate_pH2, temp_pH2, CR_pH2, nt_pH2)
!  CALL INTERPOLATION(Jmin_H  , Jmax_H  , Rate_H  , temp_H  , CR_H  , nt_H  )
!  CALL INTERPOLATION(Jmin_He , Jmax_He , Rate_He , temp_He , CR_He , nt_He )

  CALL zob(Jmin_oH2, Jmax_oH2, Rate_oH2, temp_oH2, CR_oH2, nt_oH2)
  CALL zob(Jmin_pH2, Jmax_pH2, Rate_pH2, temp_pH2, CR_pH2, nt_pH2)
  CALL zob(Jmin_H  , Jmax_H  , Rate_H  , temp_H  , CR_H  , nt_H  )
  CALL zob(Jmin_He , Jmax_He , Rate_He , temp_He , CR_He , nt_He )

!  write(*,*) 'Rate_oH2(5,4,1)=', Rate_oH2(5,4,1)

  END SUBROUTINE INTERPOL

  SUBROUTINE zob(Jmin_rate,Jmax_rate,Rate_CO,Temp_CO,CR,nt)
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
!    USE MODULE_PHYS_VAR, ONLY : Tn

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

!    write(*,*) 'Tn_interpol=', Tn_interpol

    ! determination of the nearest temperatures (higher and lower) in the array of
    ! available collision rates that brackett the kinetic temperature
    ! indexes : ilo, ihi
    it = nt+1
    DO WHILE (ilo /= it)
       it = it-1
       IF (temp_CO(it) <= Tn_interpol) THEN
          ilo = it
       END IF
    END DO
    it = 0
    DO WHILE (ihi /= it)
       it = it+1
       IF (temp_CO(it) >= Tn_interpol) THEN
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
                DE = CST1 * (ji * (ji + 1) - jf * (jf + 1)) / Tn_interpol
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
                CR(jf,ji) = Rate_CO(ji,jf,ilo) + (Tn_interpol - temp_CO(ilo)) * (Rate_CO(ji,jf,ihi) - &
                            Rate_CO(ji,jf,ilo))/ (temp_CO(ihi) - temp_CO(ilo))
                DE = CST1 * (ji * (ji + 1) - jf * (jf + 1)) / Tn_interpol
                gl = 2._DP * jjf + 1._DP
                gu = 2._DP * jji + 1._DP
                R = gl /gu
                CR(ji,jf) = CR(jf,ji) * exp(- DE) / R
             END IF
          END DO
       END DO
    END IF
    

  END SUBROUTINE zob

  
END MODULE MODULE_CO_INTERPOL
