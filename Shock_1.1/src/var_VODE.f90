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

MODULE MODULE_VAR_VODE
  !*****************************************************************************
  !** The module 'MODULE_VAR_VODE' contains variables needed for the          **
  !** DVODE package                                                           **
  !** + some variables useful as control statement.                           **
  !**                                                                         **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! variables read in READ_PARAMETERS
  INTEGER                                   :: integ_type        ! type of integration over the mute variable z (lin, log, or adaptative)
  REAL(KIND=DP)                             :: duration_max      ! max. duration of the shock (years)
  REAL(KIND=DP)                             :: length_max        ! max. length of the shock (cm)
  REAL(KIND=DP)                             :: Eps_V             ! initial accuracy (argument of DVODE)

  ! other arguments of DVODE
  REAL(KIND=DP)                             :: T0_V   = 0.0_DP   ! initial value of the integration variable
  REAL(KIND=DP)                             :: H0_V   = 1.D05    ! initial step length
  REAL(KIND=DP)                             :: Tout_V = 1.D05    ! value of T at wich output should occur
  INTEGER(KIND=LONG)                        :: MF_V = 22         ! chosen numerical procedure
  INTEGER(KIND=LONG)                        :: Itask_V = 1       ! Which Task DVODE should do
  INTEGER(KIND=LONG)                        :: Istate_V = 1      ! How was it done (2 on normal return)
  INTEGER(KIND=LONG)                        :: MAXORD_V = 0      ! Max order of scheme
  INTEGER(KIND=LONG)                        :: MXSTEP_V = 100000 ! max number of steps
  INTEGER                                   :: liw_V, lrw_V      ! Size of work arrays
  INTEGER                                   :: itol_V = 4        ! Tolerance mode
  INTEGER                                   :: iopt_V = 1        ! Optional input ?
  REAL (kind=DP), ALLOCATABLE, DIMENSION(:) :: rtol_V
  REAL (kind=DP), ALLOCATABLE, DIMENSION(:) :: atol_V
  INTEGER,                     DIMENSION(1) :: ipar_V
  REAL (kind=DP),              DIMENSION(1) :: rpar_V
  REAL (kind=DP), ALLOCATABLE, DIMENSION(:) :: rwork_V
  INTEGER,        ALLOCATABLE, DIMENSION(:) :: iwork_V

  !------------------------------------------------------------

! REAL(KIND=DP) :: T, H, Hdone, Hnext
! REAL(KIND=DP) :: T, Hdone, Hnext
  REAL(KIND=DP) :: Hdone, Hnext

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_VAR_VODE

