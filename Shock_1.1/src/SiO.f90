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

MODULE MODULE_SiO
  !*****************************************************************************
  !** The module 'MODULE_SiO' contains variables related to the SiO molecule  **
  !**    physical constants, Einstein coeff., collision rates                 **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !---------------------------
  ! physical constants
  !---------------------------
  INTEGER(KIND=LONG), PARAMETER :: Jmin_SiO=0, Jmax_SiO=20 ! min. and max. values of J
  REAL(KIND=DP), PARAMETER      :: Brot_SiO=21787.453_DP   ! rotational constant (MHz)
  REAL(KIND=DP), PARAMETER      :: Moment_SiO=3.0982_DP    ! dipolar moment (Debye)

  !--------------------------------------------
  ! Einstein coefficients
  ! calculated in EINSTEIN_COEFF_SiO
  !--------------------------------------------
  REAL(KIND=DP), DIMENSION(Jmax_SiO) :: A_SiO  ! spontaneous emission (CGS)
  REAL(KIND=DP), DIMENSION(Jmax_SiO) :: Bs_SiO ! stimulated emission (CGS)
  REAL(KIND=DP), DIMENSION(Jmax_SiO) :: Ba_SiO ! absorption (CGS)

  !-------------------------------------
  ! SiO-H2 collision rates
  ! variables read in READ_SiO_RATES
  !-------------------------------------
  ! rates (cm3/s) : Rate_SiO(i,j,k) -> i=initial state, j=final state, k=temperature
  REAL(KIND=DP),DIMENSION(Jmin_SiO:Jmax_SiO,Jmin_SiO:Jmax_SiO,8) :: Rate_SiO
  ! temperatures (K) which define Rate_SiO(i,j,k)
  REAL(KIND=DP), DIMENSION(8) :: temp_SiO

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

CONTAINS

  SUBROUTINE READ_SiO_RATES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads collision rates for SiO-H2.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    Rate_SiO, temp_SiO
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE

    INTEGER(KIND=LONG) :: i, Ji, Jf, error, Nrates
    CHARACTER(LEN=*), PARAMETER :: name_file_SiO='input/coeff_SiO.in'
    CHARACTER(LEN=*), PARAMETER :: format_SiO_temp='(8X,8F9.1)'
    INTEGER                     :: file_SiO

    ! initialization
    Rate_SiO=Zero
    temp_SiO=Zero

    ! file opening
    file_SiO = GET_FILE_NUMBER()
    OPEN(file_SiO,file=name_file_SiO,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! comments
    DO i=1,8
       READ(file_SiO,*)
    ENDDO
    ! temperatures
    READ(file_SiO,format_SiO_temp)temp_SiO
    ! comments
    DO i=1,3
       READ(file_SiO,*)
    ENDDO
    ! collision rates
    Nrates=SIZE(Rate_SiO,dim=1)*(SIZE(Rate_SiO,dim=2)-1)
    error=0
    DO i=1, Nrates
       READ(file_SiO,*,iostat=error)Ji,Jf,Rate_SiO(Ji,Jf,:)
       IF (error>0) STOP "*** WARNING, error in READ_SIO_RATES"
    ENDDO

  END SUBROUTINE READ_SiO_RATES


  SUBROUTINE EINSTEIN_COEFF_SiO
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    computes Einstein coefficients for the SiO molecule
    !    method=simple rotator molecule
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    A_SiO, Bs_SiO and Ba_SiO
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: j
    REAL(KIND=DP) :: jj, Gl, Gu

    ! initializations
    A_SiO=Zero ; Bs_SiO=Zero ; Ba_SiO=Zero

    ! computation
    DO j=1,Jmax_SiO
       jj=DBLE(j) ! conversion to real
       Gl=2._DP*(jj-1._DP)+1._DP
       Gu=2._DP*jj+1._DP

       Bs_SiO(j)=7.896D8 * moment_SiO**2 * jj/Gu ! stimulated emission
       Ba_SiO(j)=Bs_SiO(j) * Gu/Gl               ! absorption
       A_SiO(j)=11.7944D-29 * Brot_SiO**3 *jj**3 * Bs_SiO(j) ! spontaneous emission
    END DO
  END SUBROUTINE EINSTEIN_COEFF_SiO

END MODULE MODULE_SiO
