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

MODULE MODULE_READ_FE_DATA
   !-----------------------------------------------------------------------
   ! Subroutine called in INITIALIZE to read in data
   ! Data used in LINE_EXCIT, subroutine COLFEPL
   !
   ! Module to read in data required for subroutine FE_LEVEL_POPULATIONS
   ! Reads in atomic data: 
   !                       effective collision strengths (gamfepl)
   !                       radiative decay rates         (aijfepl)
   !                       upper level                   (iupfepl)
   !                       lower level                   (jlofepl)
   !----------------------------------------------------------------------
   USE MODULE_TOOLS, ONLY      : GET_FILE_NUMBER 

   IMPLICIT NONE
   INCLUDE "precision.f90"

   INTEGER (KIND=LONG),               PUBLIC  :: gam_n, aij_n
   INTEGER (KIND=LONG),               PRIVATE :: i, up, down
   REAL (KIND=DP),                    PRIVATE :: gam, A
   REAL (KIND=DP),DIMENSION(35,35),   PUBLIC  :: gamfepl
   REAL (KIND=DP),DIMENSION(256),     PUBLIC  :: aijfepl
   INTEGER, DIMENSION(256),           PUBLIC  :: iupfepl, jlofepl

   ! variables from "precision.f90" are private
   PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



   SUBROUTINE READ_FE_DATA

      ! Initialize
      gamfepl = 0.0_DP
      aijfepl = 0.0_DP
      iupfepl = 0
      jlofepl = 0

      ! Effective collisions strengths (Fe-electron collisions) for 1000 K
      gam_n=GET_FILE_NUMBER()
      OPEN (gam_n, file='input/gamma_eFe', status='old')        
      DO i = 1, 595
         READ (gam_n,*) down, up, gam
         gamfepl(down, up) = gam
      ENDDO
      CLOSE (gam_n)

      ! Radiative decay rates, upper and lower levels and transition energies
      aij_n=GET_FILE_NUMBER()           
      OPEN (aij_n, file='input/Aval_Fe', status='old') 
      DO i = 1, 256
         READ (aij_n,*) up, down, A
         aijfepl(i) = A
         iupfepl(i) = up
         jlofepl(i) = down
      ENDDO
      CLOSE (aij_n)

   END SUBROUTINE READ_FE_DATA



END MODULE MODULE_READ_FE_DATA
