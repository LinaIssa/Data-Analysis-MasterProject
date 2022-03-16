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

  !*****************************************************************************
  !** This file is bo be included in each module of the mhd shock code.       **
  !** It contains parameters used for the precision of real and integer data  **
  !** Do : real(KIND=DP) :: x                                                 **
  !**      integer(KIND=LONG) :: y                                            **
  !** The variable minus_infinity is the smallest negative number available   **
  !** The variable plus_infinity is the greatest positive number available    **
  !*****************************************************************************

  INTEGER, PARAMETER :: DP = selected_real_kind(P=15) ! real kind
  INTEGER, PARAMETER :: LONG=KIND(1)                ! integer kind
  REAL(KIND=DP),PARAMETER :: minus_infinity = -HUGE(1._DP) ! smallest < 0. number
  REAL(KIND=DP),PARAMETER :: plus_infinity  =  HUGE(1._DP) ! greatest > 0. number
