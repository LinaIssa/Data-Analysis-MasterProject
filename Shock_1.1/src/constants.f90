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

MODULE MODULE_CONSTANTS
  !*****************************************************************************
  !** The module 'MODULE_CONSTANTS' contains mathematical and physical        **
  !** constants, and unit conversion factors                                  **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !--- mathematical constants ---
  REAL(KIND=DP), PARAMETER :: pi   = 3.141592653589793238_dp
  REAL(KIND=DP), PARAMETER :: Zero = 0.0_dp
  REAL(KIND=dp), PARAMETER :: sqpi = 1.77245385090551602_dp  ! sqrt(pi)

  !--- physical constants ---
  REAL(KIND=DP), PARAMETER :: mP      = 1.672621777e-24_dp      ! mass of the proton (g)
  REAL(KIND=DP), PARAMETER :: me      = 9.10938291e-28_dp       ! mass of the electron (g)
  REAL(KIND=DP), PARAMETER :: qe      = 4.80320450571347e-10_dp ! charge of the electron (esu)
  REAL(KIND=dp), PARAMETER :: qe2     = 2.30707735237062e-19_dp ! xqe * xqe

  REAL(KIND=DP), PARAMETER :: alpha_H  = 6.67e-25_dp  ! polarisability of H (cm-3)
  REAL(KIND=DP), PARAMETER :: alpha_H2 = 7.70e-25_dp  ! polarisability of H2 (cm-3)
  REAL(KIND=DP), PARAMETER :: alpha_He = 2.10e-25_dp  ! polarisability of He (cm-3)

  REAL(KIND=DP), PARAMETER :: bohr = 5.29177e-9_dp     ! bohr radius
  REAL(KIND=DP), PARAMETER :: kB   = 1.380658e-16_dp   ! Boltzmann's constant (erg.K-1)
  REAL (KIND=dp),PARAMETER :: hp   = 6.62606957e-27_dp ! Planck constant h (erg s)
  REAL(KIND=DP), PARAMETER :: R    = 8.314510e7_dp     ! molar gas constant (erg.K-1.mol-1)
  REAL(KIND=dp), PARAMETER :: clum = 2.99792458e10_dp  ! Speed of light c (cm s-1)

  REAL(KIND=DP), PARAMETER :: wien = 2.8977729e-3_dp   ! Wien constant

  !--- units conversion ---
  REAL(KIND=DP), PARAMETER :: EVerg   = 1.60218e-12_dp! 1 eV       =  1.60218D-12 erg
  REAL(KIND=DP), PARAMETER :: amu     = 1.66054e-24_dp! 1 amu      =  1.66054D-24 g
  REAL(KIND=DP), PARAMETER :: YEARsec = 3.15569e7_dp  ! 1 year     =  3.15569D7 s
  REAL(KIND=DP), PARAMETER :: kCaleV  = 4.3363e-2_dp  ! 1 kCal/mol =  4.3363D-2 eV
  real(KIND=DP), parameter :: AUcgs   = 6.127e-9_dp   ! conversion between au and cgs
  REAL(KIND=DP), PARAMETER :: parsec  = 3.0857e18_dp  ! 1 pc = 3.0857D18 cm
  REAL(KIND=DP), PARAMETER :: arc_ster= 1/206265e0_dp ! conversion from arcsec to steradians 
  REAL(KIND=dp), PARAMETER :: hcsurk  = 1.43877695998382e+00_dp ! From Wave number (cm-1) to Temperature (Kelvin)

  !--- standard output : UNIX -> 6, MAC -> 9 ---
  INTEGER(KIND=LONG), PARAMETER :: screen = 6
  INTEGER(KIND=LONG), PUBLIC    :: iwrtmp = 13

  !--- numerical limites ---
  REAL(KIND=dp), PARAMETER  :: l_hu    = 250.0_dp

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_CONSTANTS
