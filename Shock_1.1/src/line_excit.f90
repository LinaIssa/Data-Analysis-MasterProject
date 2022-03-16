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

 MODULE MODULE_LINE_EXCIT
   !***************************************************************************
   ! Should replace part or all MODULE_FINE_STRUCTURE
   ! Sylvie Cabrit & Jacques Le Bourlot - Sept 2002
   ! Shorter version SC & GPdF - Nov 2002
   !***************************************************************************

   IMPLICIT NONE
   INCLUDE "precision.f90"

   ! ==========================================================================
   ! Atomic or molecular parameter:
   ! ==========================================================================
   !
   !  All variables names use the same scheme
   !  first 3 caracters : which parameter
   !                nlv : number of levels
   !                gst : statistical weight (usually, g = 2 J + 1)
   !                elk : energy in Kelvin

   !                ntr : number of radiative transitions
   !                aij : Aij Einstein coefficien (s-1)
   !                wlk : transition energy (in K)
   !                emi : emissivity of transition (erg . cm-3 s-1)
   !                iup : index of upper level
   !                jlo : index of lower level

   !  next 3 caracters   : which specy
   !                cat  : atomic carbon    (C I)
   !                nat  : atomic nitrogen  (N I)
   !                oat  : atomic oxygen    (O I)
   !                sat  : atomic sulfur    (S I)
   !                siat : atomic silicium  (Si I)
   !                cpl  : ionised carbon   (C II)
   !                npl  : ionised nitrogen (N II)
   !                opl  : ionised oxygen   (O II)
   !                spl  : ionised sulfur   (S II)
   !                sipl : ionised silicium (Si II)


   ! ----------------------------------------------
   ! H
   ! ----------------------------------------------
   ! Data from ??? (anderson 2000?) NIST

   ! Energy levels.
   !        1 - 1s J = 1/2
   !        2 - 2p J = 1/2
   !        3 - 2s J = 1/2
   !        4 - 2p J = 3/2
   !        5 - 3p J = 1/2
   !        6 - 3s J = 1/2
   !        7 - 3d J = 3/2
   !        8 - 3p J = 3/2
   !        9 - 3d J = 5/2

   INTEGER,        PUBLIC, PARAMETER          :: nlvhat = 9
   REAL (kind=DP), PUBLIC, DIMENSION (nlvhat) :: gsthat
   DATA gsthat                                 / 2.0_DP, 2.0_DP, 2.0_DP, 4.0_DP, 2.0_DP, 2.0_DP, & 
                                                 4.0_DP, 4.0_DP, 6.0_DP/
   REAL (kind=DP), PUBLIC, DIMENSION (nlvhat) :: elkhat
   DATA elkhat                                 / 0.0_DP, 118352.237573_DP, 118352.288342_DP, & 
                                                 118352.764005_DP, 140269.547252_DP, 140269.562361_DP, &
                                                 140269.702976_DP, 140269.703232_DP, 140269.754963_DP/
   REAL (kind=DP), PUBLIC, DIMENSION (nlvhat) :: pop_hat   ! Atomic Hydrogen populations

   ! Transitions
   !       1: 2-1 (2p J=1/2 -> 1s J=1/2) -> 1215.673645 Angstrom, Lyman alpha)
   !       2: 3-1 (2s J=1/2 -> 1s J=1/2) -> 2 photon continuum at lambda_max = 1400.000000 Angstrom
   !       3: 4-1 (2p J=3/2 -> 1s J=1/2) -> 1215.668237 Angstrom (Lyman alpha)
   !       4: 5-1 (3p J=1/2 -> 1s J=1/2) -> 1025.722966 Angstrom (Lyman beta)
   !       5: 5-3 (3p J=1/2 -> 2s J=1/2) -> 6562.771533 Angstrom
   !       6: 6-1 (3s J=1/2 -> 1s J=1/2) -> 1025.722855 Angstrom (Lyman beta)
   !       7: 6-2 (3s J=1/2 -> 2p J=1/2) -> 6562.751807 Angstrom
   !       8: 6-4 (3s J=1/2 -> 2p J=3/2) -> 6562.909442 Angstrom
   !       9: 7-1 (3d J=3/2 -> 1s J=1/2) -> 1025.721827 Angstrom (Lyman beta)
   !       10: 7-2 (3d J=3/2 -> 2p J=1/2) -> 6562.709702 Angstrom
   !       11: 7-4 (3d J=3/2 -> 2p J=3/2) -> 6562.867336 Angstrom
   !       12: 8-1 (3p J=3/2 -> 1s J=1/2) -> 1025.721825 Angstrom (Lyman beta)
   !       13: 8-3 (3p J=3/2 -> 2s J=1/2) -> 6562.724827 Angstrom
   !       14: 9-1 (3d J=5/2 -> 1s J=1/2) -> 1025.721447 Angstrom (Lyman beta)
   !       15: 9-4 (3d J=5/2 -> 2p J=3/2) -> 6562.851769 Angstrom


   INTEGER,        PUBLIC, PARAMETER          :: ntrhat = 15
   REAL (kind=DP), PUBLIC, DIMENSION (ntrhat) :: wlkhat
   DATA wlkhat                                 / 118352.23757339_DP, 102950.00000000_DP, 118352.23757339_DP, &
                                                 140269.54725244_DP, 21917.25891039_DP, 140269.56236103_DP, &
                                                 21917.32478764_DP, 21916.79835628_DP, 140269.70297558_DP, &
                                                 21917.46540219_DP, 21916.93897083_DP, 140269.70323169_DP, &
                                                 21917.41488964_DP, 140269.75496291_DP, 21916.99095815_DP/
   REAL (kind=DP), PUBLIC, DIMENSION (ntrhat) :: aijhat
   DATA aijhat                                 / 6.2649000E+08_DP, 8.2000000E-00_DP, 6.2649000E+08_DP, &
                                                 1.6725000E+08_DP, 2.2449000E+07_DP, 1.1090000E-06_DP, &
                                                 2.1046000E+06_DP, 4.2097000E+06_DP, 5.9380000E+02_DP, &
                                                 5.3877000E+07_DP, 1.0775000E+07_DP, 1.6725000E+08_DP, &
                                                 2.2448000E+07_DP, 5.9370000E+02_DP, 6.4651000E+07_DP/
   REAL (kind=DP), PUBLIC, DIMENSION (ntrhat) :: emihat, emihat_o, inthat
   INTEGER,        PUBLIC, DIMENSION (ntrhat) :: iuphat
   DATA iuphat                                 / 2, 3, 4, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 9/
   INTEGER,        PUBLIC, DIMENSION (ntrhat) :: jlohat
   DATA jlohat                                 / 1, 1, 1, 1, 3, 1, 2, 4, 1, 2, 4, 1, 3, 1, 4/
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrhat) :: namtrhat
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrhat) :: idtrhat


   ! ----------------------------------------------
   ! C
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 3P J = 0
   !        2 - 3P J = 1
   !        3 - 3P J = 2
   !        4 - 1D J = 2
   !        5 - 1S J = 0

   INTEGER,        PUBLIC, PARAMETER          :: nlvcat = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvcat) :: gstcat
   DATA gstcat                                 / 1.0_DP, 3.0_DP, 5.0_DP, 5.0_DP, 1.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvcat) :: elkcat
   DATA elkcat                                 / 0.0_DP, 23.5969_DP, 62.4454_DP, 14665.507_DP, 31147.902_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvcat) :: pop_cat   ! Atomic Carbon populations

   ! Transitions
   !       1: 2-1 -> 609.75 micron
   !       2: 3-2 -> 370.37 micron
   !       3: 3-1 -> 230.41 micron
   !       4: 4-3 -> 9850.26 Angstrom
   !       5: 4-2 -> 9824.13 Angstrom
   !       6: 4-1 -> 9808.32 Angstrom
   !       7: 5-4 -> 8727.13 Angstrom
   !       8: 5-3 -> 4627.346 Angstrom
   !       9: 5-2 -> 4621.570 Angstrom
   !      10: 5-1 

   INTEGER,        PUBLIC, PARAMETER          :: ntrcat = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntrcat) :: wlkcat
   DATA wlkcat                                 / 23.5969_DP, 38.8485_DP, 62.4454_DP, &
                                                14603.0616_DP, 14641.9101_DP, 14665.507_DP, &
                                                16482.395_DP, 31085.4566_DP, 31124.3051_DP, 31147.902_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrcat) :: aijcat
   DATA aijcat                                 / 7.93e-8_DP, 2.65e-7_DP, 1.71e-14_DP, &
                                                 2.44e-4_DP, 8.21e-5_DP, 7.77e-8_DP, &
                                                 5.28e-1_DP, 2.00e-5_DP, 2.71e-3_DP, 0.0e0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntrcat) :: emicat, emicat_o, intcat
   INTEGER,        PUBLIC, DIMENSION (ntrcat) :: iupcat
   DATA iupcat                                 / 2, 3, 3, 4, 4, 4, 5, 5, 5, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntrcat) :: jlocat
   DATA jlocat                                 / 1, 2, 1, 3, 2, 1, 4, 3, 2, 1 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrcat) :: namtrcat
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrcat) :: idtrcat


   ! ----------------------------------------------
   ! N
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)
   
   ! Energy levels.
   !        1 - 4So J = 3/2
   !        2 - 2Do J = 5/2
   !        3 - 2Do J = 3/2
   !        4 - 2Po J = 1/2
   !        5 - 2Po J = 3/2

   INTEGER,        PUBLIC, PARAMETER          :: nlvnat = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvnat) :: gstnat
   DATA gstnat                                 / 4.0_DP, 6.0_DP, 4.0_DP, 2.0_DP, 4.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvnat) :: elknat
   DATA elknat                                 / 0.0_DP, 27660.82_DP, 27673.36_DP, 41494.431_DP, 41494.986_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvnat) :: pop_nat   ! Atomic Nitrogen populations

   ! Transitions
   !       1: 3-2 -> 1.147 mm (261.21 GHz)
   !       2: 4-3 -> 1.04104 micron
   !       3: 5-3 -> 1.04100 micron
   !       4: 4-2 -> 1.04010 micron
   !       5: 5-2 -> 1.04006 micron
   !       6: 2-1 -> 5200.257 Angstrom
   !       7: 3-1 -> 5197.902 Angstrom
   !       8: 4-1 -> 3466.543 Angstrom
   !       9: 5-1 -> 3466.497 Angstrom
   !      10: 5-4

   INTEGER,        PUBLIC, PARAMETER          :: ntrnat = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntrnat) :: wlknat
   DATA wlknat                                 / 12.54_DP, 13821.071_DP, 13821.626_DP, &
                                                13833.611_DP, 13834.166_DP, 27660.82_DP, &
                                                27673.36_DP, 41494.431_DP, 41494.986_DP, 0.555_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrnat) :: aijnat
   DATA aijnat                                 / 1.27e-8_DP, 5.29e-2_DP, 2.76e-2_DP, &
                                                 3.45e-2_DP, 6.14e-2_DP, 7.27e-6_DP, &
                                                 2.02e-5_DP, 2.71e-3_DP, 6.58e-3_DP, 0.0e0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntrnat) :: eminat, eminat_o, intnat
   INTEGER,        PUBLIC, DIMENSION (ntrnat) :: iupnat
   DATA iupnat                                 / 3, 4, 5, 4, 5, 2, 3, 4, 5, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntrnat) :: jlonat
   DATA jlonat                                 / 2, 3, 3, 2, 2, 1, 1, 1, 1, 4 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrnat) :: namtrnat
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrnat) :: idtrnat


   ! ----------------------------------------------
   ! O
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 3P J = 2
   !        2 - 3P J = 1
   !        3 - 3P J = 0
   !        4 - 1D J = 2
   !        5 - 1S J = 0

   INTEGER,        PUBLIC, PARAMETER          :: nlvoat = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvoat) :: gstoat
   DATA gstoat                                 / 5.0_DP, 3.0_DP, 1.0_DP, 5.0_DP, 1.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvoat) :: elkoat
   DATA elkoat                                 / 0.0_DP, 227.717_DP, 326.582_DP, 22831.226_DP, 48621.932_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvoat) :: pop_oat   ! Atomic Oxygen populations

   ! Transitions
   !       1: 3-2 -> 145.53 micron
   !       2: 2-1 -> 63.19  micron
   !       3: 3-1 -> 44.06  micron
   !       4: 4-3 -> 6391.73 Angstrom
   !       5: 4-2 -> 6363.78 Angstrom
   !       6: 4-1 -> 6300.30 Angstrom
   !       7: 5-4 -> 5577.34 Angstrom
   !       8: 5-2 -> 2972.29 Angstrom
   !       9: 5-1 -> 2958.37 Angstrom
   !      10: 5-3

   INTEGER,        PUBLIC, PARAMETER          :: ntroat = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntroat) :: wlkoat
   DATA wlkoat                                 / 98.865_DP, 227.717_DP, 326.582_DP, &
                                                22504.644_DP, 22603.509_DP, 22831.226_DP, &
                                                25790.706_DP, 48394.215_DP, 48621.932_DP, 48295.35_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntroat) :: aijoat
   DATA aijoat                                 / 1.66e-5_dp, 8.46e-5_dp, 1.10e-10_dp, &
                                                 1.20e-6_dp, 2.17e-3_dp, 6.71e-3_dp, &
                                                 1.02e-0_dp, 7.83e-2_dp, 2.87e-4_dp, 0.0e0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntroat) :: emioat, emioat_o, intoat
   INTEGER,        PUBLIC, DIMENSION (ntroat) :: iupoat
   DATA iupoat                                 / 3, 2, 3, 4, 4, 4, 5, 5, 5, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntroat) :: jlooat
   DATA jlooat                                 / 2, 1, 1, 3, 2, 1, 4, 2, 1, 3 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntroat) :: namtroat
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntroat) :: idtroat


   ! ----------------------------------------------
   ! S
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 3P J = 2
   !        2 - 3P J = 1
   !        3 - 3P J = 0
   !        4 - 1D J = 2
   !        5 - 1S J = 0

   INTEGER,        PUBLIC, PARAMETER          :: nlvsat = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsat) :: gstsat
   DATA gstsat                                 / 5.0_DP, 3.0_DP, 1.0_DP, 5.0_DP, 1.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsat) :: elksat
   DATA elksat                                 / 0.0_DP, 569.86_DP, 825.37_DP, 13292.83_DP, 31913.28_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsat)        :: pop_sat   ! Atomic Sulfur populations

   ! Transitions
   !       1: 3-2 -> 56.31 micron
   !       2: 2-1 -> 25.25 micron
   !       3: 3-1 -> 17.43 micron
   !       4: 4-3 -> 1.15407 micron
   !       5: 4-2 -> 1.13089 micron
   !       6: 4-1 -> 1.08241 micron
   !       7: 5-4 -> 7725.05 Angstrom
   !       8: 5-2 -> 4589.26 Angstrom
   !       9: 5-1 -> 4507.31 Angstrom
   !      10: 5-3 

   INTEGER,        PUBLIC, PARAMETER          :: ntrsat = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsat) :: wlksat
   DATA wlksat                                 / 255.51_DP, 569.86_DP, 825.37_DP, &
                                                 12467.46_DP, 12722.97_DP, 13292.83_DP, &
                                                 18620.45_DP, 31343.42_DP, 31913.28_DP, 31087.91_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsat) :: aijsat
   DATA aijsat                                 / 3.02e-4_DP, 1.39e-3_DP, 6.71e-8_DP, &
                                                 3.84e-6_DP, 8.16e-3_DP, 2.78e-2_DP, &
                                                 1.53e0_DP, 3.50e-1_DP, 8.23e-3_DP, 0.0e0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsat) :: emisat, emisat_o, intsat
   INTEGER,        PUBLIC, DIMENSION (ntrsat) :: iupsat
   DATA iupsat                                 / 3, 2, 3, 4, 4, 4, 5, 5, 5, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntrsat) :: jlosat
   DATA jlosat                                 / 2, 1, 1, 3, 2, 1, 4, 2, 1, 3 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrsat) :: namtrsat
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrsat) :: idtrsat


   ! ----------------------------------------------
   ! Si
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 3P J = 0
   !        2 - 3P J = 1
   !        3 - 3P J = 2
   !        4 - 1D J = 2
   !        5 - 1S J = 0

   INTEGER,        PUBLIC, PARAMETER           :: nlvsiat = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsiat) :: gstsiat
   DATA gstsiat                                 / 1.0_DP, 3.0_DP, 5.0_DP, 5.0_DP, 1.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsiat) :: elksiat
   DATA elksiat                                 / 0.0_DP, 110.951_DP, 321.086_DP, 9062.998_DP, 22149.939_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsiat) :: pop_siat   ! Atomic Silicium populations

   ! Transitions
   !       1: 2-1 -> 129.682 micron
   !       2: 3-2 -> 68.472 micron
   !       3: 3-1 -> 44.812 micron
   !       4: 4-3 -> 1.6459 micron
   !       5: 4-2 -> 1.6073 micron
   !       6: 4-1 -> 1.5876 micron
   !       7: 5-4 -> 1.0994 micron
   !       8: 5-3 -> 6589.61 Angstrom
   !       9: 5-2 -> 6526.78 Angstrom
   !      10: 5-1

   INTEGER,        PUBLIC, PARAMETER           :: ntrsiat = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsiat) :: wlksiat
   DATA wlksiat                                 / 110.951_DP, 210.135_DP, 321.086_DP, &
                                                  8741.912_DP, 8952.047_DP, 9062.998_DP, &
                                                  13086.941_DP, 21828.853_DP, 22038.988_DP, 21828.853_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsiat) :: aijsiat
   DATA aijsiat                                 / 8.25e-6_DP, 4.21e-5_DP, 3.56e-13_DP, &
                                                  2.25e-3_DP, 7.93e-4_DP, 4.70e-7_DP, &
                                                  1.14e-0_DP, 9.02e-4_DP, 3.13e-2_DP, 0.0e0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsiat) :: emisiat, emisiat_o, intsiat
   INTEGER,        PUBLIC, DIMENSION (ntrsiat) :: iupsiat
   DATA iupsiat                                 / 2, 3, 3, 4, 4, 4, 5, 5, 5, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntrsiat) :: jlosiat
   DATA jlosiat                                 / 1, 2, 1, 3, 2, 1, 4, 3, 2, 1 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrsiat) :: namtrsiat
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrsiat) :: idtrsiat


   ! ----------------------------------------------
   ! C+
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 2P J = 1/2
   !        2 - 2P J = 3/2
   !        3 - 4P J = 1/2
   !        4 - 4P J = 3/2
   !        5 - 4P J = 5/2

   INTEGER,        PUBLIC, PARAMETER          :: nlvcpl = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvcpl) :: gstcpl
   DATA gstcpl                                 / 2.0_DP, 4.0_DP, 2.0_DP, 4.0_DP, 6.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvcpl) :: elkcpl
   DATA elkcpl                                 / 0.0_DP, 91.25_DP, 61874.63_DP, 61906.28_DP, 61947.00_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvcpl) :: pop_cpl   ! Ionised Carbon populations

   ! Transitions
   !       1: 2-1 -> 157.68 micron
   !       2: 3-2 -> 2328.12 Angstrom
   !       3: 4-2 -> 2326.93 Angstrom
   !       4: 5-2 -> 2325.40 Angstrom
   !       5: 3-1 -> 2324.69 Angstrom
   !       6: 4-1 -> 2323.50 Angstrom
   !       7: 5-1 ->  

   INTEGER,        PUBLIC, PARAMETER          :: ntrcpl = 7
   REAL (kind=DP), PUBLIC, DIMENSION (ntrcpl) :: wlkcpl
   DATA wlkcpl                                 / 91.25_DP, 61783.38_DP, 61815.03_DP, &
                                                 61855.75_DP, 61874.63_DP, 61906.28_DP, 61947.00_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrcpl) :: aijcpl
   DATA aijcpl                                 / 2.29e-6_DP, 6.55e+1_DP, 5.24e+0_DP, &
                                                 4.32e+1_DP, 5.53e+1_DP, 1.71e-0_DP, 0.0e0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntrcpl) :: emicpl, emicpl_o, intcpl
   INTEGER,        PUBLIC, DIMENSION (ntrcpl) :: iupcpl
   DATA iupcpl                                 / 2, 3, 4, 5, 3, 4, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntrcpl) :: jlocpl
   DATA jlocpl                                 / 1, 2, 2, 2, 1, 1, 1 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrcpl) :: namtrcpl
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrcpl) :: idtrcpl


   ! ----------------------------------------------
   ! N+
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 3P J = 0
   !        2 - 3P J = 1
   !        3 - 3P J = 2
   !        4 - 1D J = 2
   !        5 - 1S J = 0

   INTEGER,        PUBLIC, PARAMETER          :: nlvnpl = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvnpl) :: gstnpl
   DATA gstnpl                                 / 1.0_DP, 3.0_DP, 5.0_DP, 5.0_DP, 1.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvnpl) :: elknpl
   DATA elknpl                                 / 0.0_DP, 70.07_DP, 188.20_DP, 22037.48_DP, 47033.77_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvnpl) :: pop_npl   ! Ionised Nitrogen populations

   ! Transitions
   !       1: 2-1 -> 205.34 micron
   !       2: 3-2 -> 121.80 micron
   !       3: 3-1 -> 76.45 micron
   !       4: 4-3 -> 6583.45 Angstrom
   !       5: 4-2 -> 6548.05 Angstrom
   !       6: 4-1 -> 6527.23 Angstrom
   !       7: 5-4 -> 5754.59 Angstrom
   !       8: 5-3 -> 3070.55 Angstrom
   !       9: 5-2 -> 3062.83 Angstrom
   !      10: 5-1 -> 

   INTEGER,        PUBLIC, PARAMETER          :: ntrnpl = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntrnpl) :: wlknpl
   DATA wlknpl                                 / 70.07_DP, 118.13_DP, 188.20_DP, &
                                                 21849.28_DP, 21967.41_DP, 22037.48_DP, &
                                                 24996.29_DP, 46845.57_DP, 46963.7_DP, 47033.77_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrnpl) :: aijnpl
   DATA aijnpl                                 / 2.08e-6_DP, 7.46e-6_DP, 1.16e-12_DP, &
                                                 2.99e-3_DP, 1.01e-3_DP, 5.35e-7_DP, &
                                                 1.12e-0_DP, 1.51e-4_DP, 3.38e-2_DP, 0.0e0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntrnpl) :: eminpl, eminpl_o, intnpl
   INTEGER,        PUBLIC, DIMENSION (ntrnpl) :: iupnpl
   DATA iupnpl                                 / 2, 3, 3, 4, 4, 4, 5, 5, 5, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntrnpl) :: jlonpl
   DATA jlonpl                                 / 1, 2, 1, 3, 2, 1, 4, 3, 2, 1 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrnpl) :: namtrnpl
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrnpl) :: idtrnpl


   ! ----------------------------------------------
   ! O+
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 4So J = 3/2
   !        2 - 2Do J = 5/2
   !        3 - 2Do J = 3/2
   !        4 - 2Po J = 3/2
   !        5 - 2Po J = 1/2

   INTEGER,        PUBLIC, PARAMETER          :: nlvopl = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvopl) :: gstopl
   DATA gstopl                                 / 4.0_DP, 6.0_DP, 4.0_DP, 4.0_DP, 2.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvopl) :: elkopl
   DATA elkopl                                 / 0.0_DP, 38575.94_DP, 38604.75_DP, 58226.77_DP, 58229.63_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvopl) :: pop_opl   ! Ionised Oxygen populations

   ! Transitions
   !       1: 5-4 -> 5.03 mm (59.66 GHz)
   !       2: 3-2 -> 0.4995 mm (600.18 GHz)
   !       3: 4-3 -> 7330.73 Angstrom
   !       4: 5-3 -> 7329.67 Angstrom
   !       5: 4-2 -> 7319.99 Angstrom
   !       6: 5-2 -> 7318.92 Angstrom
   !       7: 2-1 -> 3728.815 Angstrom
   !       8: 3-1 -> 3726.032 Angstrom
   !       9: 4-1 -> 2470.341 Angstrom
   !      10: 5-1 -> 2470.219 Angstrom

   INTEGER,        PUBLIC, PARAMETER          :: ntropl = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntropl) :: wlkopl
   DATA wlkopl                                 / 2.86_DP, 28.81_DP, 19622.02_DP, &
                                                 19624.88_DP, 19650.83_DP, 19653.69_DP, &
                                                 38575.94_DP, 38604.75_DP, 58226.77_DP, &
                                                 58229.63_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntropl) :: aijopl
   DATA aijopl                                 / 2.08e-11_DP, 1.20e-7_DP, 6.14e-2_DP, &
                                                 1.02e-1_DP, 1.17e-1_DP, 6.15e-2_DP, &
                                                 3.82e-5_DP, 1.65e-4_DP, 5.64e-2_DP, &
                                                 2.32e-2_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntropl) :: emiopl, emiopl_o, intopl
   INTEGER,        PUBLIC, DIMENSION (ntropl) :: iupopl
   DATA iupopl                                 / 5, 3, 4, 5, 4, 5, 2, 3, 4, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntropl) :: jloopl
   DATA jloopl                                 / 4, 2, 3, 3, 2, 2, 1, 1, 1, 1 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntropl) :: namtropl
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntropl) :: idtropl


   ! ----------------------------------------------
   ! S+
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 4So J = 3/2
   !        2 - 2Do J = 3/2
   !        3 - 2Do J = 5/2
   !        4 - 2Po J = 1/2
   !        5 - 2Po J = 3/2

   INTEGER,        PUBLIC, PARAMETER          :: nlvspl = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvspl) :: gstspl
   DATA gstspl                                 / 4.0_DP, 4.0_DP, 6.0_DP, 2.0_DP, 4.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvspl) :: elkspl
   DATA elkspl                                 / 0.0_DP, 21370.92_DP, 21416.66_DP, 35287.17_DP, 35354.38_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvspl)        :: pop_spl   ! Ionised Sulfur populations

   ! Transitions
   !       1: 3-2 -> 0.31456 mm (953.04 GHz)
   !       2: 5-4 -> 0.21409 mm (1400.33 GHz)
   !       3: 4-3 -> 1.0373 micron
   !       4: 4-2 -> 1.0339 micron
   !       5: 5-3 -> 1.0323 micron
   !       6: 5-2 -> 1.0290 micron
   !       7: 2-1 -> 6730.82 Angstrom
   !       8: 3-1 -> 6716.44 Angstrom
   !       9: 4-1 -> 4076.35 Angstrom
   !      10: 5-1 -> 4068.60 Angstrom

   INTEGER,        PUBLIC, PARAMETER          :: ntrspl = 10
   REAL (kind=DP), PUBLIC, DIMENSION (ntrspl) :: wlkspl
   DATA wlkspl                                 / 45.74_DP, 67.21_DP, 13870.51_DP, &
                                                 13916.25_DP, 13937.72_DP, 13983.46_DP, &
                                                 21370.92_DP, 21416.66_DP, 35287.17_DP, &
                                                 35354.38_DP /
 
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrspl) :: aijspl
   DATA aijspl                                 / 3.35e-7_DP, 1.03e-6_DP, 7.79e-2, &
                                                 1.63e-1_DP, 1.79e-1_DP, 1.33e-1_DP, &
                                                 8.82e-4_DP, 2.60e-4_DP, 9.06e-2_DP, &
                                                 2.25e-1_DP /
 
   REAL (kind=DP), PUBLIC, DIMENSION (ntrspl) :: emispl, emispl_o, intspl
   INTEGER,        PUBLIC, DIMENSION (ntrspl) :: iupspl
   DATA iupspl                                 / 3, 5, 4, 4, 5, 5, 2, 3, 4, 5 /
   INTEGER,        PUBLIC, DIMENSION (ntrspl) :: jlospl
   DATA jlospl                                 / 2, 4, 3, 2, 3, 2, 1, 1, 1, 1 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrspl) :: namtrspl
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrspl) :: idtrspl


   ! ----------------------------------------------
   ! Si+
   ! ----------------------------------------------
   ! Data from NIST (12 Sept 2002)

   ! Energy levels.
   !        1 - 2P J = 1/2
   !        2 - 2P J = 3/2
   !        3 - 4P J = 1/2
   !        4 - 4P J = 3/2
   !        5 - 4P J = 5/2

   INTEGER,        PUBLIC, PARAMETER           :: nlvsipl = 5
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsipl) :: gstsipl
   DATA gstsipl                                 / 2.0_DP, 4.0_DP, 2.0_DP, 4.0_DP, 6.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsipl) :: elksipl
   DATA elksipl                                 / 0.0_DP, 413.29_DP, 61617.06_DP, 61772.93_DP, 62025.14_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvsipl) :: pop_sipl   ! Ionised Silicium populations

   ! Transitions
   !       1: 2-1 -> 34.814 micron
   !       2: 3-2 -> 2350.17 Angstrom
   !       3: 4-2 -> 2344.20 Angstrom
   !       4: 5-2 -> 2334.60 Angstrom
   !       5: 3-1 -> 2334.40 Angstrom
   !       6: 4-1 -> 2328.52 Angstrom

   INTEGER,        PUBLIC, PARAMETER           :: ntrsipl = 6
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsipl) :: wlksipl
   DATA wlksipl                                 / 413.29_DP, 61203.77_DP, 61359.64_DP, &
                                                  61611.85_DP, 61617.06_DP, 61772.93_DP /
   ! Revision DRF - Nov 2002
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsipl) :: aijsipl
   DATA aijsipl                                 / 2.2e-4_DP, 4.9e+3_DP, 1.7e+3_DP, &
                                                 2.7e+3_DP, 6.3e+3_DP, 2.0e+1_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (ntrsipl) :: emisipl, emisipl_o, intsipl
   INTEGER,        PUBLIC, DIMENSION (ntrsipl) :: iupsipl
   DATA iupsipl                                 / 2, 3, 4, 5, 3, 4 /
   INTEGER,        PUBLIC, DIMENSION (ntrsipl) :: jlosipl
   DATA jlosipl                                 / 1, 2, 2, 2, 1, 1 /
   CHARACTER(LEN=6), PUBLIC, DIMENSION (ntrsipl) :: namtrsipl
   CHARACTER(LEN=3), PUBLIC, DIMENSION (ntrsipl) :: idtrsipl


   ! ----------------------------------------------
   ! Fe+
   ! ----------------------------------------------
   INTEGER (KIND=LONG), PUBLIC, PARAMETER       :: nlvfepl=35
   REAL (kind=DP), PUBLIC, DIMENSION (nlvfepl)  :: gstfepl
   DATA gstfepl                  / 10.0_DP, 8.0_DP, 6.0_DP, 4.0_DP, 2.0_DP, & 
                                   10.0_DP, 8.0_DP, 6.0_DP, 4.0_DP, 8.0_DP, &
                                   6.0_DP, 4.0_DP, 2.0_DP, 6.0_DP, 4.0_DP,  & 
                                   2.0_DP, 6.0_DP, 4.0_DP, 2.0_DP, 14.0_DP, &
                                   12.0_DP, 10.0_DP, 8.0_DP, 10.0_DP, 8.0_DP, &
                                   6.0_DP, 4.0_DP, 12.0_DP, 10.0_DP, 8.0_DP, &
                                   6.0_DP, 8.0_DP, 6.0_DP, 4.0_DP, 2.0_DP /
   REAL (kind=DP), PUBLIC, DIMENSION (nlvfepl)  :: elkfepl
   DATA elkfepl                       / 0.0_DP, 553.95083_DP, 959.550155_DP, &
             1240.471089_DP, 1404.604668_DP, 2692.421982_DP, 3494.151388_DP, &
           4081.244575_DP, 4483.687486_DP, 11440.426112_DP, 12068.552694_DP, & 
          12483.62126_DP, 12723.508803_DP, 19378.809799_DP, 19664.465355_DP, &
          19997.46714_DP, 29957.534620_DP, 31370.030326_DP, 32228.575202_DP, & 
          30563.50_DP, 30820.75_DP, 31038.54_DP, 31224.77_DP, 32556.77_DP, &
          32804.55_DP, 32990.78_DP, 33123.35_DP, 36570.14_DP, 37113.05_DP, &
          37365.56_DP, 37471.30_DP, 45278.68_DP, 45141.37_DP, 45106.65_DP, &
          45112.96_DP /

   REAL (kind=DP), PUBLIC, DIMENSION (nlvfepl)  :: pop_fepl

   ! Observed transitions (microns).    
   !  1: 12-3 -> 1.248   12: 13-8 -> 1.664
   !  2: 10-1 -> 1.257   13: 11-7 -> 1.677
   !  3: 13-5 -> 1.271   14: 12-8 -> 1.711
   !  4: 12-4 -> 1.279   15: 13-9 -> 1.745
   !  5: 11-3 -> 1.295   16: 12-9 -> 1.798
   !  6: 12-5 -> 1.298   17: 11-8 -> 1.800
   !  7: 10-2 -> 1.321   18: 10-7 -> 1.810
   !  8: 11-4 -> 1.328   19: 7-6  -> 17.936
   !  9: 11-6 -> 1.534   20: 2-1  -> 25.988
   ! 10: 12-7 -> 1.600   21: 9-8  -> 35.777
   ! 11: 10-6 -> 1.644

   INTEGER (KIND=LONG), PUBLIC, PARAMETER      :: ntrfepl=256
   REAL (KIND=DP), PUBLIC, DIMENSION(ntrfepl)  :: wlkfepl
   REAL(KIND=DP),PUBLIC,DIMENSION(ntrfepl)     :: emifepl,emifepl_o,intfepl

   ! emissivity (erg/cm3/s), emissivity at last call to DRIVE, integrated intensity (erg/cm2/s/sr)

   ! ----------------------------------------------
   ! End of public database on levels and transitions
   ! ----------------------------------------------


   ! Coupling terms between fluids. Give global energy transfer by collisions
   !          cool_X < 0 if the fluid X looses energy
   REAL (kind=DP), PUBLIC  :: Cool_n, Cool_i, Cool_e  ! erg cm-3 s-1

   ! Minimum collision rate (used when nothing else is known)
   REAL (kind=DP), PUBLIC, PARAMETER :: colr_min = 1.0e-30_DP

   ! changes made in Nov.2002:

   ! * Collision matrix defined only inside of the subroutine that calculates
   ! the level populations
   ! => Avoids adding a new set of variables for each new atom/ion.
   ! => routine can easily be cut and pasted for each new atom
   ! * Teff's are declared here instead of inside each atomic routine
   ! * alpha's and Teff's will be calculated automatically (condition on charge_cool)
   !
   ! These are private variables, used in LINE_THERMAL_BALANCE and THERMAL_LOSS

   REAL (kind=DP) :: alp_H, alp_H2, alp_He, alp_Hp
   REAL (kind=DP) :: Teff_H, Teff_H2, Teff_He, Teff_Hp, Teff_e
   REAL (kind=DP) :: Dens_cool, mass_cool, charge_cool
   ! INTEGER (KIND=LONG), PRIVATE :: i


   !--------------------------------------------
   ! atomic cooling (of neutral fluid)
   !--------------------------------------------
   REAL(KIND=DP) :: total_n_cool   = 0.0_DP      ! total cooling term of neutrals  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Hat  = 0.0_DP      ! cooling rate of neutrals by H   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Cat  = 0.0_DP      ! cooling rate of neutrals by C   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Oat  = 0.0_DP      ! cooling rate of neutrals by O   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Nat  = 0.0_DP      ! cooling rate of neutrals by N   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Sat  = 0.0_DP      ! cooling rate of neutrals by S   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Siat = 0.0_DP      ! cooling rate of neutrals by Si  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Cp   = 0.0_DP      ! cooling rate of neutrals by C+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Op   = 0.0_DP      ! cooling rate of neutrals by O+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Np   = 0.0_DP      ! cooling rate of neutrals by N+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Sp   = 0.0_DP      ! cooling rate of neutrals by S+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Sip  = 0.0_DP      ! cooling rate of neutrals by Si+ (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_n_Fep  = 0.0_DP      ! cooling rate of neutrals by Fe+ (erg cm-3 s-1)

   !--------------------------------------------
   ! atomic cooling (of positively charged fluid)
   !--------------------------------------------
   REAL(KIND=DP) :: total_i_cool   = 0.0_DP      ! total cooling term of positive ions  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Hat  = 0.0_DP      ! cooling rate of positive ions by H   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Cat  = 0.0_DP      ! cooling rate of positive ions by C   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Oat  = 0.0_DP      ! cooling rate of positive ions by O   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Nat  = 0.0_DP      ! cooling rate of positive ions by N   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Sat  = 0.0_DP      ! cooling rate of positive ions by S   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Siat = 0.0_DP      ! cooling rate of positive ions by Si  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Cp   = 0.0_DP      ! cooling rate of positive ions by C+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Op   = 0.0_DP      ! cooling rate of positive ions by O+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Np   = 0.0_DP      ! cooling rate of positive ions by N+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Sp   = 0.0_DP      ! cooling rate of positive ions by S+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Sip  = 0.0_DP      ! cooling rate of positive ions by Si+ (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_i_Fep  = 0.0_DP      ! cooling rate of positive ions by Fe+ (erg cm-3 s-1)

   !--------------------------------------------
   ! atomic cooling (of negative fluid)
   !--------------------------------------------
   REAL(KIND=DP) :: total_e_cool   = 0.0_DP      ! total cooling term of negative ions  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Hat  = 0.0_DP      ! cooling rate of negative ions by H   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Cat  = 0.0_DP      ! cooling rate of negative ions by C   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Oat  = 0.0_DP      ! cooling rate of negative ions by O   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Nat  = 0.0_DP      ! cooling rate of negative ions by N   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Sat  = 0.0_DP      ! cooling rate of negative ions by S   (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Siat = 0.0_DP      ! cooling rate of negative ions by Si  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Cp   = 0.0_DP      ! cooling rate of negative ions by C+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Op   = 0.0_DP      ! cooling rate of negative ions by O+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Np   = 0.0_DP      ! cooling rate of negative ions by N+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Sp   = 0.0_DP      ! cooling rate of negative ions by S+  (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Sip  = 0.0_DP      ! cooling rate of negative ions by Si+ (erg cm-3 s-1)
   REAL(KIND=DP) :: cooling_e_Fep  = 0.0_DP      ! cooling rate of negative ions by Fe+ (erg cm-3 s-1)



CONTAINS



  SUBROUTINE BUILD_TR_NAME
     IMPLICIT none
     CHARACTER(LEN=1) :: dum1
     CHARACTER(LEN=1) :: dum2
     INTEGER          :: i
     DO i = 1, ntrhat 
        WRITE(dum1,'(I1)') iuphat(i)
        WRITE(dum2,'(I1)') jlohat(i)
        namtrhat(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrhat(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrcat 
        WRITE(dum1,'(I1)') iupcat(i)
        WRITE(dum2,'(I1)') jlocat(i)
        namtrcat(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrcat(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrnat 
        WRITE(dum1,'(I1)') iupnat(i)
        WRITE(dum2,'(I1)') jlonat(i)
        namtrnat(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrnat(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntroat 
        WRITE(dum1,'(I1)') iupoat(i)
        WRITE(dum2,'(I1)') jlooat(i)
        namtroat(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtroat(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrsat 
        WRITE(dum1,'(I1)') iupsat(i)
        WRITE(dum2,'(I1)') jlosat(i)
        namtrsat(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrsat(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrsiat
        WRITE(dum1,'(I1)') iupsiat(i)
        WRITE(dum2,'(I1)') jlosiat(i)
        namtrsiat(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrsiat(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrcpl 
        WRITE(dum1,'(I1)') iupcpl(i)
        WRITE(dum2,'(I1)') jlocpl(i)
        namtrcpl(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrcpl(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrnpl 
        WRITE(dum1,'(I1)') iupnpl(i)
        WRITE(dum2,'(I1)') jlonpl(i)
        namtrnpl(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrnpl(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntropl 
        WRITE(dum1,'(I1)') iupopl(i)
        WRITE(dum2,'(I1)') jloopl(i)
        namtropl(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtropl(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrspl 
        WRITE(dum1,'(I1)') iupspl(i)
        WRITE(dum2,'(I1)') jlospl(i)
        namtrspl(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrspl(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
     DO i = 1, ntrsipl
        WRITE(dum1,'(I1)') iupsipl(i)
        WRITE(dum2,'(I1)') jlosipl(i)
        namtrsipl(i) = TRIM(dum1)//' -> '//TRIM(dum2)
        idtrsipl(i)  = TRIM(dum1)//'_'//TRIM(dum2)
     ENDDO
  END SUBROUTINE BUILD_TR_NAME



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLHAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  use tex
  implicit none

  REAL (kind=DP) :: eff_col_str                           ! effective collision strength
  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

! ec matrices will be multiplied by collisioner density in THERMAL_LOSS
! units: cm3 s-1

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

!--- effective temperatures with electrons : Te

! collisions with electrons
  eff_col_str = 9.479e-4_DP * Te ** 0.5582_DP
  ec_e(1,2) = 8.629e-6_DP * eff_col_str / gsthat(2) / sqrt(Te)
  eff_col_str = 8.827e-2_DP * Te ** 0.1293_DP
  ec_e(1,3) = 8.629e-6_DP * eff_col_str / gsthat(3) / sqrt(Te)
  eff_col_str = 1.931e-3_DP * Te ** 0.5565_DP
  ec_e(1,4) = 8.629e-6_DP * eff_col_str / gsthat(4) / sqrt(Te)
  eff_col_str = 5.424e-4_DP * Te ** 0.4647_DP

  ec_e(1,5) = 8.629e-6_DP * eff_col_str / gsthat(5) / sqrt(Te)
  eff_col_str = 3.458e-2_DP * Te ** 0.3459_DP
  ec_e(2,5) = 8.629e-6_DP * eff_col_str / gsthat(5) / sqrt(Te)
  eff_col_str = 2.029e-3_DP * Te ** 0.6646_DP
  ec_e(3,5) = 8.629e-6_DP * eff_col_str / gsthat(5) / sqrt(Te)
  eff_col_str = 6.915e-3_DP * Te ** 0.3459_DP
  ec_e(4,5) = 8.629e-6_DP * eff_col_str / gsthat(5) / sqrt(Te)

  eff_col_str = 2.499e-2_DP * Te ** 0.1094_DP
  ec_e(1,6) = 8.629e-6_DP * eff_col_str / gsthat(6) / sqrt(Te)
  eff_col_str = 2.085e-1_DP * Te ** 0.1353_DP
  ec_e(2,6) = 8.629e-6_DP * eff_col_str / gsthat(6) / sqrt(Te)
  eff_col_str = 1.689e-2_DP * Te ** 0.4757_DP
  ec_e(3,6) = 8.629e-6_DP * eff_col_str / gsthat(6) / sqrt(Te)
  eff_col_str = 4.112e-1_DP * Te ** 0.1368_DP
  ec_e(4,6) = 8.629e-6_DP * eff_col_str / gsthat(6) / sqrt(Te)

  eff_col_str = 3.447e-3_DP * Te ** 0.2171_DP
  ec_e(1,7) = 8.629e-6_DP * eff_col_str / gsthat(7) / sqrt(Te)
  eff_col_str = 3.276e-3_DP * Te ** 0.7035_DP
  ec_e(2,7) = 8.629e-6_DP * eff_col_str / gsthat(7) / sqrt(Te)
  eff_col_str = 3.961e-3_DP * Te ** 0.6131_DP
  ec_e(3,7) = 8.629e-6_DP * eff_col_str / gsthat(7) / sqrt(Te)
  eff_col_str = 6.663e-3_DP * Te ** 0.7019_DP
  ec_e(4,7) = 8.629e-6_DP * eff_col_str / gsthat(7) / sqrt(Te)

  eff_col_str = 1.085e-3_DP * Te ** 0.4647_DP
  ec_e(1,8) = 8.629e-6_DP * eff_col_str / gsthat(8) / sqrt(Te)
  eff_col_str = 6.915e-2_DP * Te ** 0.3459_DP
  ec_e(2,8) = 8.629e-6_DP * eff_col_str / gsthat(8) / sqrt(Te)
  eff_col_str = 4.009e-3_DP * Te ** 0.6654_DP
  ec_e(3,8) = 8.629e-6_DP * eff_col_str / gsthat(8) / sqrt(Te)
  eff_col_str = 1.387e-1_DP * Te ** 0.3455_DP
  ec_e(4,8) = 8.629e-6_DP * eff_col_str / gsthat(8) / sqrt(Te)

  eff_col_str = 5.198e-2_DP * Te ** 0.2167_DP
  ec_e(1,9) = 8.629e-6_DP * eff_col_str / gsthat(9) / sqrt(Te)
  eff_col_str = 5.023e-3_DP * Te ** 0.7013_DP
  ec_e(2,9) = 8.629e-6_DP * eff_col_str / gsthat(9) / sqrt(Te)
  eff_col_str = 5.976e-3_DP * Te ** 0.6128_DP
  ec_e(3,9) = 8.629e-6_DP * eff_col_str / gsthat(9) / sqrt(Te)
  eff_col_str = 9.939e-3_DP * Te ** 0.7025_DP
  ec_e(4,9) = 8.629e-6_DP * eff_col_str / gsthat(9) / sqrt(Te)


  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  end subroutine COLHAT


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLCAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  use tex
  implicit none

  REAL (kind=DP) :: auxt4                           ! auxiliary factor
  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

! ec matrices will be multiplied by collisioner density in THERMAL_LOSS
! units: cm3 s-1

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do


! Build collision matrix for atomic carbon
! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

! collisions with H
! Fine structure: Ref unknown

  ec_H(1,2) = 1.01e-10_DP * Teff_H**0.117_DP
  ec_H(1,3) = 4.49e-11_DP * Teff_H**0.194_DP
  ec_H(2,3) = 1.06e-10_DP * Teff_H**0.234_DP

! 1D - 3P
! DRF : Do not use O I + H rates

! collisions with H2
! fit from Schroder, Staemmler, Smith, Flower & Jaquet (1991)

  ec_paraH2(1,2)  = 0.80e-10_DP
  ec_orthoH2(1,2) = 0.75e-10_DP
  ec_paraH2(1,3)  = 0.90e-10_DP
  ec_orthoH2(1,3) = 3.54e-11_DP * Teff_H2**0.167_DP
  ec_paraH2(2,3)  = 2.00e-10_DP
  ec_orthoH2(2,3) = 5.25e-11_DP * Teff_H2**0.244_DP

! collisions with He
! Fine structure: Staemmler & Flower (1991)

  ec_He(1,2) = 8.38e-12_DP * Teff_He**0.159_DP
  ec_He(1,3) = 5.98e-11_DP * Teff_He**0.078_DP
  ec_He(2,3) = 3.68e-11_DP * Teff_He**0.041_DP

! collisions avec H+
! from AA 236, 515 (1990)  Roueff & Le Bourlot (fit DRF)
  ec_Hp(1,2) = tpow(-10.359_DP + 0.7959_DP*log10(Teff_Hp) - 0.08748_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(1,3) = tpow(-13.232_DP + 2.4171_DP*log10(Teff_Hp) - 0.29151_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(2,3) = tpow(-11.290_DP + 1.7915_DP*log10(Teff_Hp) - 0.23010_DP*(log10(Teff_Hp))**2._DP)

!--- effective temperatures with electrons : Te

! collisions with electrons
! Should we use Pequignot 1990...
! From Pequignot & Aldrovandi, A&A 50, 141, 1976 as compiled by Mendoza (1983)

  auxt4 = 1.83e-4_DP * Te ** 0.444_DP
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstcat(4) * gstcat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstcat(4) * gstcat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstcat(4) * gstcat(3) / 9.0_DP
  auxt4 = 9.86e-5_DP * Te ** 0.343_DP
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstcat(5) * gstcat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstcat(5) * gstcat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstcat(5) * gstcat(3) / 9.0_DP
  auxt4 = 2.77e-3_DP
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstcat(5)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  end subroutine COLCAT


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLOAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2
  REAL (kind=DP)                      :: toto


  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic oxygen
! Energy levels.
!        1 - 3P J = 2
!        2 - 3P J = 1
!        3 - 3P J = 0
!        4 - 1D J = 2
!        5 - 1S J = 0


! All collision rates should be checked against Pequignot 1990

! collisions with H

! Fine structure

!!$  ec_H(1,2) = 4.37e-12_dp * Teff_H**0.66_dp
!!$  ec_H(1,3) = 1.06e-12_dp * Teff_H**0.80_dp
!!$  ec_H(2,3) = 1.35e-11_dp * Teff_H**0.45_dp

!  fit from Abrahamsson et al  2007, ApJ654, 1171 (done by Evelyne - XII 08)
  ec_H(1,2) = 5.9601e-11_dp * Teff_H**0.39561_dp
  ec_H(1,3) = 6.1719e-11_dp * Teff_H**0.36291_dp
  IF (Teff_H <= 1000.0_dp) THEN
     ec_H(2,3) = (2.5544e-11_dp - 1.232e-14_dp * Teff_H) * Teff_H**0.62996_dp
  ELSE
     ec_H(2,3) =  1.063e-9_dp
  ENDIF


! 1D - 3P
! fit by JLB, from datas of Federman & Shipsey, Ap J, 1983, 269, 791 (using potential V2)

  ec_H(1,4) = 9.85e-14_DP * Teff_H**0.245_DP
  ec_H(2,4) = 9.85e-14_DP * Teff_H**0.245_DP
  ec_H(3,4) = 9.85e-14_DP * Teff_H**0.245_DP

! collisions with H2
! fit from Jaquet, Staemmler, Smith & Flower, J.Phys.B (1991)

  ec_paraH2(1,2)  = 3.46e-11_DP * Teff_H2**0.316_DP
  ec_orthoH2(1,2) = 2.70e-11_DP * Teff_H2**0.362_DP
  ec_paraH2(1,3)  = 7.07e-11_DP * Teff_H2**0.268_DP
  ec_orthoH2(1,3) = 5.49e-11_DP * Teff_H2**0.317_DP
  ec_paraH2(2,3)  = 1.44e-14_DP * Teff_H2**1.109_DP
  ec_orthoH2(2,3) = 4.64e-14_DP * Teff_H2**0.976_DP


! collisions with He
! fit from Monteiro & Flower, MNRAS 228, 101 (1987)

  ec_He(1,2) = 1.55e-12_DP * Teff_He**0.694_DP
  ec_He(1,3) = 2.52e-12_DP * Teff_He**0.681_DP
  ec_He(2,3) = 2.00e-15_DP * Teff_He**1.428_DP

! collisions avec H+

! guess (Evelyne)
! ec_Hp(1,2) = 1.0e-8_DP
! ec_Hp(1,3) = 1.0e-8_DP
! ec_Hp(2,3) = 1.0e-8_DP

! Pequignot A&A 231, 499 (1990) + erratum A&A 313, 1026 (1996)

  if (Teff_Hp < 194.664_dp) then
    toto = 1.40e-3_dp * (Teff_Hp*1.0e-2_dp)**0.90_dp
  else if (Teff_Hp < 3686.2414_dp) then
    toto = 2.14e-2_dp * (Teff_Hp*1.0e-3_dp)**1.30_dp
  else
    toto = 2.78e-1_dp * (Teff_Hp*1.0e-4_dp)**0.87_dp
  endif
  ec_Hp(1,2) = toto * 8.629e-6_dp / (sqrt(Teff_Hp) * gstoat(2))

  if (Teff_Hp < 511.9563_dp) then
    toto = 1.12e-4_dp * (Teff_Hp*1.0e-2_dp)**1.60_dp
  else if (Teff_Hp < 7510.155_dp) then
    toto = 3.90e-3_dp * (Teff_Hp*1.0e-3_dp)**1.40_dp
  else
    toto = 8.25e-2_dp * (Teff_Hp*1.0e-4_dp)**0.80_dp
  endif
  ec_Hp(1,3) = toto * 8.629e-6_dp / sqrt(Teff_Hp)

  if (Teff_Hp < 2090.2558_dp) then
    toto = 3.10e-4_dp * (Teff_Hp*1.0e-2_dp)**1.06_dp
  else
    toto = 2.29e-2_dp * (Teff_Hp*1.0e-4_dp)**0.69_dp
  endif
  ec_Hp(2,3) = toto * 8.629e-6_dp / sqrt(Teff_Hp)

!--- effective temperatures with electrons : Te

! collisions with electrons
! From Pequignot, 1990, A&A, 231, 499 and Berrington, 1998, MNRAS, 293, L83

  auxt4 = 4.28e-1_DP * (Te/1.e4_DP) ** 1.43_DP / (1.00_DP + 6.05e-1_DP * (Te/1.e4_DP) ** 1.105_DP)
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstoat(4) / sqrt(Te) * gstoat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstoat(4) / sqrt(Te) * gstoat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstoat(4) / sqrt(Te) * gstoat(3) / 9.0_DP
  auxt4 = 5.88e-2_DP * (Te/1.e4_DP) ** 1.50_DP / (1.00_DP + 8.00e-1_DP * (Te/1.e4_DP) ** 1.125_DP)
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstoat(5) / sqrt(Te) * gstoat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstoat(5) / sqrt(Te) * gstoat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstoat(5) / sqrt(Te) * gstoat(3) / 9.0_DP
  auxt4 = 1.16e-1_DP * (Te/1.e4_DP) ** 0.53_DP / (1.00_DP + 1.11e-1_DP * (Te/1.e4_DP) ** 0.160_DP)
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstoat(5) / sqrt(Te)

  auxt4 = 1.49e-4_DP * Te ** 0.4565_DP
  ec_e(1,2) = 8.629e-6_DP * auxt4 / gstoat(2) / sqrt(Te)
  auxt4 = 4.98e-5_DP * Te ** 0.4955_DP
  ec_e(1,3) = 8.629e-6_DP * auxt4 / gstoat(3) / sqrt(Te)
  auxt4 = 1.83e-9_DP * Te ** 1.347_DP
  ec_e(2,3) = 8.629e-6_DP * auxt4 / gstoat(3) / sqrt(Te)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLOAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLNAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic nitrogen
! Energy levels.
!        1 - 4So J = 3/2
!        2 - 2Do J = 5/2
!        3 - 2Do J = 3/2
!        4 - 2Po J = 1/2
!        5 - 2Po J = 3/2

! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_orthoH2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions avec H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperatures with electrons : Te

! collisions with electrons
! From compilation of Mendoza (1983)

  auxt4 = 7.92e-6_DP * Te**0.589_DP
  ec_e(1,3) = 8.629e-6_DP * auxt4 / gstnat(3)
  auxt4 = 1.5_DP * auxt4
  ec_e(1,2) = 8.629e-6_DP * auxt4 / gstnat(2)
  auxt4 = 1.74e-6_DP * Te**0.621_DP
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstnat(4)
  auxt4 = 2.0_DP * auxt4
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 4.78e-5_DP * Te**0.431_DP
  ec_e(2,3) = 8.629e-6_DP * auxt4 / gstnat(3)
  auxt4 = 2.61e-6_DP * Te**0.609_DP
  ec_e(4,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 3.59e-4_DP * Te**0.217_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 1.13e-4_DP * Te**0.279_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstnat(5)
  auxt4 = 6.82e-5_DP * Te**0.301_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstnat(4)
  auxt4 = 1.65e-4_DP * Te**0.193_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstnat(4)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLNAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLSAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic sulfur

! Energy levels.
!        1 - 3P J = 3
!        2 - 3P J = 1
!        3 - 3P J = 0
!        4 - 1D J = 2
!        5 - 1S J = 0

! All collision rates as for O, but CT with H is excluded, as it is non-resonant

! collisions with H

! Fine structure

  ec_H(1,2) = 4.37e-12_DP * Teff_H**0.66_DP
  ec_H(1,3) = 1.06e-12_DP * Teff_H**0.80_DP
  ec_H(2,3) = 1.35e-11_DP * Teff_H**0.45_DP

! collisions with H2
! fit from Jaquet, Staemmler, Smith & Flower, J.Phys.B (1991)

  ec_paraH2(1,2) = 3.46e-11_DP * Teff_H2**0.316_DP
  ec_orthoH2(1,2) = 2.70e-11_DP * Teff_H2**0.362_DP
  ec_paraH2(1,3) = 7.07e-11_DP * Teff_H2**0.268_DP
  ec_orthoH2(1,3) = 5.49e-11_DP * Teff_H2**0.317_DP
  ec_paraH2(2,3) = 1.44e-14_DP * Teff_H2**1.109_DP
  ec_orthoH2(2,3) = 4.64e-14_DP * Teff_H2**0.976_DP

! collisions with He
! fit from Monteiro & Flower, MNRAS 228, 101 (1987)

  ec_He(1,2) = 1.55e-12_DP * Teff_He**0.694_DP
  ec_He(1,3) = 2.52e-12_DP * Teff_He**0.681_DP
  ec_He(2,3) = 2.00e-15_DP * Teff_He**1.428_DP

! collisions avec H+
! guess (Evelyne)

  ec_Hp(1,2) = 1.0e-8_DP
  ec_Hp(1,3) = 1.0e-8_DP
  ec_Hp(2,3) = 1.0e-8_DP

!--- effective temperatures with electrons : Te

! collisions with electrons
! From Pequignot, 1990, A&A, 231, 499 and Berrington, 1998, MNRAS, 293, L83

  auxt4 = 4.28e-1_DP * (Te/1.e4_DP) ** 1.43_DP / (1.00_DP + 6.05e-1_DP * (Te/1.e4_DP) ** 1.105_DP)
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstsat(4) / sqrt(Te) * gstsat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstsat(4) / sqrt(Te) * gstsat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstsat(4) / sqrt(Te) * gstsat(3) / 9.0_DP
  auxt4 = 5.88e-2_DP * (Te/1.e4_DP) ** 1.50_DP / (1.00_DP + 8.00e-1_DP * (Te/1.e4_DP) ** 1.125_DP)
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstsat(5) / sqrt(Te) * gstsat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstsat(5) / sqrt(Te) * gstsat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstsat(5) / sqrt(Te) * gstsat(3) / 9.0_DP
  auxt4 = 1.16e-1_DP * (Te/1.e4_DP) ** 0.53_DP / (1.00_DP + 1.11e-1_DP * (Te/1.e4_DP) ** 0.160_DP)
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstsat(5) / sqrt(Te)

  auxt4 = 1.49e-4_DP * Te ** 0.4565_DP
  ec_e(1,2) = 8.629e-6_DP * auxt4 / gstsat(2) / sqrt(Te)
  auxt4 = 4.98e-5_DP * Te ** 0.4955_DP
  ec_e(1,3) = 8.629e-6_DP * auxt4 / gstsat(3) / sqrt(Te)
  auxt4 = 1.83e-9_DP * Te ** 1.347_DP
  ec_e(2,3) = 8.629e-6_DP * auxt4 / gstsat(3) / sqrt(Te)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLSAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLSiAT(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  use tex
  implicit none

  REAL (kind=DP) :: auxt4
  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv)= colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for atomic silicium
! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

! All collision rates as for C

! collisions with H
! Fine structure: Launay & Roueff, 1977, A&A, 56, 289

  ec_H(1,2) = 1.01e-10_DP * Tn**0.117_DP
  ec_H(1,3) = 4.49e-11_DP * Tn**0.194_DP
  ec_H(2,3) = 1.06e-10_DP * Tn**0.234_DP

! collisions with H2
! fit from Schroder, Staemmler, Smith, Flower & Jaquet (1991)

  ec_paraH2(1,2) = 0.80e-10_DP
  ec_orthoH2(1,2)= 0.75e-10_DP
  ec_paraH2(1,3) = 0.90e-10_DP
  ec_orthoH2(1,3)= 3.54e-11_DP * Teff_H2**0.167_DP
  ec_paraH2(2,3) = 2.00e-10_DP
  ec_orthoH2(2,3)= 5.25e-11_DP * Teff_H2**0.244_DP

! collisions with He
! Fine structure: Staemmler & Flower (1991)

  ec_He(1,2) = 8.38e-12_DP * Teff_He**0.159_DP
  ec_He(1,3) = 5.98e-11_DP * Teff_He**0.078_DP
  ec_He(2,3) = 3.68e-11_DP * Teff_He**0.041_DP

! collisions avec H+

  ec_Hp(1,2) = tpow(-10.359_DP + 0.7959_DP*log10(Teff_Hp) - 0.08748_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(1,3) = tpow(-13.232_DP + 2.4171_DP*log10(Teff_Hp) - 0.29151_DP*(log10(Teff_Hp))**2._DP)
  ec_Hp(2,3) = tpow(-11.290_DP + 1.7915_DP*log10(Teff_Hp) - 0.23010_DP*(log10(Teff_Hp))**2._DP)

!--- effective temperatures with electrons : Te

! collisions with electrons
! Should we use Pequignot 1990...
! From Pequignot & Aldrovandi, A&A 50, 141, 1976 as compiled by Mendoza (1983)

  auxt4 = 1.83e-4_DP * Te ** 0.444_DP
  ec_e(1,4) = 8.629e-6_DP * auxt4 / gstsiat(4) * gstsiat(1) / 9.0_DP
  ec_e(2,4) = 8.629e-6_DP * auxt4 / gstsiat(4) * gstsiat(2) / 9.0_DP
  ec_e(3,4) = 8.629e-6_DP * auxt4 / gstsiat(4) * gstsiat(3) / 9.0_DP
  auxt4 = 9.86e-5_DP * Te ** 0.343_DP
  ec_e(1,5) = 8.629e-6_DP * auxt4 / gstsiat(5) * gstsiat(1) / 9.0_DP
  ec_e(2,5) = 8.629e-6_DP * auxt4 / gstsiat(5) * gstsiat(2) / 9.0_DP
  ec_e(3,5) = 8.629e-6_DP * auxt4 / gstsiat(5) * gstsiat(3) / 9.0_DP
  auxt4 = 2.77e-3_DP
  ec_e(4,5) = 8.629e-6_DP *auxt4 / gstsiat(5)

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  end subroutine COLSiAT

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLCPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2
  REAL (kind=DP)                      :: toto

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv)  = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionised carbon
! Energy levels.
!        1 - 2P J = 1/2
!        2 - 2P J = 3/2
!        3 - 4P J = 1/2
!        4 - 4P J = 3/2
!        5 - 4P J = 5/2


! collisions with H

  ec_H(1,2) = 8.86e-10_dp

! collisions with H2

  ec_paraH2(1,2) = 3.94e-10_dp
  ec_orthoH2(1,2)= 4.92e-10_dp

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions avec H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperatures with electrons : Te

! collisions with electrons
! k21 : Hayes & Nussbaumer, 1984, A&A 134, 193
! Dec 02, CM & DRF
! k21 changed to power law fit to Wilson & Bell, 2001, MNRAS, 37, 1027
! other : Lennon et al., 1985, AJ 294, 200

  toto = 0.65582_dp * Te**(0.13244_DP)
  ec_e(1,2) = 8.629e-6_DP * min(toto,2.5_dp) / (4.0_DP*sqrt(Te))
  ec_e(1,3) = 8.629e-6_DP * (0.288_DP - 5.89e-7_DP*Te) / (2.0_DP * sqrt(Te))
  ec_e(1,4) = 8.629e-6_DP * (0.424_DP - 8.35e-7_DP*Te) / (4.0_DP * sqrt(Te))
  ec_e(1,5) = 8.629e-6_DP * (0.257_DP - 3.95e-7_DP*Te) / (6.0_DP * sqrt(Te))
  ec_e(2,3) = 8.629e-6_DP * (0.197_DP - 3.15e-7_DP*Te) / (2.0_DP * sqrt(Te))
  ec_e(2,4) = 8.629e-6_DP * (0.544_DP - 9.84e-7_DP*Te) / (4.0_DP * sqrt(Te))
  ec_e(2,5) = 8.629e-6_DP * (1.200_DP - 2.40e-6_DP*Te) / (6.0_DP * sqrt(Te))

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLCPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLNPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2  = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv)  = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionised nitrogen
! Energy levels.
!        1 - 3P J = 0
!        2 - 3P J = 1
!        3 - 3P J = 2
!        4 - 1D J = 2
!        5 - 1S J = 0

! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_H2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions avec H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperatures with electrons : Te

! collisions with electrons
! Mendoza, 1983, International astrophysical union, symposium 103

  ec_e(1,2) =  8.629e-6_DP * (0.138_DP * Te**0.119_DP)   / (3.0_DP * sqrt(Te))
  ec_e(1,3) =  8.629e-6_DP * (0.0834_DP * Te**0.130_DP)  / (5.0_DP * sqrt(Te))
  ec_e(2,3) =  8.629e-6_DP * (0.359_DP*Te**0.125_DP)     / (5.0_DP * sqrt(Te))
  ec_e(1,4) =  8.629e-6_DP * (0.203_DP * Te**0.0414_DP)  / (5.0_DP * sqrt(Te))
  ec_e(2,4) =  8.629e-6_DP * (0.609_DP * Te**0.0414_DP)  / (5.0_DP * sqrt(Te))
  ec_e(3,4) =  8.629e-6_DP * (1.02_DP * Te**0.0414_DP)   / (5.0_DP * sqrt(Te))
  ec_e(4,5) =  8.629e-6_DP * (3.14_DP*Te**(-0.142_DP))     / (1.0_DP * sqrt(Te))
  ec_e(3,5) =  8.629e-6_DP * (0.158_DP + Te*4.26e-7_DP)  / (1.0_DP * sqrt(Te))
  ec_e(2,5) =  8.629e-6_DP * (0.0949_DP + Te*2.56e-7_DP) / (1.0_DP * sqrt(Te))
  ec_e(1,5) =  8.629e-6_DP * (0.0316_DP + Te*8.52e-8_DP) / (1.0_DP * sqrt(Te))

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLNPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLOPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2 = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv)  = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionised oxygen
! Energy levels.
!        1 - 4So J = 3/2
!        2 - 2Do J = 5/2
!        3 - 2Do J = 3/2
!        4 - 2Po J = 3/2
!        5 - 2Po J = 1/2


! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_H2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions avec H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperatures with electrons : Te

! collisions with electrons

! ec_e(1,2) = 0.0_DP

  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLOPL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine COLSPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MODULE_PHYS_VAR
  implicit none

  INTEGER        :: i

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
  REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

  call ALPHATEFF

  ec_H  = 0.0_DP
  ec_orthoH2 = 0.0_DP
  ec_paraH2 = 0.0_DP
  ec_He = 0.0_DP
  ec_Hp = 0.0_DP
  ec_e  = 0.0_DP
  do i=1,nlv-1
    ec_H(i,i+1:nlv)  = colr_min
    ec_orthoH2(i,i+1:nlv) = colr_min
    ec_paraH2(i,i+1:nlv) = colr_min
    ec_He(i,i+1:nlv) = colr_min
    ec_Hp(i,i+1:nlv) = colr_min
    ec_e (i,i+1:nlv) = colr_min
  end do

! Build collision matrix for ionised sulfur
! Energy levels.
!        1 - 4So J = 3/2
!        2 - 2Do J = 3/2
!        3 - 2Do J = 5/2
!        4 - 2Po J = 1/2
!        5 - 2Po J = 3/2

! collisions with H

! ec_H(1,2) = 0.0_DP

! collisions with H2

! ec_H2(1,2) = 0.0_DP

! collisions with He

! ec_He(1,2) = 0.0_DP

! collisions avec H+

! ec_Hp(1,2) = 0.0_DP

!--- effective temperatures with electrons : Te

! collisions with electrons
! Keenan et al. 1996, MNRAS, 281, 1073

  ec_e(1,2) =  8.629e-6_DP * (12.195_DP * Te**(-0.16174_DP)) /(4.0_DP*sqrt(Te))
  ec_e(1,3) =  8.629e-6_DP * (18.311_DP * Te**(-0.16186_DP)) /(6.0_DP*sqrt(Te))
  ec_e(1,4) =  8.629e-6_DP * (1.2736_DP - (1.0000e-5_DP*Te)) /(2.0_DP*sqrt(Te))
  ec_e(1,5) =  8.629e-6_DP * (2.5629_DP - (2.1143e-5_DP*Te)) /(4.0_DP*sqrt(Te))
  ec_e(2,3) =  8.629e-6_DP * (8.3293_DP - (8.8206e-5_DP*Te)) /(6.0_DP*sqrt(Te))
  ec_e(2,4) =  8.629e-6_DP * (8.0121_DP * Te**(-0.16233_DP)) /(2.0_DP*sqrt(Te))
  ec_e(2,5) =  8.629e-6_DP * (9.3033_DP * Te**(-0.12316_DP)) /(4.0_DP*sqrt(Te))
  ec_e(3,4) =  8.629e-6_DP * (6.4954_DP * Te**(-0.11783_DP)) /(2.0_DP*sqrt(Te))
  ec_e(3,5) =  8.629e-6_DP * (19.378_DP * Te**(-0.14738_DP)) /(4.0_DP*sqrt(Te))
  ec_e(4,5) =  8.629e-6_DP * (11.584_DP * Te**(-0.15835_DP)) /(4.0_DP*sqrt(Te))


  call THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  return
  end subroutine COLSPL

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COLSiPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     USE MODULE_PHYS_VAR
     IMPLICIT none

     INTEGER        :: i

     INTEGER       , intent(in)                     :: ntr, nlv
     INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
     REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
     REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
     REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
     REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
     REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
     REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2

     CALL ALPHATEFF

     ec_H  = 0.0_DP
     ec_orthoH2 = 0.0_DP
     ec_paraH2 = 0.0_DP
     ec_He = 0.0_DP
     ec_Hp = 0.0_DP
     ec_e  = 0.0_DP
     DO i=1,nlv-1
       ec_H(i,i+1:nlv)  = colr_min
       ec_orthoH2(i,i+1:nlv) = colr_min
       ec_paraH2(i,i+1:nlv) = colr_min
       ec_He(i,i+1:nlv) = colr_min
       ec_Hp(i,i+1:nlv) = colr_min
       ec_e (i,i+1:nlv) = colr_min
     ENDDO

     ! Build collision matrix for ionised silicium
     ! Energy levels.
     !        1 - 2P J = 1/2
     !        2 - 2P J = 3/2
     !        3 - 4P J = 1/2
     !        4 - 4P J = 3/2
     !        5 - 4P J = 5/2

     ! collisions with H
      ec_H(1,2) = 6.00e-10_DP

     ! collisions with H2
     ! Same as C+
     ec_paraH2(1,2) = 3.94e-10_DP
     ec_orthoH2(1,2)= 4.92e-10_DP

     ! collisions with He
     ! ec_He(1,2) = 0.0_DP

     ! collisions avec H+
     ! ec_Hp(1,2) = 0.0_DP
     
     !--- effective temperatures with electrons : Te

     ! collisions with electrons
     ! Dufton & Kingston, 1991, MNRAS 248, 827
     ec_e(1,2) =  8.629e-6_DP * 5.6_DP / (4.0_DP*SQRT(Te))

     CALL THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                          ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  END SUBROUTINE COLSiPL

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  SUBROUTINE COLFEPL(ntr,nlv,gst,elk,wlk,aij,iup,jlo,emi,pop)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! The collision strengths of Zhang and Pradhan (95) were used to calculate
  ! rates due to electrons.  Radiative decay rates used are those of
  ! Nussbaumer and Storey (88) for the first four terms and Quinet, Le
  ! Dourneuf and Zeippen (96) for the rest.
  !
  ! Deexcitation rates for collisions with H, H2, He, H+ calculated using the
  ! orbiting approximation and assuming a probability of 1/2 for each transition.

     USE MODULE_PHYS_VAR
     USE MODULE_READ_FE_DATA, ONLY : gamfepl
     USE MODULE_CONSTANTS
     USE MODULE_CHEMICAL_SPECIES
     IMPLICIT none

     INTEGER        :: i

     INTEGER       , intent(in)                     :: ntr, nlv
     INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
     REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: aij
     REAL (kind=DP),              DIMENSION (ntr)   :: wlk
     REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
     REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi
     REAL (kind=DP), intent(out), DIMENSION (nlv)   :: pop
     REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He, ec_Hp, ec_e
     REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2,ec_paraH2
     REAL (kind=DP) :: H2_red, H_red, He_red  ! reduced masses (AU)


     DO i=1, ntr 
        wlk(i) = elk(iup(i)) - elk(jlo(i))
     ENDDO

     ! ec matrices will be multiplied by collisioner density in THERMAL_LOSS
     ! units: cm3 s-1

     ec_H  = 0.0_DP
     ec_orthoH2 = 0.0_DP
     ec_paraH2  = 0.0_DP
     ec_He = 0.0_DP
     ec_Hp = 0.0_DP
     ec_e  = 0.0_DP
     DO i=1,nlv-1
        ec_H(i,i+1:nlv)  = colr_min
        ec_orthoH2(i,i+1:nlv)= colr_min
        ec_paraH2(i,i+1:nlv) = colr_min
        ec_He(i,i+1:nlv) = colr_min
        ec_Hp(i,i+1:nlv) = colr_min
        ec_e (i,i+1:nlv) = colr_min
     ENDDO

     ! Reduced masses
     H2_red =(mP/me) * (mass_Feplus/amu*mass_H2/amu) / (mass_Feplus/amu+mass_H2/amu)  
     H_red  =(mP/me) * (mass_Feplus/amu*mass_H/amu) / (mass_Feplus/amu+mass_H/amu)
     He_red =(mP/me) * (mass_Feplus/amu*mass_He/amu) / (mass_Feplus/amu+mass_He/amu)

    ! H2, He
     DO i=1,4
        ec_orthoH2(i,i+1:5) = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red)
        ec_paraH2(i,i+1:5)  = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red) 
        ec_He(i,i+1:5)      = AUcgs*pi* sqrt((alpha_He/bohr**3)/He_red)
     ENDDO

     DO i=6,nlv-1
         ec_orthoH2(i,i+1:nlv) = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red) 
         ec_paraH2(i,i+1:nlv)  = AUcgs*pi* sqrt((alpha_H2/bohr**3)/H2_red) 
           ec_He(i,i+1:nlv)    = AUcgs*pi* sqrt((alpha_He/bohr**3)/He_red) 
     ENDDO


     DO i=1,nlv-1
        !  H                                                                      
        ec_H(i,i+1:nlv) = AUcgs*pi* sqrt((alpha_H/bohr**3)/H_red)   

        ! Electron collisional deexcitation rates (cm3 s-1)
        ec_e(i,i+1:nlv) = 8.629D-6 *  gamfepl(i,i+1:nlv) / (gst(i)*sqrt(Te)) 
     ENDDO

     CALL THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                     ec_paraH2,ec_He,ec_Hp, ec_e,emi,pop,ntr,nlv)

  END SUBROUTINE COLFEPL

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE LINE_THERMAL_BALANCE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Nov. 2002: calls routines that compute the atomic/ionic
  ! line emissivities and the corresponding thermal loss terms
  !
  ! Principle: - one routine COLxxx for each atom, where specific
  !  collision deexcitation rates are calculated; this routine
  !  calls:  ALPHATEFF (calculates alph's and Teff's)
  !          THERMAL_LOSS (matrix inversion and cooling terms)
  ! To add a new coolant xxx:
  !      - set Dens_cool, mass_cool, and charge_cool for specy xxx
  !      - copy routine COLCAT into a new routine COLxxx. only
  !  numerical expressions for the rates need to be changed.
  !      - replace 'cat' with 'xxx' in arguments to COLxxx.
  !      - give atomic parameters at beginning of this module
  !      - that's it !

     USE module_phys_var,only:cool_KN
     USE MODULE_CHEMICAL_SPECIES
     USE MODULE_READ_FE_DATA
     USE MODULE_MOLECULAR_COOLING
     USE MODULE_VAR_TH_BALANCE
     IMPLICIT none

     REAL(KIND=DP) :: prev_cool_n
     REAL(KIND=DP) :: prev_cool_i
     REAL(KIND=DP) :: prev_cool_e

     prev_cool_n = 0.0_DP
     prev_cool_i = 0.0_DP
     prev_cool_e = 0.0_DP

     Cool_n = 0.0_DP
     Cool_i = 0.0_DP
     Cool_e = 0.0_DP

     IF (cool_KN==2) THEN
        ! Skip this entirely for analytical cooling, which will be computed
        ! in compute_molecular_cooling
        RETURN
     ENDIF

     ! ---------------------------
     ! Excitation of atomic hydrogen
     ! ---------------------------
     Dens_cool = Dens_H
     mass_cool = mass_H
     charge_cool = 0.0_DP

     CALL COLHAT (ntrhat,nlvhat,gsthat,elkhat,wlkhat,aijhat, &
               iuphat,jlohat,emihat,pop_hat)
     cooling_n_Hat = cool_n - prev_cool_n
     cooling_i_Hat = cool_i - prev_cool_i
     cooling_e_Hat = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of atomic carbon
     ! ---------------------------
     Dens_cool = Dens_C
     mass_cool = mass_C
     charge_cool = 0.0_DP

     CALL COLCAT (ntrcat,nlvcat,gstcat,elkcat,wlkcat,aijcat, &
               iupcat,jlocat,emicat,pop_cat)
     cooling_n_Cat = cool_n - prev_cool_n
     cooling_i_Cat = cool_i - prev_cool_i
     cooling_e_Cat = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of atomic oxygen
     ! ---------------------------
     Dens_cool = Dens_O
     mass_cool = mass_O
     charge_cool = 0.0_DP

     CALL COLOAT (ntroat,nlvoat,gstoat,elkoat,wlkoat,aijoat, &
                  iupoat,jlooat,emioat,pop_oat)
     cooling_n_Oat = cool_n - prev_cool_n
     cooling_i_Oat = cool_i - prev_cool_i
     cooling_e_Oat = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of atomic nitrogen
     ! ---------------------------
     Dens_cool = Dens_N
     mass_cool = mass_N
     charge_cool = 0.0_DP

     CALL COLNAT (ntrnat,nlvnat,gstnat,elknat,wlknat,aijnat, &
                  iupnat,jlonat,eminat,pop_nat)
     cooling_n_Nat = cool_n - prev_cool_n
     cooling_i_Nat = cool_i - prev_cool_i
     cooling_e_Nat = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of atomic sulfur
     ! ---------------------------
     Dens_cool = Dens_S
     mass_cool = mass_S
     charge_cool = 0.0_DP

     CALL COLSAT (ntrsat,nlvsat,gstsat,elksat,wlksat,aijsat, &
                  iupsat,jlosat,emisat,pop_sat)
     cooling_n_Sat = cool_n - prev_cool_n
     cooling_i_Sat = cool_i - prev_cool_i
     cooling_e_Sat = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of atomic silicium
     ! ---------------------------
     Dens_cool = Dens_Si
     mass_cool = mass_Si
     charge_cool = 0.0_DP

     CALL COLSiAT (ntrsiat,nlvsiat,gstsiat,elksiat,wlksiat,aijsiat, &
                   iupsiat,jlosiat,emisiat,pop_siat)
     cooling_n_Siat = cool_n - prev_cool_n
     cooling_i_Siat = cool_i - prev_cool_i
     cooling_e_Siat = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of ionised carbon
     ! ---------------------------
     Dens_cool = Dens_Cplus
     mass_cool = mass_Cplus
     charge_cool = 1.0_DP

     CALL COLCPL (ntrcpl,nlvcpl,gstcpl,elkcpl,wlkcpl,aijcpl, &
                  iupcpl,jlocpl,emicpl,pop_cpl)
     cooling_n_Cp = cool_n - prev_cool_n
     cooling_i_Cp = cool_i - prev_cool_i
     cooling_e_Cp = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of ionised nitrogen
     ! ---------------------------
     Dens_cool = Dens_Nplus
     mass_cool = mass_Nplus
     charge_cool = 1.0_DP

     CALL COLNPL (ntrnpl,nlvnpl,gstnpl,elknpl,wlknpl,aijnpl, &
                  iupnpl,jlonpl,eminpl,pop_npl)
     cooling_n_Np = cool_n - prev_cool_n
     cooling_i_Np = cool_i - prev_cool_i
     cooling_e_Np = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of ionised oxygen
     ! ---------------------------
     Dens_cool = Dens_Oplus
     mass_cool = mass_Oplus
     charge_cool = 1.0_DP

     CALL COLOPL (ntropl,nlvopl,gstopl,elkopl,wlkopl,aijopl, &
                  iupopl,jloopl,emiopl,pop_opl)
     cooling_n_Op = cool_n - prev_cool_n
     cooling_i_Op = cool_i - prev_cool_i
     cooling_e_Op = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of ionised sulfur
     ! ---------------------------
     Dens_cool = Dens_Splus
     mass_cool = mass_Splus
     charge_cool = 1.0_DP

     CALL COLSPL (ntrspl,nlvspl,gstspl,elkspl,wlkspl,aijspl, &
                  iupspl,jlospl,emispl,pop_spl)
     cooling_n_Sp = cool_n - prev_cool_n
     cooling_i_Sp = cool_i - prev_cool_i
     cooling_e_Sp = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

     ! ---------------------------
     ! Excitation of ionised silicium
     ! ---------------------------
     Dens_cool = Dens_Siplus
     mass_cool = mass_Siplus
     charge_cool = 1.0_DP

     CALL COLSiPL (ntrsipl,nlvsipl,gstsipl,elksipl,wlksipl,aijsipl, &
                   iupsipl,jlosipl,emisipl,pop_sipl)
     cooling_n_Sip = cool_n - prev_cool_n
     cooling_i_Sip = cool_i - prev_cool_i
     cooling_e_Sip = cool_e - prev_cool_e
     prev_cool_n = cool_n
     prev_cool_i = cool_i
     prev_cool_e = cool_e

    ! Commented by Pierre L. to avoid stability problems. 22/4/2010
    ! ---------------------------
    ! Excitation of ionized iron
    ! ---------------------------
    !!$
    !!$  Dens_cool = Dens_Feplus
    !!$  mass_cool = mass_Feplus
    !!$  charge_cool = 1.0_DP
    !!$
    !!$  call COLFEPL (ntrfepl,nlvfepl,gstfepl,elkfepl,wlkfepl,aijfepl, &
    !!$                iupfepl,jlofepl,emifepl,pop_fepl)
    !!$

    !-------------------------------------------------------
    ! Save for outputs- Tabone 05/2019
    ! Convention : line_rad_... are signed quantities
    !              < 0 means cooling
    !              > 0 means heating
    !-------------------------------------------------------
    line_rad_n_Hat   = cooling_n_Hat
    line_rad_n_Oat   = cooling_n_Oat
    line_rad_n_Op    = cooling_n_Op
    line_rad_n_Cat   = cooling_n_Cat
    line_rad_n_Cp    = cooling_n_Cp
    line_rad_n_Nat   = cooling_n_Nat
    line_rad_n_Np    = cooling_n_Np
    line_rad_n_Sat   = cooling_n_Sat
    line_rad_n_Sp    = cooling_n_Sp
    line_rad_n_Siat  = cooling_n_Siat
    line_rad_n_Sip   = cooling_n_Sip

    line_rad_i_Hat   = cooling_i_Hat
    line_rad_i_Oat   = cooling_i_Oat
    line_rad_i_Op    = cooling_i_Op
    line_rad_i_Cat   = cooling_i_Cat
    line_rad_i_Cp    = cooling_i_Cp
    line_rad_i_Nat   = cooling_i_Nat
    line_rad_i_Np    = cooling_i_Np
    line_rad_i_Sat   = cooling_i_Sat
    line_rad_i_Sp    = cooling_i_Sp
    line_rad_i_Siat  = cooling_i_Siat
    line_rad_i_Sip   = cooling_i_Sip

    line_rad_e_Hat   = cooling_e_Hat
    line_rad_e_Oat   = cooling_e_Oat
    line_rad_e_Op    = cooling_e_Op
    line_rad_e_Cat   = cooling_e_Cat
    line_rad_e_Cp    = cooling_e_Cp
    line_rad_e_Nat   = cooling_e_Nat
    line_rad_e_Np    = cooling_e_Np
    line_rad_e_Sat   = cooling_e_Sat
    line_rad_e_Sp    = cooling_e_Sp
    line_rad_e_Siat  = cooling_e_Siat
    line_rad_e_Sip   = cooling_e_Sip

    line_rad_n_atoms = line_rad_n_Hat  &
                     + line_rad_n_Oat  &
                     + line_rad_n_Op   &
                     + line_rad_n_Cat  &
                     + line_rad_n_Cp   &
                     + line_rad_n_Nat  &
                     + line_rad_n_Np   &
                     + line_rad_n_Sat  &
                     + line_rad_n_Sp   &
                     + line_rad_n_Siat &
                     + line_rad_n_Sip
    line_rad_i_atoms = line_rad_i_Hat  &
                     + line_rad_i_Oat  &
                     + line_rad_i_Op   &
                     + line_rad_i_Cat  &
                     + line_rad_i_Cp   &
                     + line_rad_i_Nat  &
                     + line_rad_i_Np   &
                     + line_rad_i_Sat  &
                     + line_rad_i_Sp   &
                     + line_rad_i_Siat &
                     + line_rad_i_Sip
    line_rad_e_atoms = line_rad_e_Hat  &
                     + line_rad_e_Oat  &
                     + line_rad_e_Op   &
                     + line_rad_e_Cat  &
                     + line_rad_e_Cp   &
                     + line_rad_e_Nat  &
                     + line_rad_e_Np   &
                     + line_rad_e_Sat  &
                     + line_rad_e_Sp   &
                     + line_rad_e_Siat &
                     + line_rad_e_Sip

  END SUBROUTINE LINE_THERMAL_BALANCE

! %%%%%%%%%%%%%%%%%%%%
  subroutine ALPHATEFF
! %%%%%%%%%%%%%%%%%%%%

! Nov. 2002: initializes automatically the alpha's and Teff's
!            according to the charge and mass of the coolant

  USE MODULE_CHEMICAL_SPECIES
  USE MODULE_PHYS_VAR
  USE MODULE_CONSTANTS

  implicit none

  REAL (kind=DP) :: T1, T2     ! Partial effective temperatures

! Neutral atom

  if (charge_cool == 0.0_DP) then

    alp_H  = 1.0_DP
    alp_H2 = 1.0_DP
    alp_He = 1.0_DP
    alp_Hp = mass_Hplus / (mass_Hplus + mass_cool)

!--- effective temperatures with neutrals : Tn

    Teff_H  = Tn
    Teff_H2 = Tn
    Teff_He = Tn

!--- effective temperatures for H+-neutral:

    T1 = ABS_DeltaV * ABS_DeltaV * mass_Hplus * mass_cool &
     / (3.0_DP * kB * (mass_Hplus + mass_cool))
    T2 = (mass_cool * Tn + mass_Hplus * Ti) / (mass_Hplus + mass_cool)
    Teff_Hp = T1 + T2

  else

! Positive Ion

    alp_H  = mass_cool / (mass_H + mass_cool)
    alp_H2 = mass_cool / (mass_H2 + mass_cool)
    alp_He = mass_cool / (mass_He + mass_cool)
    alp_Hp = 0.0_DP

!--- effective temperatures with neutrals :

    T1 = ABS_DeltaV * ABS_DeltaV * mass_H * mass_cool &
     / (3.0_DP * kB * (mass_H + mass_cool))
    T2 = (mass_cool * Ti + mass_H * Tn) / (mass_H + mass_cool)
    Teff_H  = T1 + T2
    T1 = ABS_DeltaV * ABS_DeltaV * mass_H2 * mass_cool &
     / (3.0_DP * kB * (mass_H2 + mass_cool))
    T2 = (mass_cool * Ti + mass_H2 * Tn) / (mass_H2 + mass_cool)
    Teff_H2 = T1 + T2
    T1 = ABS_DeltaV * ABS_DeltaV * mass_He * mass_cool &
     / (3.0_DP * kB * (mass_He + mass_cool))
    T2 = (mass_cool * Ti + mass_He * Tn) / (mass_He + mass_cool)
    Teff_He = T1 + T2

!--- effective temperatures for H+-ion : Ti

    Teff_Hp = Ti

  end if

!--- effective temperature for electron collisions : Te

  Teff_e = Te

  end subroutine ALPHATEFF

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine THERMAL_LOSS (gst,elk,wlk,aij,iup,jlo,ec_H,ec_orthoH2, &
                           ec_paraH2,ec_He,ec_Hp,ec_e,emi,pop,ntr,nlv)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Nov 2002:
! * finish calculating statistical equilibrium matrix
!   using de-excitation rates computed in subroutine COLxxx
!       - multiply by collisioner density
!       - Derive excitation rates from detailed balance
!       - compute diagonal terms
!       - add radiative transitions
! * solve for level populations using LAPACK routine dgesvx (Pierre L. 18/3/2010)
! * Compute radiative emissivities in each line => emi
! * compute thermal losses of the 3 fluids due to
!   collisional excitation of these lines => Cool_n, Cool_i, Cool_e

  USE MODULE_CONSTANTS
  USE MODULE_CHEMICAL_SPECIES
  USE MODULE_H2
  use tex
  implicit none

  INTEGER       , intent(in)                     :: ntr, nlv
  INTEGER       , intent(in),  DIMENSION (ntr)   :: iup, jlo
  REAL (kind=DP), intent(in),  DIMENSION (ntr)   :: wlk, aij
  REAL (kind=DP), intent(in),  DIMENSION (nlv)   :: gst, elk
  REAL (kind=DP), intent(out), DIMENSION (ntr)   :: emi

  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H, ec_He,ec_Hp, ec_e
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_orthoH2, ec_paraH2
  REAL (kind=DP), DIMENSION (nlv,nlv) :: ec_H2, aa ! total H2, total matrix
  REAL (kind=DP), DIMENSION (nlv)     :: pop       ! fractional level populations

!  Used by LAPACK

  INTEGER                            :: info
   INTEGER                            :: nrhs = 1
   REAL (kind=DP), dimension (nlv)    :: indx

  INTEGER        :: i, j, itr
  REAL (kind=DP) :: radi, boltz, lasum
  REAL (kind=DP) :: bth_H, bth_H2, bth_He, bth_Hp, bth_e

! Warning for radiative losses
  logical,save::warning=.false.

  bth_H  = 0.0_DP
  bth_H2 = 0.0_DP
  bth_He = 0.0_DP
  bth_Hp = 0.0_DP
  bth_e  = 0.0_DP

! finish calculating the collisional excitation matrix:
!       - multiply by collisioner density
! 	- Derive excitation rates from detailed balance
!       - compute diagonal terms
!       - add radiative transitions

  ec_H = ec_H * Dens_H
  ec_He = ec_He * Dens_He
  ec_Hp = ec_Hp * Dens_Hplus
  ec_e  = ec_e * Dens_e
  ec_paraH2  = ec_paraH2 * Dens_paraH2
  ec_orthoH2 = ec_orthoH2 * Dens_orthoH2

  ec_H2 = ec_paraH2 + ec_orthoH2

  do i=2,nlv
    do j=i-1,1,-1
      boltz = texp(-(elk(i)-elk(j)) / Teff_H) * gst(i) / gst(j)
      ec_H(i,j)  = ec_H(j,i)  * boltz
      boltz = texp(-(elk(i)-elk(j)) / Teff_H2) * gst(i) / gst(j)
      ec_H2(i,j) = ec_H2(j,i) * boltz
      boltz = texp(-(elk(i)-elk(j)) / Teff_He) * gst(i) / gst(j)
      ec_He(i,j) = ec_He(j,i) * boltz
      boltz = texp(-(elk(i)-elk(j)) / Teff_Hp) * gst(i) / gst(j)
      ec_Hp(i,j) = ec_Hp(j,i) * boltz
      boltz = texp(-(elk(i)-elk(j)) / Teff_e) * gst(i) / gst(j)
      ec_e(i,j) = ec_e(j,i) * boltz
    end do
  end do

  do i=1,nlv
    ec_H(i,i)  = - SUM(ec_H(:,i))
    ec_H2(i,i) = - SUM(ec_H2(:,i))
    ec_He(i,i) = - SUM(ec_He(:,i))
    ec_Hp(i,i) = - SUM(ec_Hp(:,i))
    ec_e(i,i)  = - SUM(ec_e(:,i))
  end do

! Total collision matrix

  aa = ec_H + ec_H2 + ec_He + ec_Hp + ec_e

! add radiative decay terms

  do itr=1,ntr
    i = iup(itr)
    j = jlo(itr)
    aa(i,i) = aa(i,i) - aij(itr)
    aa(j,i) = aa(j,i) + aij(itr)
  end do

! Replace last row of matrix by conservation, and set right hand side

  aa(nlv,:) = 1.0_DP
  pop = 0.0_DP
  pop(nlv) = 1.0_DP

!!$! Solve with LAPACK routine dgesvx (wrapper by Pierre L. 18/3/2010)
!!$  call square_solve(aa,nlv,pop)
! Solve with LAPACK routine dgesv
  call dgesv (nlv, nrhs, aa, nlv, indx, pop, nlv, info)

  lasum = 0.0_DP
  do i=1,nlv
! Pierre L. comments the following line:
    pop(i) = max(pop(i),0.0_DP)
! (Otherwise, populations get unstable and go wrong. 18/3/2010)
    lasum = lasum + max(pop(i),0.0_DP)
  end do
  if (lasum.ne.0d0) pop = pop / lasum
  ! Pierre L. adds this warning:
  do i=1,nlv
    if (pop(I)<-1d-1) then
       print*,'pop(',i,')<0:',pop(i)
    endif
 enddo
! print *, '    level energies (K):', elk
! print *, '    fractional populations:', pop

! calculate net thermal energy contributed by each collisioner in excitation

  do j=1,nlv-1
    do i=j+1,nlv
      bth_H  = bth_H  + (ec_H(j,i)  * pop(i) - ec_H(i,j)  * pop(j)) * (elk(i) - elk(j))
      bth_H2 = bth_H2 + (ec_H2(j,i) * pop(i) - ec_H2(i,j) * pop(j)) * (elk(i) - elk(j))
      bth_He = bth_He + (ec_He(j,i) * pop(i) - ec_He(i,j) * pop(j)) * (elk(i) - elk(j))
      bth_Hp = bth_Hp + (ec_Hp(j,i) * pop(i) - ec_Hp(i,j) * pop(j)) * (elk(i) - elk(j))
      bth_e  = bth_e  + (ec_e(j,i)  * pop(i) - ec_e(i,j)  * pop(j)) * (elk(i) - elk(j))
    end do
  end do


! calculate total radiative losses

  radi = 0.0_DP
  do i=1,ntr
    emi(i) = wlk(i) * aij(i) * pop(iup(i))
    radi = radi + emi(i)
  end do

!  print *, '  total radiative losses =', radi

! Check that total thermal loss = total radiative emission within 1e-1

  if (abs(bth_H+bth_H2+bth_He+bth_Hp+bth_e+radi) > 1.0e-1*abs(radi).and..not.warning) then
     print*,'WARNING: Thermal losses are rather not like radiative emission..'
     print *, bth_H, bth_H2, bth_He, bth_Hp, bth_e
     print *, bth_H+bth_H2+bth_He+bth_Hp+bth_e, radi, mass_cool/amu 
     warning=.true.
     print*,'This warning is now disabled'
!!$     print*,'But I prefer to stop...'
!!$     stop
  endif

  Cool_n = Cool_n + kB * (bth_H * alp_H + bth_H2 * alp_H2 &
                  + bth_He * alp_He + bth_Hp * alp_Hp) * Dens_cool
  Cool_i = Cool_i + kB * (bth_H * (1.0_DP-alp_H) + bth_H2 * (1.0_DP-alp_H2) &
                  + bth_He * (1.0_DP-alp_He) + bth_Hp * (1.0_DP-alp_Hp)) * Dens_cool
  Cool_e = Cool_e + kB * bth_e * Dens_cool

! Nov 2002: multiply emi by Dens_cool to get erg/s/cm3 (not per atom)

  emi = emi * kB * Dens_cool
  end subroutine THERMAL_LOSS

! %%%%%%%%%%%%%%%%%%%%%
  subroutine LINE_INTEG
! %%%%%%%%%%%%%%%%%%%%%

  USE MODULE_CONSTANTS
  USE MODULE_PHYS_VAR
  USE MODULE_CHEMICAL_SPECIES
  implicit none

  inthat = inthat + (emihat_o + emihat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intcat = intcat + (emicat_o + emicat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intnat = intnat + (eminat_o + eminat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intoat = intoat + (emioat_o + emioat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intsat = intsat + (emisat_o + emisat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intsiat = intsiat + (emisiat_o + emisiat) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intcpl = intcpl + (emicpl_o + emicpl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intnpl = intnpl + (eminpl_o + eminpl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intopl = intopl + (emiopl_o + emiopl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intspl = intspl + (emispl_o + emispl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intsipl = intsipl + (emisipl_o + emisipl) * &
         0.5_DP * dist_step / (4.0_DP*pi)
  intfepl = intfepl + (emifepl_o + emifepl) * &
         0.5_DP * dist_step / (4._DP*pi) 

  emihat_o  = emihat
  emicat_o  = emicat
  eminat_o  = eminat
  emioat_o  = emioat
  emisat_o  = emisat
  emisiat_o = emisiat
  emicpl_o  = emicpl
  eminpl_o  = eminpl
  emiopl_o  = emiopl
  emispl_o  = emispl
  emisipl_o = emisipl
  emifepl_o = emifepl

  end subroutine LINE_INTEG

END MODULE MODULE_LINE_EXCIT
