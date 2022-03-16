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

MODULE MODULE_MOLECULAR_COOLING
  !****************************************************************************
  !** The module 'MODULE_MOLECULAR_COOLING contains all variables and        **
  !** subroutines needed to calculate the radiative cooling rates due to     **
  !** several molecules, but not H2 (see MODULE_H2).                         **
  !**                                                                        **
  !** The cooling from H2O and CO is taken from Neufeld & Kaufman 1993.      **
  !**                                                                        **
  !** note : the only public routine is MOLECULAR_COOLING                    **
  !****************************************************************************
  USE MODULE_TECHCONFIG
  USE MODULE_VAR_TH_BALANCE
  
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !----------------------------------------------------------
  ! cooling rate (erg/cm3/s) calculated in MOLECULAR COOLING
  !----------------------------------------------------------
  REAL(KIND=DP) :: cooling_13CO      = 0.0_DP  ! cooling rate for 13CO
  REAL(KIND=DP) :: cooling_vib_13CO  = 0.0_DP  ! cooling rate for vibrationa 13CO added by Tabone 02/2018
  REAL(KIND=DP) :: cooling_rot_13CO  = 0.0_DP  ! cooling rate for vibrationa 13CO added by Tabone 02/2018
  REAL(KIND=DP) :: cooling_OH        = 0.0_DP  ! cooling rate for OH
  REAL(KIND=DP) :: cooling_NH3       = 0.0_DP  ! cooling rate for NH3

  REAL(KIND=DP) :: cooling_rot_o_H2O = 0.0_DP  ! cooling rate for rotational ortho-H2O
  REAL(KIND=DP) :: cooling_rot_p_H2O = 0.0_DP  ! cooling rate for rotational para-H2O
  REAL(KIND=DP) :: cooling_vib_H2O   = 0.0_DP  ! cooling rate for vibrational H2O
  REAL(KIND=DP) :: cooling_H2O       = 0.0_DP  ! total cooling rate for H2O

  REAL(KIND=DP) :: cooling_rot_CO    = 0.0_DP  ! cooling rate for rotational CO
  REAL(KIND=DP) :: cooling_vib_CO    = 0.0_DP  ! cooling rate for vibrational CO
  REAL(KIND=DP) :: cooling_CO        = 0.0_DP  ! total cooling rate for CO

  REAL(KIND=DP) :: B_photo_grain_n   = 0.0_DP  ! Photo-electric heating

  REAL(KIND=DP) :: molec_cool        = 0.0_DP  ! total molecular (not H2) cooling rate

  INTEGER(KIND=LONG)         :: f_err_cool                  ! file number for error messages
  CHARACTER(len=lenfilename) :: n_err_cool = 'err_cool.out' ! file name for error messages
  LOGICAL                    :: err_cool   = .FALSE.        ! tells if there was errors

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity
  PRIVATE :: ROTATIONAL_CO, VIBRATIONAL_CO
  PRIVATE :: ROTATIONAL_o_H2O, ROTATIONAL_p_H2O, VIBRATIONAL_H2O


CONTAINS


  SUBROUTINE COMPUTE_MOLECULAR_COOLING
    !---------------------------------------------------------------------------
    ! called by :
    !     SOURCE
    ! purpose :
    !     Calculate the cooling rates (erg/cm3/s) due to 13CO, H2O, OH, NH3.
    ! subroutine/function needed :
    !     ROTATIONAL_CO
    !     VIBRATIONAL_CO
    !     ROTATIONAL_o_H2O
    !     ROTATIONAL_p_H2O
    !     VIBRATIONAL_H2O
    ! input variables :
    ! output variables :
    ! results :
    !     cooling_13CO
    !     cooling_OH
    !     cooling_NH3
    !     cooling_rot_CO, cooling_vib_CO, cooling_CO
    !     cooling_rot_o_H2O, cooling_rot_p_H2O, cooling_vib_H2O, cooling_H2O
    !     molec_cool
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY         : Tn, Vgrad, Cool_KN,nH
    USE MODULE_CHEMICAL_SPECIES, ONLY : Dens_H, Dens_H2, Dens_CO, &
                                        Dens_H2O, Dens_OH, Dens_NH3, Dens_13CO ! Tabone 02/2018
    USE MODULE_CONSTANTS, ONLY        : Zero
    USE MODULE_DEBUG_JLB

    IMPLICIT NONE
    REAL(KIND=DP)            :: Dens_H_H2
    REAL(KIND=DP)            :: op_H2O
    REAL(KIND=DP)            :: Dens_ortho_H2O
    REAL(KIND=DP)            :: Dens_para_H2O
    REAL(KIND=DP), PARAMETER :: G = 1._DP      ! geometrical factor
    INTEGER(KIND=LONG)       :: error
    REAL(KIND=DP), PARAMETER :: Tmin = 0._DP

    cooling_rot_13CO  = Zero ! Tabone 02/2018
    cooling_vib_13CO  = Zero ! Tabone 02/2018
    cooling_13CO      = Zero
    cooling_OH        = Zero
    cooling_NH3       = Zero
    cooling_rot_CO    = Zero
    cooling_vib_CO    = Zero
    cooling_CO        = Zero
    cooling_rot_o_H2O = Zero
    cooling_rot_p_H2O = Zero
    cooling_vib_H2O   = Zero
    cooling_H2O       = Zero
    IF (cool_KN==2) THEN
       molec_cool=analytical_cooling(nH,Tn)
       RETURN
    ENDIF
    !---------------------------------------
    ! "simple" cooling for 13CO, OH and NH3 ! modified by Tabone to include a K&N calculation for 13CO
    !---------------------------------------
    Dens_H_H2 = Dens_H + Dens_H2 / sqrt(2.0_DP)

    ! H2O cooling is now calculated more precisely, see below
    ! cooling_H2O  = 1.0D-26 * Tn * EXP(-35.0_DP/Tn) * Dens_H2O * Dens_H_H2
    ! *** Warning Tabone *** - expressions below not valid at large densities ?
    cooling_OH   = 1.0D-26 * Tn * EXP(-120.0_DP/Tn) * Dens_OH &
                 * (10.0_DP * Dens_H + Dens_H2 / 1.414_DP)
    cooling_NH3  = 1.0D-26 * Tn * EXP(-40.0_DP/Tn) * Dens_NH3 * Dens_H_H2

    !----------------------------------------
    ! "more precise" cooling for H2O and CO
    ! taken from Neufeld & Kaufman 1993
    ! note : this cooling is calculated only
    !        for Tn > Tmin
    !----------------------------------------

    ! 31 mai 2001 - JLB : switch off (temporarily?) kaufmann & Neufeld cooling

    IF (Tn > Tmin .AND. Cool_KN == 1) THEN
       !--- CO ---
       ! rotational
       CALL ROTATIONAL_CO (Dens_H2, G*Dens_CO/Vgrad, Tn, cooling_rot_CO, error)
       IF (error /= 0) CALL W_ERR_COOL(error,'ROTATIONAL_CO')
       ! change to erg/cm3/s
       cooling_rot_CO = cooling_rot_CO * Dens_H2 * Dens_CO

       ! vibrational
       CALL VIBRATIONAL_CO (Dens_H2, G*Dens_CO/Vgrad, Tn, cooling_vib_CO, error)
       IF (error /= 0) CALL W_ERR_COOL(error,'VIBRATIONAL_CO')
       ! change to erg/cm3/s
       cooling_vib_CO = cooling_vib_CO * Dens_H2 * Dens_CO

       ! total
       cooling_CO = cooling_rot_CO + cooling_vib_CO


      !--- 13CO --- added by Tabone 02/2018 from MAMOJ
       Dens_13CO = Dens_CO / 90._DP

       ! rotational
       CALL ROTATIONAL_CO(Dens_H2,G*Dens_13CO/Vgrad,Tn,cooling_rot_13CO,error)
       IF (error /= 0) CALL W_ERR_COOL(error,'ROTATIONAL_CO')
       ! change to erg/cm3/s
       cooling_rot_13CO = cooling_rot_13CO * Dens_H2 * Dens_13CO

       ! vibrational
       CALL VIBRATIONAL_CO(Dens_H2,G*Dens_13CO/Vgrad,Tn,cooling_vib_13CO,error)
       IF (error /= 0) CALL W_ERR_COOL(error,'VIBRATIONAL_CO')
       ! change to erg/cm3/s
       cooling_vib_13CO = cooling_vib_13CO * Dens_H2 * Dens_13CO

       ! total
       cooling_13CO = cooling_rot_13CO + cooling_vib_13CO

       !--- H2O ---
       ! rotational (para and ortho)
       !ccccccccccccccccccccccccccccccccccccc
       !ccccccccccccccccccccccccccccccccccccc
       op_H2O = 3._DP ! ortho/para ratio for H2O
       !ccccccccccccccccccccccccccccccccccccc
       !ccccccccccccccccccccccccccccccccccccc
       Dens_para_H2O  = Dens_H2O / (1._DP + op_H2O)
       Dens_ortho_H2O = Dens_H2O * op_H2O / (1._DP + op_H2O)

       CALL ROTATIONAL_p_H2O (Dens_H2, G*Dens_para_H2O/Vgrad, Tn, cooling_rot_p_H2O, error)
       IF (error /= 0) CALL W_ERR_COOL(error,'ROTATIONAL_p_H2O')
       ! change to erg/cm3/s
       cooling_rot_p_H2O = cooling_rot_p_H2O * Dens_H2 * Dens_para_H2O

       CALL ROTATIONAL_o_H2O (Dens_H2, G*Dens_ortho_H2O/Vgrad, Tn, cooling_rot_o_H2O, error)
       IF (error /= 0) CALL W_ERR_COOL(error,'ROTATIONAL_o_H2O')
       ! change to erg/cm3/s
       cooling_rot_o_H2O = cooling_rot_o_H2O * Dens_H2 * Dens_ortho_H2O

       ! vibrational
       CALL VIBRATIONAL_H2O (Dens_H2, G*Dens_H2O/Vgrad, Tn, cooling_vib_H2O, error)
       IF (error /= 0) CALL W_ERR_COOL(error,'VIBRATIONAL_H2O')
       ! change to erg/cm3/s
       cooling_vib_H2O = cooling_vib_H2O * Dens_H2 * Dens_H2O

       ! total
       cooling_H2O = cooling_rot_o_H2O + cooling_rot_p_H2O + cooling_vib_H2O

    ELSE
       cooling_CO = 0.0_DP
       !---------- Tabone 01-05-2019: fitted formula valid for any density --------
       cooling_13CO =                               ( Dens_CO*5.5D-24*( Tn**(2.2_DP) / &
                         (1._DP+(Tn/1200._DP)**2.05_DP) )*EXP(-2.0_DP/Tn)/(90.0_DP)) * &
                                 ((Dens_H_H2/1.0D6)**4/((Dens_H_H2/1.0D6)**4+1._DP)) + &
                  ( 1.0D-26 * Tn * EXP(-5.0_DP /Tn) * Dens_CO * Dens_H_H2 / 90.0_DP) * &
                                                       EXP(-(Dens_H_H2/1.0D6)*3._DP)
       
       cooling_H2O =                     ( Dens_H2O*5.5D-14*( (Tn/300._DP)**(2.8_DP) / &
                                    (1._DP+(Tn/300._DP)**2.35_DP) )*EXP(-30._DP/Tn)) * &
                               ((Dens_H_H2/1.0D10)**4/((Dens_H_H2/1.0D10)**4+1._DP)) + & 
                            (1.0D-26 * Tn * EXP(-35.0_DP/Tn) * Dens_H2O * Dens_H_H2) * & 
                                                      EXP(-(Dens_H_H2/1.0D10)*3._DP) 
       !---- Former low density H2O and 13CO cooling ----
       ! cooling_H2O = 1.0D-26 * Tn * EXP(-35.0_DP/Tn) * Dens_H2O * Dens_H_H2
       ! cooling_13CO = 1.0D-26 * Tn * EXP(-5.0_DP /Tn) * Dens_CO * Dens_H_H2 / 90.0_DP
    END IF

    !------------------------------------
    ! sum of the different contributions
    !------------------------------------

    molec_cool =  cooling_13CO &
               +  cooling_OH   &
               +  cooling_NH3  &
               +  cooling_CO   &
               +  cooling_H2O

    !---------------------------------------
    !--------Save in output variables-------
    !---------------------------------------
    line_rad_n_OH    = -cooling_OH
    line_rad_n_NH3   = -cooling_NH3
    line_rad_n_13CO  = -cooling_13CO      ! sum over rot + vib
    
    line_rad_n_rCO   = -cooling_rot_CO  
    line_rad_n_vCO   = -cooling_vib_CO 
    line_rad_n_CO    = -cooling_CO        ! sum over rot + vib
    
    line_rad_n_roH2O = -cooling_rot_o_H2O
    line_rad_n_rpH2O = -cooling_rot_p_H2O
    line_rad_n_vH2O  = -cooling_vib_H2O
    line_rad_n_H2O   = -cooling_H2O       ! sum over vib + rot_o + rot_p
 
    line_rad_n_molec = line_rad_n_13CO &
                     + line_rad_n_OH   &
                     + line_rad_n_NH3  &
                     + line_rad_n_CO   &
                     + line_rad_n_H2O  &
                     + line_rad_n_H2      ! add cooling by H2 saved in SOURCE (evolution.f90)
    line_rad_i_molec = line_rad_i_H2      ! term computed in SOURCE (evolution.f90)
    line_rad_e_molec = line_rad_e_H2   &
                     + line_rad_e_CO      ! terms computed in SOURCE (evolution.f90)

    line_rad_n   = line_rad_n_atoms + line_rad_n_molec
    line_rad_i   = line_rad_i_atoms + line_rad_i_molec
    line_rad_e   = line_rad_e_atoms + line_rad_e_molec
    line_rad_tot = line_rad_n + line_rad_i + line_rad_e

  END SUBROUTINE COMPUTE_MOLECULAR_COOLING




  SUBROUTINE ROTATIONAL_CO(DENS,COL,TEMP,COOL,IER)
    !---------------------------------------------------------------------------
    ! called by :
    !     COMPUTE_MOLECULAR_COOLING
    ! purpose :
    !     Calculate the cooling rate (erg.cm3.s-1) due to rotational CO. Adapted
    !     from Neufeld & Kaufman 1993.
    !     Warning : to obtain the cooling rate in erg/cm3/s, multiply by the
    !               number density (cm-3) of the two reactants (H2 and CO).
    ! subroutine/function needed :
    !     SPLIE2 (from Numerical Recipies in F90)
    !     SPLIN2 (from Numerical Recipies in F90)
    ! input variables :
    !     DENS -> H2 number density (cm-3)
    !     COL  -> optical depth parameter = G * n(CO)/Vgrad (cm-2.km-1.s)
    !             where G is the geometrical factor from Neufeld & Melnick 1991
    !     TEMP -> kinetic temperature of neutrals (K)
    ! output variables :
    !     COOL -> cooling rate (erg.cm3.s-1)
    !     IER  -> warning flag :
    !                0  => NORMAL EXECUTION
    !                1  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !               11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO LARGE, > 2000)
    !              -11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO SMALL, < 10)
    !               10  => TEMP OUT OF RANGE    (TOO LARGE, > 2000)
    !              -10  => TEMP OUT OF RANGE    (TOO SMALL, < 10)
    !
    !             WHEN IER .NE. 0 THE VALUE RETURNED FOR COOL IS
    !             AN EXTRAPOLATION OF DUBIOUS ACCURACY
    ! results :
    !---------------------------------------------------------------------------
    USE NUM_REC, ONLY          : SPLIE2, SPLIN2
    USE MODULE_CONSTANTS, ONLY : Zero

    IMPLICIT NONE

    REAL(KIND=DP),      INTENT(in)           :: DENS
    REAL(KIND=DP),      INTENT(in)           :: COL
    REAL(KIND=DP),      INTENT(in)           :: TEMP
    REAL(KIND=DP),      INTENT(out)          :: COOL
    INTEGER(KIND=LONG), INTENT(out)          :: IER

    INTEGER(KIND=LONG), PARAMETER            :: M = 11
    INTEGER(KIND=LONG), PARAMETER            :: N = 10
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: TKIN
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: CMAX
    REAL(KIND=DP),      SAVE, DIMENSION(N)   :: COLUMN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: F
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: DENSHF
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: COOLFN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: ABEST
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: Y2A
    INTEGER(KIND=LONG)                       :: I
    INTEGER(KIND=LONG)                       :: J
    INTEGER(KIND=LONG), SAVE                 :: Ncall = 0
    REAL(KIND=DP)                            :: DH
    REAL(KIND=DP)                            :: DCR
    REAL(KIND=DP)                            :: FACTOR
    REAL(KIND=DP)                            :: T
    REAL(KIND=DP)                            :: C
    REAL(KIND=DP)                            :: TMINUS
    REAL(KIND=DP)                            :: TPLUS
    REAL(KIND=DP)                            :: CPLUS
    REAL(KIND=DP)                            :: T1
    REAL(KIND=DP)                            :: T2
    REAL(KIND=DP)                            :: COOL1
    REAL(KIND=DP)                            :: COOL2
    REAL(KIND=DP)                            :: DFDT
    REAL(KIND=DP)                            :: C1
    REAL(KIND=DP)                            :: C2
    REAL(KIND=DP)                            :: DFDC

    Ncall = Ncall + 1_LONG
    ! initialize the arrays only on first call to subroutine
    IF (Ncall == 1_LONG) THEN
       TKIN = &
            (/1.0000, 1.3010, 1.4771, 1.6990, 1.9031, 2.0000 &
            , 2.4771, 2.7782, 3.0000, 3.1761, 3.3010 /)
       COLUMN = &
            (/14.5000,15.0000,15.5000,16.0000,16.5000,17.0000 &
            , 17.5000,18.0000,18.5000,19.0000 /)
       CMAX = &
            (/24.9553,24.5569,24.3560,24.1165,23.8846,23.8013 &
            , 23.4032,23.0729,22.8126,22.6073,22.4738 /)
       F = RESHAPE( &
            (/21.1462,20.3891,19.9734,19.4701,19.0225,18.8142 &
            , 17.8243,17.2328,16.8220,16.5203,16.3264,21.2500 &
            , 20.4489,20.0154,19.4963,19.0392,18.8277,17.8287 &
            , 17.2350,16.8233,16.5211,16.3270,21.4471,20.5821 &
            , 20.1170,19.5655,19.0861,18.8663,17.8425,17.2419 &
            , 16.8273,16.5237,16.3289,21.7345,20.8077,20.3065 &
            , 19.7114,19.1956,18.9604,17.8815,17.2624,16.8397 &
            , 16.5318,16.3349,22.0842,21.1097,20.5793,19.9452 &
            , 19.3916,19.1382,17.9762,17.3180,16.8751,16.5558 &
            , 16.3529,22.4730,21.4622,20.9103,20.2476,19.6653 &
            , 19.3972,18.1535,17.4413,16.9626,16.6196,16.4028 &
            , 22.8873,21.8470,21.2788,20.5953,19.9926,19.7142 &
            , 18.4100,17.6494,17.1310,16.7558,16.5173,23.3163 &
            , 22.2534,21.6722,20.9727,20.3550,20.0691,18.7232 &
            , 17.9287,17.3799,16.9773,16.7188,23.7607,22.6752 &
            , 22.0830,21.3704,20.7410,20.4496,19.0737,18.2564 &
            , 17.6877,17.2677,16.9979,24.2121,23.1083,22.5066 &
            , 21.7830,21.1441,20.8482,19.4497,18.6163,18.0346 &
            , 17.6052,17.3338 /),(/M,N/))
       DENSHF = RESHAPE( &
            (/ 3.3116, 3.7168, 3.9577, 4.2410, 4.3926, 4.5149 &
            ,  4.9462, 5.2445, 5.4266, 5.5730, 5.6778, 3.1147 &
            ,  3.5864, 3.8652, 4.1653, 4.3542, 4.4671, 4.9263 &
            ,  5.2340, 5.4206, 5.5692, 5.6750, 2.7964, 3.3351 &
            ,  3.6439, 4.0065, 4.2329, 4.3634, 4.8831, 5.2030 &
            ,  5.4168, 5.5572, 5.6815, 2.3659, 2.9467, 3.3111 &
            ,  3.7333, 3.9902, 4.1308, 4.7536, 5.1330, 5.3607 &
            ,  5.5354, 5.6526, 1.8982, 2.5016, 2.8742, 3.3324 &
            ,  3.6089, 3.7829, 4.5023, 4.9539, 5.2388, 5.4361 &
            ,  5.5755, 1.3965, 2.0233, 2.4039, 2.8620, 3.1680 &
            ,  3.3420, 4.1179, 4.6470, 4.9971, 5.2489, 5.4130 &
            ,  0.9008, 1.5183, 1.9023, 2.3783, 2.6913, 2.8668 &
            ,  3.6742, 4.2478, 4.6324, 4.9184, 5.1323, 0.4036 &
            ,  1.0215, 1.4063, 1.8856, 2.1867, 2.3765, 3.1962 &
            ,  3.7738, 4.1802, 4.5055, 4.7296,-0.0941, 0.5236 &
            ,  0.9087, 1.3900, 1.6902, 1.8811, 2.6930, 3.2881 &
            ,  3.7032, 4.0255, 4.2613,-0.5915, 0.0258, 0.4109 &
            ,  0.8786, 1.1926, 1.3843, 2.1962, 2.7944, 3.2014 &
            ,  3.5375, 3.7765 /),(/M,N/))
       ABEST = RESHAPE( &
            (/ 0.4220, 0.3880, 0.3780, 0.3890, 0.4180, 0.4310 &
            ,  0.3780, 0.3750, 0.3690, 0.3590, 0.3460, 0.4130 &
            ,  0.3770, 0.3660, 0.3790, 0.4080, 0.4220, 0.3720 &
            ,  0.3710, 0.3670, 0.3570, 0.3450, 0.4180, 0.3750 &
            ,  0.3620, 0.3710, 0.3970, 0.4100, 0.3610, 0.3630 &
            ,  0.3600, 0.3520, 0.3410, 0.4360, 0.3860, 0.3720 &
            ,  0.3750, 0.3940, 0.4050, 0.3480, 0.3500, 0.3490 &
            ,  0.3420, 0.3320, 0.4600, 0.4260, 0.4090, 0.4050 &
            ,  0.4060, 0.4120, 0.3510, 0.3390, 0.3350, 0.3290 &
            ,  0.3190, 0.4900, 0.4650, 0.4500, 0.4440, 0.4480 &
            ,  0.4520, 0.3680, 0.3450, 0.3290, 0.3190, 0.3080 &
            ,  0.5210, 0.4950, 0.4830, 0.4840, 0.4930, 0.4900 &
            ,  0.3980, 0.3730, 0.3500, 0.3220, 0.3070, 0.5480 &
            ,  0.5250, 0.5180, 0.5190, 0.5250, 0.5210, 0.4300 &
            ,  0.4020, 0.3780, 0.3540, 0.3260, 0.5660, 0.5480 &
            ,  0.5450, 0.5540, 0.5480, 0.5440, 0.4580, 0.4330 &
            ,  0.4100, 0.3830, 0.3570, 0.5830, 0.5680, 0.5650 &
            ,  0.5750, 0.5660, 0.5610, 0.4810, 0.4620, 0.4380 &
            ,  0.4140, 0.3870 /),(/M,N/))
    END IF

    DO I = 1, M
       DO J = 1, N
          DH = 10._DP**DENSHF(I,J)
          DCR = 10._DP**(CMAX(I) - F(I,J))
          FACTOR = 1._DP + (1._DP - DH/DCR) * (DENS / DH)**ABEST(I,J) + DENS / DCR
          COOLFN(I,J) = LOG10(FACTOR) + CMAX(I)
       END DO
    END DO

    CALL SPLIE2(TKIN,COLUMN,COOLFN,Y2A)
    ! [2-D SPLINE FIT FROM NUMERICAL RECIPES]

    T = LOG10(TEMP + 1._DP)
    C = LOG10(COL + 1._DP)
    IF (C <= COLUMN(1)) C = COLUMN(1)

    TMINUS = TKIN(1) - T
    TPLUS  = T - TKIN(M)
    CPLUS  = C - COLUMN(N)

    IF (CPLUS <= Zero) THEN
       !--------------------------
       !A) COLUMN is within range
       !--------------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------------
          !A.1) TEMPERATURE > TMAX, COLUMN WITHIN RANGE
          !---------------------------------------------
          T1 = TKIN(M-1)
          T2 = TKIN(M)
          COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
          COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
          DFDT = (COOL2 - COOL1) / (T2 - T1)
          COOL = COOL2 + TPLUS * DFDT
          IER = 10
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------------
             !A.2) TEMPERATURE < TMIN, COLUMN WITHIN RANGE
             !---------------------------------------------
             T1 = TKIN(1)
             T2 = TKIN(2)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
             DFDT = (COOL2 - COOL1) / (T2 - T1)
             COOL = COOL1 - TMINUS * DFDT
             IER = -10
          ELSE
             !-----------------------------------------
             !A.3) TEMPERATURE AND COLUMN WITHIN RANGE
             !-----------------------------------------
             COOL = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C)
             IER = 0
          END IF
       END IF

    ELSE
       !--------------------
       !B) COLUMN is > CMAX
       !--------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------
          !B.1) TEMPERATURE > TMAX, COLUMN > CMAX
          !---------------------------------------
          T1 = TKIN(M-1)
          T2 = TKIN(M)
          C1 = COLUMN(N-1)
          C2 = COLUMN(N)
          DFDT = (COOLFN(M,N) - COOLFN(M-1,N)) / (T2 - T1)
          DFDC = (COOLFN(M,N) - COOLFN(M,N-1)) / (C2 - C1)
          COOL = COOLFN(M,N) + TPLUS * DFDT + CPLUS * DFDC
          IER = 11
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------
             !B.2) TEMPERATURE < TMIN, COLUMN > CMAX
             !---------------------------------------
             T1 = TKIN(1)
             T2 = TKIN(2)
             C1 = COLUMN(N-1)
             C2 = COLUMN(N)
             DFDT = (COOLFN(2,N) - COOLFN(1,N)) / (T2 - T1)
             DFDC = (COOLFN(1,N) - COOLFN(1,N-1)) / (C2 - C1)
             COOL = COOLFN(1,N) - TMINUS * DFDT + CPLUS * DFDC
             IER = -11
          ELSE
             !---------------------------------------------
             !B.3) COLUMN > CMAX, TEMPERATURE WITHIN RANGE
             !---------------------------------------------
             C1 = COLUMN(N-1)
             C2 = COLUMN(N)
             COOL1 = SPLIN2(TKIN, COLUMN, COOLFN, Y2A, T, C1)
             COOL2 = SPLIN2(TKIN, COLUMN, COOLFN, Y2A, T, C2)
             DFDC = (COOL2 - COOL1) / (C2 - C1)
             COOL = COOL2 + CPLUS * DFDC
             IER = 1
          END IF
       END IF

    ENDIF

    COOL = 10._DP**(-COOL)

    RETURN
  END SUBROUTINE ROTATIONAL_CO



  SUBROUTINE VIBRATIONAL_CO(DENS,COL,TEMP,COOL,IER)
    !---------------------------------------------------------------------------
    ! called by :
    !     COMPUTE_MOLECULAR_COOLING
    ! purpose :
    !     Calculate the cooling rate (erg.cm3.s-1) due to vibrational CO. Adapted
    !     from Neufeld & Kaufman 1993.
    !     Warning : to obtain the cooling rate in erg/cm3/s, multiply by the
    !               number density (cm-3) of the two reactants (H2 and CO).
    ! subroutine/function needed :
    !     SPLIE2 (from Numerical Recipies in F90)
    !     SPLIN2 (from Numerical Recipies in F90)
    ! input variables :
    !     DENS -> H2 number density (cm-3)
    !     COL  -> optical depth parameter = G * n(CO)/Vgrad (cm-2.km-1.s)
    !             where G is the geometrical factor from Neufeld & Melnick 1991
    !     TEMP -> kinetic temperature of neutrals (K)
    ! output variables :
    !     COOL -> cooling rate (erg.cm3.s-1)
    !     IER  -> warning flag :
    !                0  => NORMAL EXECUTION
    !                1  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !               11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO LARGE, > 2000)
    !              -11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO SMALL, < 10)
    !               10  => TEMP OUT OF RANGE    (TOO LARGE, > 2000)
    !              -10  => TEMP OUT OF RANGE    (TOO SMALL, < 10)
    !
    !             WHEN IER .NE. 0 THE VALUE RETURNED FOR COOL IS
    !             AN EXTRAPOLATION OF DUBIOUS ACCURACY
    ! results :
    !---------------------------------------------------------------------------
    USE NUM_REC, ONLY          : SPLIE2, SPLIN2
    USE MODULE_CONSTANTS, ONLY : Zero

    IMPLICIT NONE

    REAL(KIND=DP),      INTENT(in)           :: DENS
    REAL(KIND=DP),      INTENT(in)           :: COL
    REAL(KIND=DP),      INTENT(in)           :: TEMP
    REAL(KIND=DP),      INTENT(out)          :: COOL
    INTEGER(KIND=LONG), INTENT(out)          :: IER

    INTEGER(KIND=LONG), PARAMETER            :: M = 6
    INTEGER(KIND=LONG), PARAMETER            :: N = 13
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: TKIN
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: CMAX
    REAL(KIND=DP),      SAVE, DIMENSION(N)   :: COLUMN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: F
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: COOLFN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: Y2A
    INTEGER(KIND=LONG)                       :: I
    INTEGER(KIND=LONG)                       :: J
    INTEGER(KIND=LONG), SAVE                 :: Ncall = 0
    REAL(KIND=DP)                            :: ARG
    REAL(KIND=DP)                            :: DCR
    REAL(KIND=DP)                            :: FACTOR
    REAL(KIND=DP)                            :: T
    REAL(KIND=DP)                            :: C
    REAL(KIND=DP)                            :: TMINUS
    REAL(KIND=DP)                            :: TPLUS
    REAL(KIND=DP)                            :: CPLUS
    REAL(KIND=DP)                            :: T1
    REAL(KIND=DP)                            :: T2
    REAL(KIND=DP)                            :: COOL1
    REAL(KIND=DP)                            :: COOL2
    REAL(KIND=DP)                            :: DFDT
    REAL(KIND=DP)                            :: C1
    REAL(KIND=DP)                            :: C2
    REAL(KIND=DP)                            :: DFDC

    Ncall = Ncall + 1_LONG
    ! initialize the arrays only on first call to subroutine
    IF (Ncall == 1_LONG) THEN
       TKIN = (/ 2.0000, 2.3010, 2.6021, 3.0000, 3.3010, 3.6021 /)
       COLUMN = &
            (/13.0000,13.5000,14.0000,14.5000,15.0000,15.5000 &
            , 16.0000,16.5000,17.0000,17.5000,18.0000,18.5000 &
            , 19.0000 /)
       F = RESHAPE( &
            (/10.8332,10.8250,10.8207,10.8004,10.7676,10.7946 &
            , 10.8343,10.8257,10.8213,10.8007,10.7677,10.7946 &
            , 10.8377,10.8281,10.8230,10.8016,10.7681,10.7948 &
            , 10.8485,10.8357,10.8283,10.8044,10.7694,10.7954 &
            , 10.8817,10.8591,10.8448,10.8132,10.7734,10.7974 &
            , 10.9771,10.9288,10.8950,10.8403,10.7858,10.8035 &
            , 11.2040,11.1101,11.0336,10.9192,10.8237,10.8227 &
            , 11.5669,11.4463,11.3281,11.1102,10.9300,10.8807 &
            , 11.9894,11.8611,11.7309,11.4176,11.1701,11.0403 &
            , 12.4374,12.3067,12.1648,11.7715,11.5437,11.3661 &
            , 12.8940,12.7656,12.6068,12.1752,11.9528,11.7480 &
            , 13.3559,13.2305,13.0346,12.6002,12.3092,12.0552 &
            , 13.8216,13.6996,13.4403,13.0312,12.6270,12.3148 /) &
            ,(/M,N/))
       DO I=1,M
          ARG = -68.0_DP * 10._DP**(-TKIN(I)/3._DP)
          CMAX(I) = -LOG10(1.83E-26 * (10._DP**TKIN(I)) * EXP(ARG))
       END DO
    END IF

    DO I = 1, M
       DO J = 1, N
          DCR = 10._DP**(CMAX(I) - F(I,J))
          FACTOR = 1._DP + DENS / DCR
          COOLFN(I,J) = LOG10(FACTOR) + CMAX(I)
       END DO
    END DO

    CALL SPLIE2(TKIN,COLUMN,COOLFN,Y2A)
    !        [2-D SPLINE FIT FROM NUMERICAL RECIPES]

    T = LOG10(TEMP + 1._DP)
    C = LOG10(COL + 1._DP)
    IF (C <= COLUMN(1)) C = COLUMN(1)

    TMINUS = TKIN(1) - T
    TPLUS  = T - TKIN(M)
    CPLUS  = C - COLUMN(N)

    IF (CPLUS <= Zero) THEN
       !--------------------------
       !A) COLUMN is within range
       !--------------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------------
          !A.1) TEMPERATURE > TMAX, COLUMN WITHIN RANGE
          !---------------------------------------------
          T1=TKIN(M-1)
          T2=TKIN(M)
          COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
          COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
          DFDT = (COOL2 - COOL1) / (T2 - T1)
          COOL = COOL2 + TPLUS * DFDT
          IER = 10
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------------
             !A.2) TEMPERATURE < TMIN, COLUMN WITHIN RANGE
             !---------------------------------------------
             T1 = TKIN(1)
             T2 = TKIN(2)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
             DFDT = (COOL2 - COOL1) / (T2 - T1)
             COOL = COOL1 - TMINUS * DFDT
             IER = -10
          ELSE
             !-----------------------------------------
             !A.3) TEMPERATURE AND COLUMN WITHIN RANGE
             !-----------------------------------------
             COOL = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C)
             IER = 0
          END IF
       END IF

    ELSE
       !--------------------
       !B) COLUMN is > CMAX
       !--------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------
          !B.1) TEMPERATURE > TMAX, COLUMN > CMAX
          !---------------------------------------
          T1 = TKIN(M-1)
          T2 = TKIN(M)
          C1 = COLUMN(N-1)
          C2 = COLUMN(N)
          DFDT = (COOLFN(M,N) - COOLFN(M-1,N)) / (T2 - T1)
          DFDC = (COOLFN(M,N) - COOLFN(M,N-1)) / (C2 - C1)
          COOL = COOLFN(M,N) + TPLUS * DFDT + CPLUS * DFDC
          IER = 11
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------
             !B.2) TEMPERATURE < TMIN, COLUMN > CMAX
             !---------------------------------------
             T1 = TKIN(1)
             T2 = TKIN(2)
             C1 = COLUMN(N-1)
             C2 = COLUMN(N)
             DFDT = (COOLFN(2,N) - COOLFN(1,N)) / (T2 - T1)
             DFDC = (COOLFN(1,N) - COOLFN(1,N-1)) / (C2 - C1)
             COOL = COOLFN(1,N) - TMINUS * DFDT + CPLUS * DFDC
             IER = -11
          ELSE
             !---------------------------------------------
             !B.3) COLUMN > CMAX, TEMPERATURE WITHIN RANGE
             !---------------------------------------------
             C1 = COLUMN(N-1)
             C2 = COLUMN(N)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C1)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C2)
             DFDC = (COOL2 - COOL1) / (C2 - C1)
             COOL = COOL2 + CPLUS * DFDC
             IER = 1
          END IF
       END IF

    ENDIF

    COOL = 10._DP**(-COOL) * EXP(-3080. / TEMP)

    RETURN
  END SUBROUTINE VIBRATIONAL_CO



  SUBROUTINE ROTATIONAL_o_H2O(DENS,COL,TEMP,COOL,IER)
    !---------------------------------------------------------------------------
    ! called by :
    !     COMPUTE_MOLECULAR_COOLING
    ! purpose :
    !     Calculate the cooling rate (erg.cm3.s-1) due to rotational ortho H2O.
    !     Adapted from Neufeld & Kaufman 1993.
    !     Warning : to obtain the cooling rate in erg/cm3/s, multiply by the
    !               number density (cm-3) of the two reactants (H2, ortho-H2O).
    ! subroutine/function needed :
    !     SPLIE2 (from Numerical Recipies in F90)
    !     SPLIN2 (from Numerical Recipies in F90)
    ! input variables :
    !     DENS -> H2 number density (cm-3)
    !     COL  -> optical depth parameter = G * n(ortho-H2O)/Vgrad (cm-2.km-1.s)
    !             where G is the geometrical factor from Neufeld & Melnick 1991
    !     TEMP -> kinetic temperature of neutrals (K)
    ! output variables :
    !     COOL -> cooling rate (erg.cm3.s-1)
    !     IER  -> warning flag :
    !                0  => NORMAL EXECUTION
    !                1  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !               11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO LARGE, > 2000)
    !              -11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO SMALL, < 10)
    !               10  => TEMP OUT OF RANGE    (TOO LARGE, > 2000)
    !              -10  => TEMP OUT OF RANGE    (TOO SMALL, < 10)
    !
    !             WHEN IER .NE. 0 THE VALUE RETURNED FOR COOL IS
    !             AN EXTRAPOLATION OF DUBIOUS ACCURACY
    ! results :
    !---------------------------------------------------------------------------
    USE NUM_REC, ONLY          : SPLIE2, SPLIN2
    USE MODULE_CONSTANTS, ONLY : Zero

    IMPLICIT NONE

    REAL(KIND=DP),      INTENT(in)           :: DENS
    REAL(KIND=DP),      INTENT(in)           :: COL
    REAL(KIND=DP),      INTENT(in)           :: TEMP
    REAL(KIND=DP),      INTENT(out)          :: COOL
    INTEGER(KIND=LONG), INTENT(out)          :: IER

    INTEGER(KIND=LONG), PARAMETER            :: M = 11
    INTEGER(KIND=LONG), PARAMETER            :: N = 10
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: TKIN
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: CMAX
    REAL(KIND=DP),      SAVE, DIMENSION(N)   :: COLUMN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: F
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: DENSHF
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: COOLFN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: ABEST
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: Y2A
    INTEGER(KIND=LONG)                       :: I
    INTEGER(KIND=LONG)                       :: J
    INTEGER(KIND=LONG), SAVE                 :: Ncall = 0
    REAL(KIND=DP)                            :: DH
    REAL(KIND=DP)                            :: DCR
    REAL(KIND=DP)                            :: FACTOR
    REAL(KIND=DP)                            :: T
    REAL(KIND=DP)                            :: C
    REAL(KIND=DP)                            :: TMINUS
    REAL(KIND=DP)                            :: TPLUS
    REAL(KIND=DP)                            :: CPLUS
    REAL(KIND=DP)                            :: T1
    REAL(KIND=DP)                            :: T2
    REAL(KIND=DP)                            :: COOL1
    REAL(KIND=DP)                            :: COOL2
    REAL(KIND=DP)                            :: DFDT
    REAL(KIND=DP)                            :: C1
    REAL(KIND=DP)                            :: C2
    REAL(KIND=DP)                            :: DFDC

    Ncall = Ncall + 1_LONG
    ! initialize the arrays only on first call to subroutine
    IF (Ncall == 1_LONG) THEN
       TKIN = &
            (/ 1.0000, 1.3010, 1.4771, 1.6990, 1.9031, 2.0000 &
            ,  2.3010, 2.6021, 3.0000, 3.3010, 3.6021 /)
       COLUMN = &
            (/ 10.0000,11.0000,12.0000,13.0000,14.0000,15.0000 &
            , 16.0000,17.0000,18.0000,19.0000 /)
       CMAX = &
            (/26.8133,25.8779,25.4274,24.9556,24.5820,24.4095 &
            , 23.9075,23.4457,22.8830,22.4955,22.1327 /)
       F = RESHAPE( &
            (/17.9430,16.7079,16.0839,15.4130,14.8546,14.6010 &
            , 13.8550,13.1585,12.3099,11.8499,11.6265,17.9647 &
            , 16.7240,16.0947,15.4183,14.8569,14.6026,13.8555 &
            , 13.1587,12.3099,11.8499,11.6265,18.1372,16.8552 &
            , 16.1854,15.4660,14.8784,14.6174,13.8604,13.1604 &
            , 12.3103,11.8501,11.6267,18.7671,17.3607,16.5781 &
            , 15.7181,15.0250,14.7287,13.9041,13.1763,12.3145 &
            , 11.8522,11.6281,19.6995,18.1087,17.2472,16.2714 &
            , 15.4680,15.1142,14.1295,13.2889,12.3515,11.8723 &
            , 11.6415,20.6659,18.9338,18.0505,17.0116,16.1183 &
            , 15.7238,14.6107,13.6399,12.5498,12.0183,11.7515 &
            , 21.5792,19.8109,18.9124,17.8186,16.8770,16.4531 &
            , 15.2414,14.1789,13.0024,12.4860,12.1796,22.5254 &
            , 20.7085,19.8006,18.6861,17.6987,17.2497,15.9570 &
            , 14.8195,13.6346,13.1447,12.8265,23.5005,21.6410 &
            , 20.7233,19.5788,18.5638,18.0994,16.7349,15.5296 &
            , 14.3564,13.8718,13.5436,24.4090,22.5793,21.6545 &
            , 20.4941,19.4566,18.9882,17.5593,16.2934,15.1454 &
            , 14.6656,14.3246 /), (/M,N/))
       DENSHF = RESHAPE( &
            (/ 8.8146, 8.9046, 9.0329, 9.1395, 9.1941, 9.2045 &
            ,  9.2142, 9.2875, 9.5279, 9.6410, 9.5598, 8.7860 &
            ,  8.8834, 9.0142, 9.1261, 9.2018, 9.1987, 9.2115 &
            ,  9.2861, 9.5275, 9.6407, 9.5597, 8.5971, 8.7280 &
            ,  8.8790, 9.0306, 9.1277, 9.1546, 9.1940, 9.2728 &
            ,  9.5233, 9.6388, 9.5583, 7.9562, 8.1364, 8.3130 &
            ,  8.5537, 8.7818, 8.8484, 9.0048, 9.1722, 9.4851 &
            ,  9.6198, 9.5450, 7.0116, 7.2052, 7.3961, 7.6862 &
            ,  7.9962, 8.1062, 8.3807, 8.6956, 9.2521, 9.4831 &
            ,  9.4483, 6.0247, 6.2048, 6.4102, 6.7077, 7.0394 &
            ,  7.1694, 7.4782, 7.8482, 8.5825, 8.9688, 9.0183 &
            ,  5.0340, 5.2088, 5.4174, 5.7028, 6.0567, 6.1900 &
            ,  6.4966, 6.8840, 7.6846, 8.1364, 8.2468, 4.0185 &
            ,  4.2137, 4.4058, 4.7068, 5.0519, 5.1894, 5.5405 &
            ,  5.9310, 6.7523, 7.2340, 7.3716, 3.0239, 3.2201 &
            ,  3.4100, 3.7112, 4.0584, 4.1951, 4.5393, 4.9793 &
            ,  5.8272, 6.3430, 6.4944, 2.0317, 2.2080, 2.4158 &
            ,  2.7040, 3.0633, 3.2001, 3.5415, 3.9836, 4.9060 &
            ,  5.4459, 5.6036 /), (/M,N/))
       ABEST = RESHAPE( &
            (/ 0.7050, 0.4860, 0.4760, 0.4550, 0.4520, 0.4540 &
            ,  0.4410, 0.4050, 0.3780, 0.3500, 0.3400, 0.6430 &
            ,  0.4870, 0.4750, 0.4540, 0.4500, 0.4530, 0.4400 &
            ,  0.4050, 0.3780, 0.3500, 0.3400, 0.5660, 0.4970 &
            ,  0.4740, 0.4460, 0.4390, 0.4410, 0.4330, 0.4000 &
            ,  0.3750, 0.3480, 0.3390, 0.5590, 0.5330, 0.4780 &
            ,  0.4410, 0.4180, 0.4140, 0.4040, 0.3770, 0.3590 &
            ,  0.3380, 0.3330, 0.5890, 0.5970, 0.5270, 0.4670 &
            ,  0.4470, 0.4290, 0.3870, 0.3480, 0.3340, 0.3240 &
            ,  0.3110, 0.7160, 0.6430, 0.5750, 0.5170, 0.4880 &
            ,  0.4700, 0.4010, 0.3530, 0.3240, 0.3030, 0.2950 &
            ,  0.8450, 0.6790, 0.6110, 0.5580, 0.5230, 0.5000 &
            ,  0.4220, 0.3570, 0.3250, 0.2890, 0.2910, 0.8570 &
            ,  0.7000, 0.6320, 0.5800, 0.5500, 0.5280, 0.4430 &
            ,  0.3710, 0.3170, 0.2760, 0.2770, 0.9270, 0.7170 &
            ,  0.6550, 0.6050, 0.5720, 0.5480, 0.4580, 0.3860 &
            ,  0.3100, 0.2640, 0.2670, 0.8700, 0.7270, 0.6720 &
            ,  0.6220, 0.5920, 0.5640, 0.4750, 0.3980, 0.3090 &
            ,  0.2570, 0.2580 /), (/M,N/))
    END IF

    DO I=1,M
       DO J=1,N
          DH=10._DP**DENSHF(I,J)
          DCR=10._DP**(CMAX(I)-F(I,J))
          FACTOR=1._DP+(1._DP-DH/DCR)*(DENS/DH)**ABEST(I,J)+DENS/DCR
          COOLFN(I,J)=LOG10(FACTOR)+CMAX(I)
       END DO
    END DO

    CALL SPLIE2(TKIN,COLUMN,COOLFN,Y2A)
    !        [2-D SPLINE FIT FROM NUMERICAL RECIPES]

    T=LOG10(TEMP+1._DP)
    C=LOG10(COL+1._DP)
    IF (C <= COLUMN(1)) C=COLUMN(1)

    TMINUS=TKIN(1)-T
    TPLUS=T-TKIN(M)
    CPLUS=C-COLUMN(N)

    IF (CPLUS <= Zero) THEN
       !--------------------------
       !A) COLUMN is within range
       !--------------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------------
          !A.1) TEMPERATURE > TMAX, COLUMN WITHIN RANGE
          !---------------------------------------------
          T1=TKIN(M-1)
          T2=TKIN(M)
          COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
          COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
          DFDT=(COOL2-COOL1)/(T2-T1)
          COOL=COOL2+TPLUS*DFDT
          IER=10
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------------
             !A.2) TEMPERATURE < TMIN, COLUMN WITHIN RANGE
             !---------------------------------------------
             T1=TKIN(1)
             T2=TKIN(2)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
             DFDT=(COOL2-COOL1)/(T2-T1)
             COOL=COOL1-TMINUS*DFDT
             IER=-10
          ELSE
             !-----------------------------------------
             !A.3) TEMPERATURE AND COLUMN WITHIN RANGE
             !-----------------------------------------
             COOL = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C)
             IER=0
          END IF
       END IF

    ELSE
       !--------------------
       !B) COLUMN is > CMAX
       !--------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------
          !B.1) TEMPERATURE > TMAX, COLUMN > CMAX
          !---------------------------------------
          T1=TKIN(M-1)
          T2=TKIN(M)
          C1=COLUMN(N-1)
          C2=COLUMN(N)
          DFDT=(COOLFN(M,N)-COOLFN(M-1,N))/(T2-T1)
          DFDC=(COOLFN(M,N)-COOLFN(M,N-1))/(C2-C1)
          COOL=COOLFN(M,N)+TPLUS*DFDT+CPLUS*DFDC
          IER=11
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------
             !B.2) TEMPERATURE < TMIN, COLUMN > CMAX
             !---------------------------------------
             T1=TKIN(1)
             T2=TKIN(2)
             C1=COLUMN(N-1)
             C2=COLUMN(N)
             DFDT=(COOLFN(2,N)-COOLFN(1,N))/(T2-T1)
             DFDC=(COOLFN(1,N)-COOLFN(1,N-1))/(C2-C1)
             COOL=COOLFN(1,N)-TMINUS*DFDT+CPLUS*DFDC
             IER=-11
          ELSE
             !---------------------------------------------
             !B.3) COLUMN > CMAX, TEMPERATURE WITHIN RANGE
             !---------------------------------------------
             C1=COLUMN(N-1)
             C2=COLUMN(N)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C1)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C2)
             DFDC=(COOL2-COOL1)/(C2-C1)
             COOL=COOL2+CPLUS*DFDC
             IER=1
          END IF
       END IF

    ENDIF

    COOL=10._DP**(-COOL)

    RETURN
  END SUBROUTINE ROTATIONAL_o_H2O



  SUBROUTINE ROTATIONAL_p_H2O(DENS,COL,TEMP,COOL,IER)
    !---------------------------------------------------------------------------
    ! called by :
    !     COMPUTE_MOLECULAR_COOLING
    ! purpose :
    !     Calculate the cooling rate (erg.cm3.s-1) due to rotational para H2O.
    !     Adapted from Neufeld & Kaufman 1993.
    !     Warning : to obtain the cooling rate in erg/cm3/s, multiply by the
    !               number density (cm-3) of the two reactants (H2, para-H2O).
    ! subroutine/function needed :
    !     SPLIE2 (from Numerical Recipies in F90)
    !     SPLIN2 (from Numerical Recipies in F90)
    ! input variables :
    !     DENS -> H2 number density (cm-3)
    !     COL  -> optical depth parameter = G * n(para-H2O)/Vgrad (cm-2.km-1.s)
    !             where G is the geometrical factor from Neufeld & Melnick 1991
    !     TEMP -> kinetic temperature of neutrals (K)
    ! output variables :
    !     COOL -> cooling rate (erg.cm3.s-1)
    !     IER  -> warning flag :
    !                0  => NORMAL EXECUTION
    !                1  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !               11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO LARGE, > 2000)
    !              -11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO SMALL, < 10)
    !               10  => TEMP OUT OF RANGE    (TOO LARGE, > 2000)
    !              -10  => TEMP OUT OF RANGE    (TOO SMALL, < 10)
    !
    !             WHEN IER .NE. 0 THE VALUE RETURNED FOR COOL IS
    !             AN EXTRAPOLATION OF DUBIOUS ACCURACY
    ! results :
    !---------------------------------------------------------------------------
    USE NUM_REC, ONLY          : SPLIE2, SPLIN2
    USE MODULE_CONSTANTS, ONLY : Zero

    IMPLICIT NONE

    REAL(KIND=DP),      INTENT(in)           :: DENS
    REAL(KIND=DP),      INTENT(in)           :: COL
    REAL(KIND=DP),      INTENT(in)           :: TEMP
    REAL(KIND=DP),      INTENT(out)          :: COOL
    INTEGER(KIND=LONG), INTENT(out)          :: IER

    INTEGER(KIND=LONG), PARAMETER            :: M = 11
    INTEGER(KIND=LONG), PARAMETER            :: N = 10
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: TKIN
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: CMAX
    REAL(KIND=DP),      SAVE, DIMENSION(N)   :: COLUMN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: F
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: DENSHF
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: COOLFN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: ABEST
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: Y2A
    INTEGER(KIND=LONG)                       :: I
    INTEGER(KIND=LONG)                       :: J
    INTEGER(KIND=LONG), SAVE                 :: Ncall = 0
    REAL(KIND=DP)                            :: DH
    REAL(KIND=DP)                            :: DCR
    REAL(KIND=DP)                            :: FACTOR
    REAL(KIND=DP)                            :: T
    REAL(KIND=DP)                            :: C
    REAL(KIND=DP)                            :: TMINUS
    REAL(KIND=DP)                            :: TPLUS
    REAL(KIND=DP)                            :: CPLUS
    REAL(KIND=DP)                            :: T1
    REAL(KIND=DP)                            :: T2
    REAL(KIND=DP)                            :: COOL1
    REAL(KIND=DP)                            :: COOL2
    REAL(KIND=DP)                            :: DFDT
    REAL(KIND=DP)                            :: C1
    REAL(KIND=DP)                            :: C2
    REAL(KIND=DP)                            :: DFDC

    Ncall = Ncall + 1_LONG
    ! initialize the arrays only on first call to subroutine
    IF (Ncall == 1_LONG) THEN
       TKIN = &
            (/ 1.0000, 1.3010, 1.4771, 1.6990, 1.9031, 2.0000 &
            ,  2.3010, 2.6021, 3.0000, 3.3010, 3.6021 /)
       COLUMN = &
            (/ 10.0000,11.0000,12.0000,13.0000,14.0000,15.0000 &
            ,  16.0000,17.0000,18.0000,19.0000 /)
       CMAX = &
            (/ 27.0084,25.7266,25.2370,24.7535,24.3835,24.2204 &
            ,  23.7683,23.3615,22.8559,22.5105,22.1685 /)
       F = RESHAPE( &
            (/17.7188,16.5987,16.1178,15.4290,14.8554,14.6011 &
            , 13.8546,13.1583,12.3143,11.8709,11.6589,17.7649 &
            , 16.6308,16.1317,15.4337,14.8576,14.6026,13.8551 &
            , 13.1584,12.3143,11.8709,11.6589,18.0667,16.8504 &
            , 16.2372,15.4766,14.8791,14.6176,13.8600,13.1601 &
            , 12.3148,11.8712,11.6590,18.8277,17.4079,16.6075 &
            , 15.7195,15.0271,14.7299,13.9038,13.1761,12.3190 &
            , 11.8734,11.6605,19.6829,18.1066,17.2624,16.2778 &
            , 15.4655,15.1121,14.1293,13.2888,12.3564,11.8944 &
            , 11.6748,20.5034,18.9375,18.0628,17.0066,16.1166 &
            , 15.7227,14.6106,13.6399,12.5576,12.0473,11.7905 &
            , 21.3702,19.8284,18.9297,17.8209,16.8716,16.4486 &
            , 15.2402,14.1788,13.0206,12.5312,12.2342,22.2794 &
            , 20.7513,19.8098,18.6870,17.6972,17.2468,15.9559 &
            , 14.8195,13.6690,13.2048,12.8955,23.2186,21.6919 &
            , 20.7282,19.5788,18.5595,18.0941,16.7348,15.5299 &
            , 14.4035,13.9381,13.6169,24.1895,22.5979,21.6502 &
            , 20.4940,19.4548,18.9730,17.5575,16.2959,15.2021 &
            , 14.7392,14.4037/), (/M,N/))
       DENSHF = RESHAPE( &
            (/ 9.3047, 9.1061, 8.9438, 8.7043, 8.5504, 8.5013 &
            ,  8.5510, 8.9173, 9.4095, 9.7067, 9.6889, 9.2489 &
            ,  9.0644, 8.9094, 8.6866, 8.5401, 8.4932, 8.5474 &
            ,  8.9153, 9.4090, 9.7065, 9.6888, 8.9411, 8.8157 &
            ,  8.7081, 8.5623, 8.4589, 8.4262, 8.5111, 8.8938 &
            ,  9.4043, 9.7045, 9.6875, 8.1687, 8.1647, 8.1607 &
            ,  8.1240, 8.1016, 8.1051, 8.2903, 8.7469, 9.3605 &
            ,  9.6861, 9.6758, 7.2136, 7.2444, 7.2844, 7.2992 &
            ,  7.3193, 7.3469, 7.6304, 8.2053, 9.1021, 9.5561 &
            ,  9.5805, 6.2043, 6.2434, 6.3141, 6.3313, 6.3485 &
            ,  6.3817, 6.6865, 7.3206, 8.4106, 9.0582, 9.1958 &
            ,  5.2092, 5.2498, 5.3043, 5.3238, 5.3537, 5.3878 &
            ,  5.6949, 6.3408, 7.4979, 8.2310, 8.4485, 4.2163 &
            ,  4.2571, 4.3087, 4.3281, 4.3630, 4.4025, 4.7395 &
            ,  5.3800, 6.5922, 7.3399, 7.6012, 3.2262, 3.2438 &
            ,  3.3149, 3.3340, 3.3661, 3.4124, 3.7878, 4.4911 &
            ,  5.6413, 6.4750, 6.7446, 2.2079, 2.2486, 2.3032 &
            ,  2.3249, 2.3670, 2.4128, 2.7899, 3.5141, 4.7516 &
            ,  5.5530, 5.8545 /), (/M,N/))
       ABEST =  RESHAPE( &
            (/ 0.4900, 0.7190, 0.6890, 0.5280, 0.4620, 0.4410 &
            ,  0.4000, 0.3730, 0.3520, 0.3570, 0.3710, 0.5200 &
            ,  0.6460, 0.6640, 0.5250, 0.4600, 0.4400, 0.3990 &
            ,  0.3720, 0.3510, 0.3560, 0.3710, 0.4900, 0.6540 &
            ,  0.6280, 0.5090, 0.4500, 0.4300, 0.3920, 0.3670 &
            ,  0.3480, 0.3550, 0.3700, 0.5740, 0.7250, 0.6460 &
            ,  0.4910, 0.4340, 0.4120, 0.3700, 0.3460, 0.3390 &
            ,  0.3460, 0.3630, 0.7520, 0.6760, 0.6440, 0.5180 &
            ,  0.4450, 0.4180, 0.3550, 0.3230, 0.3200, 0.3250 &
            ,  0.3390, 0.7560, 0.6950, 0.6760, 0.5500, 0.4690 &
            ,  0.4320, 0.3530, 0.3190, 0.3120, 0.3050, 0.3220 &
            ,  0.7880, 0.7330, 0.6940, 0.5750, 0.4780, 0.4380 &
            ,  0.3620, 0.3240, 0.3100, 0.2870, 0.3130, 0.8030 &
            ,  0.7530, 0.7110, 0.5830, 0.4950, 0.4610, 0.3840 &
            ,  0.3370, 0.3050, 0.2740, 0.2980, 0.8190, 0.7680 &
            ,  0.7270, 0.6020, 0.5170, 0.4820, 0.4050, 0.3590 &
            ,  0.2940, 0.2680, 0.2900, 0.9060, 0.7970, 0.7360 &
            ,  0.6170, 0.5340, 0.5010, 0.4230, 0.3730, 0.2970 &
            ,  0.2610, 0.2830 /), (/M,N/))
    END IF

    DO I=1,M
       DO J=1,N
          DH=10._DP**DENSHF(I,J)
          DCR=10._DP**(CMAX(I)-F(I,J))
          FACTOR=1._DP+(1._DP-DH/DCR)*(DENS/DH)**ABEST(I,J)+DENS/DCR
          COOLFN(I,J)=LOG10(FACTOR)+CMAX(I)
       END DO
    END DO

    CALL SPLIE2(TKIN,COLUMN,COOLFN,Y2A)
    !        [2-D SPLINE FIT FROM NUMERICAL RECIPES]


    T=LOG10(TEMP+1._DP)
    C=LOG10(COL+1._DP)
    IF (C <= COLUMN(1)) C=COLUMN(1)

    TMINUS=TKIN(1)-T
    TPLUS=T-TKIN(M)
    CPLUS=C-COLUMN(N)

    IF (CPLUS <= Zero) THEN
       !--------------------------
       !A) COLUMN is within range
       !--------------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------------
          !A.1) TEMPERATURE > TMAX, COLUMN WITHIN RANGE
          !---------------------------------------------
          T1=TKIN(M-1)
          T2=TKIN(M)
          COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
          COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
          DFDT=(COOL2-COOL1)/(T2-T1)
          COOL=COOL2+TPLUS*DFDT
          IER=10
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------------
             !A.2) TEMPERATURE < TMIN, COLUMN WITHIN RANGE
             !---------------------------------------------
             T1=TKIN(1)
             T2=TKIN(2)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
             DFDT=(COOL2-COOL1)/(T2-T1)
             COOL=COOL1-TMINUS*DFDT
             IER=-10
          ELSE
             !-----------------------------------------
             !A.3) TEMPERATURE AND COLUMN WITHIN RANGE
             !-----------------------------------------
             COOL = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C)
             IER=0
          END IF
       END IF

    ELSE
       !--------------------
       !B) COLUMN is > CMAX
       !--------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------
          !B.1) TEMPERATURE > TMAX, COLUMN > CMAX
          !---------------------------------------
          T1=TKIN(M-1)
          T2=TKIN(M)
          C1=COLUMN(N-1)
          C2=COLUMN(N)
          DFDT=(COOLFN(M,N)-COOLFN(M-1,N))/(T2-T1)
          DFDC=(COOLFN(M,N)-COOLFN(M,N-1))/(C2-C1)
          COOL=COOLFN(M,N)+TPLUS*DFDT+CPLUS*DFDC
          IER=11
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------
             !B.2) TEMPERATURE < TMIN, COLUMN > CMAX
             !---------------------------------------
             T1=TKIN(1)
             T2=TKIN(2)
             C1=COLUMN(N-1)
             C2=COLUMN(N)
             DFDT=(COOLFN(2,N)-COOLFN(1,N))/(T2-T1)
             DFDC=(COOLFN(1,N)-COOLFN(1,N-1))/(C2-C1)
             COOL=COOLFN(1,N)-TMINUS*DFDT+CPLUS*DFDC
             IER=-11
          ELSE
             !---------------------------------------------
             !B.3) COLUMN > CMAX, TEMPERATURE WITHIN RANGE
             !---------------------------------------------
             C1=COLUMN(N-1)
             C2=COLUMN(N)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C1)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C2)
             DFDC=(COOL2-COOL1)/(C2-C1)
             COOL=COOL2+CPLUS*DFDC
             IER=1
          END IF
       END IF

    ENDIF

    COOL=10._DP**(-COOL)

    RETURN
  END SUBROUTINE ROTATIONAL_p_H2O



  SUBROUTINE VIBRATIONAL_H2O(DENS,COL,TEMP,COOL,IER)
    !---------------------------------------------------------------------------
    ! called by :
    !     COMPUTE_MOLECULAR_COOLING
    ! purpose :
    !     Calculate the cooling rate (erg.cm3.s-1) due to vibrational H2O.
    !     Adapted from Neufeld & Kaufman 1993.
    !     Warning : to obtain the cooling rate in erg/cm3/s, multiply by the
    !               number density (cm-3) of the two reactants (H2, H2O).
    ! subroutine/function needed :
    !     SPLIE2 (from Numerical Recipies in F90)
    !     SPLIN2 (from Numerical Recipies in F90)
    ! input variables :
    !     DENS -> H2 number density (cm-3)
    !     COL  -> optical depth parameter = G * n(H2O)/Vgrad (cm-2.km-1.s)
    !             where G is the geometrical factor from Neufeld & Melnick 1991
    !     TEMP -> kinetic temperature of neutrals (K)
    ! output variables :
    !     COOL -> cooling rate (erg.cm3.s-1)
    !     IER  -> warning flag :
    !                0  => NORMAL EXECUTION
    !                1  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !               11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO LARGE, > 2000)
    !              -11  => COL OUT OF RANGE     (TOO LARGE, > 1E19)
    !                      & TEMP OUT OF RANGE  (TOO SMALL, < 10)
    !               10  => TEMP OUT OF RANGE    (TOO LARGE, > 2000)
    !              -10  => TEMP OUT OF RANGE    (TOO SMALL, < 10)
    !
    !             WHEN IER .NE. 0 THE VALUE RETURNED FOR COOL IS
    !             AN EXTRAPOLATION OF DUBIOUS ACCURACY
    ! results :
    !---------------------------------------------------------------------------
    USE NUM_REC, ONLY          : SPLIE2, SPLIN2
    USE MODULE_CONSTANTS, ONLY : Zero

    IMPLICIT NONE

    REAL(KIND=DP),      INTENT(in)           :: DENS
    REAL(KIND=DP),      INTENT(in)           :: COL
    REAL(KIND=DP),      INTENT(in)           :: TEMP
    REAL(KIND=DP),      INTENT(out)          :: COOL
    INTEGER(KIND=LONG), INTENT(out)          :: IER

    INTEGER(KIND=LONG), PARAMETER            :: M = 6
    INTEGER(KIND=LONG), PARAMETER            :: N = 13
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: TKIN
    REAL(KIND=DP),      SAVE, DIMENSION(M)   :: CMAX
    REAL(KIND=DP),      SAVE, DIMENSION(N)   :: COLUMN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: F
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: COOLFN
    REAL(KIND=DP),      SAVE, DIMENSION(M,N) :: Y2A
    INTEGER(KIND=LONG)                       :: I
    INTEGER(KIND=LONG)                       :: J
    INTEGER(KIND=LONG), SAVE                 :: Ncall = 0
    REAL(KIND=DP)                            :: ARG
    REAL(KIND=DP)                            :: DCR
    REAL(KIND=DP)                            :: FACTOR
    REAL(KIND=DP)                            :: T
    REAL(KIND=DP)                            :: C
    REAL(KIND=DP)                            :: TMINUS
    REAL(KIND=DP)                            :: TPLUS
    REAL(KIND=DP)                            :: CPLUS
    REAL(KIND=DP)                            :: T1
    REAL(KIND=DP)                            :: T2
    REAL(KIND=DP)                            :: COOL1
    REAL(KIND=DP)                            :: COOL2
    REAL(KIND=DP)                            :: DFDT
    REAL(KIND=DP)                            :: C1
    REAL(KIND=DP)                            :: C2
    REAL(KIND=DP)                            :: DFDC

    Ncall = Ncall + 1_LONG
    ! initialize the arrays only on first call to subroutine
    IF (Ncall == 1_LONG) THEN
       TKIN = (/ 2.0000, 2.3010, 2.6021, 3.0000, 3.3010, 3.6021 /)
       COLUMN = &
            (/13.0000,13.5000,14.0000,14.5000,15.0000,15.5000 &
            , 16.0000,16.5000,17.0000,17.5000,18.0000,18.5000 &
            , 19.0000 /)
       F = RESHAPE( &
            (/10.9764,11.0466,11.0684,10.8760,10.6678,10.6729 &
            , 10.9788,11.0478,11.0689,10.8761,10.6679,10.6729 &
            , 10.9862,11.0512,11.0706,10.8766,10.6681,10.6730 &
            , 11.0079,11.0619,11.0758,10.8781,10.6687,10.6732 &
            , 11.0633,11.0922,11.0914,10.8828,10.6706,10.6739 &
            , 11.1769,11.1666,11.1340,10.8965,10.6761,10.6760 &
            , 11.3605,11.3111,11.2310,10.9324,10.6900,10.6814 &
            , 11.6107,11.5330,11.4034,11.0113,10.7228,10.6939 &
            , 11.9106,11.8195,11.6492,11.1510,10.7934,10.7227 &
            , 12.2401,12.1485,11.9488,11.3545,10.9214,10.7855 &
            , 12.5862,12.4982,12.2799,11.6126,11.1120,10.9008 &
            , 12.9403,12.8519,12.6292,11.9132,11.3569,11.0728 &
            , 13.2898,13.2006,12.9932,12.2484,11.6474,11.2954/) &
            ,(/M,N/))
       DO I=1,M
          ARG=-47.5_DP*10._DP**(-TKIN(I)/3._DP)
          CMAX(I)=-LOG10(1.03E-26*(10._DP**TKIN(I))*EXP(ARG))
       END DO
    END IF

    DO I=1,M
       DO J=1,N
          DCR=10._DP**(CMAX(I)-F(I,J))
          FACTOR=1._DP+DENS/DCR
          COOLFN(I,J)=LOG10(FACTOR)+CMAX(I)
       END DO
    END DO

    CALL SPLIE2(TKIN,COLUMN,COOLFN,Y2A)
    !        [2-D SPLINE FIT FROM NUMERICAL RECIPES]


    T=LOG10(TEMP+1._DP)
    C=LOG10(COL+1._DP)
    IF (C <= COLUMN(1)) C=COLUMN(1)

    TMINUS=TKIN(1)-T
    TPLUS=T-TKIN(M)
    CPLUS=C-COLUMN(N)

    IF (CPLUS <= Zero) THEN
       !--------------------------
       !A) COLUMN is within range
       !--------------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------------
          !A.1) TEMPERATURE > TMAX, COLUMN WITHIN RANGE
          !---------------------------------------------
          T1=TKIN(M-1)
          T2=TKIN(M)
          COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
          COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
          DFDT=(COOL2-COOL1)/(T2-T1)
          COOL=COOL2+TPLUS*DFDT
          IER=10
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------------
             !A.2) TEMPERATURE < TMIN, COLUMN WITHIN RANGE
             !---------------------------------------------
             T1=TKIN(1)
             T2=TKIN(2)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T1,C)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T2,C)
             DFDT=(COOL2-COOL1)/(T2-T1)
             COOL=COOL1-TMINUS*DFDT
             IER=-10
          ELSE
             !-----------------------------------------
             !A.3) TEMPERATURE AND COLUMN WITHIN RANGE
             !-----------------------------------------
             COOL = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C)
             IER=0
          END IF
       END IF

    ELSE
       !--------------------
       !B) COLUMN is > CMAX
       !--------------------
       IF (TPLUS > Zero) THEN
          !---------------------------------------
          !B.1) TEMPERATURE > TMAX, COLUMN > CMAX
          !---------------------------------------
          T1=TKIN(M-1)
          T2=TKIN(M)
          C1=COLUMN(N-1)
          C2=COLUMN(N)
          DFDT=(COOLFN(M,N)-COOLFN(M-1,N))/(T2-T1)
          DFDC=(COOLFN(M,N)-COOLFN(M,N-1))/(C2-C1)
          COOL=COOLFN(M,N)+TPLUS*DFDT+CPLUS*DFDC
          IER=11
       ELSE
          IF (TMINUS > Zero) THEN
             !---------------------------------------
             !B.2) TEMPERATURE < TMIN, COLUMN > CMAX
             !---------------------------------------
             T1=TKIN(1)
             T2=TKIN(2)
             C1=COLUMN(N-1)
             C2=COLUMN(N)
             DFDT=(COOLFN(2,N)-COOLFN(1,N))/(T2-T1)
             DFDC=(COOLFN(1,N)-COOLFN(1,N-1))/(C2-C1)
             COOL=COOLFN(1,N)-TMINUS*DFDT+CPLUS*DFDC
             IER=-11
          ELSE
             !---------------------------------------------
             !B.3) COLUMN > CMAX, TEMPERATURE WITHIN RANGE
             !---------------------------------------------
             C1=COLUMN(N-1)
             C2=COLUMN(N)
             COOL1 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C1)
             COOL2 = SPLIN2(TKIN,COLUMN,COOLFN,Y2A,T,C2)
             DFDC=(COOL2-COOL1)/(C2-C1)
             COOL=COOL2+CPLUS*DFDC
             IER=1
          END IF
       END IF

    ENDIF

    COOL=10._DP**(-COOL)*EXP(-2325._DP/TEMP)

    RETURN
  END SUBROUTINE VIBRATIONAL_H2O



  SUBROUTINE W_ERR_COOL(error,procedure)
    !---------------------------------------------------------------------------
    ! called by :
    !     COMPUTE_MOLECULAR_COOLING
    ! purpose :
    !     Write an error message (if occured in calculation of the cooling due
    !     to CO or H2O) in the error file.
    !     Because the same error occurs several times, only the first one is
    !     written.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR

    IMPLICIT NONE

    INTEGER(KIND=LONG), INTENT(in) :: error
    CHARACTER(LEN=*),   INTENT(in) :: procedure
    INTEGER(KIND=LONG), SAVE       :: last_error = 0
    INTEGER(KIND=LONG), SAVE       :: last_last_error = 0

    IF (error /= last_error .AND. error /= last_last_error) THEN
       err_cool = .TRUE.
       WRITE(f_err_cool,'("step= ",I5," Tn= ",ES8.2," *** IER =",I3," in ",A)',&
            ADVANCE='NO') counter, Tn, error, procedure
       SELECT CASE (error)
       CASE (1)
          WRITE(f_err_cool,'("-> COL too large")')
       CASE (11)
          WRITE(f_err_cool,'("-> COL too large and Tn too large")')
       CASE (-11)
          WRITE(f_err_cool,'("-> COL too large and Tn too small")')
       CASE (10)
          WRITE(f_err_cool,'("-> Tn too large")')
       CASE (-10)
          WRITE(f_err_cool,'("-> Tn too small")')
       END SELECT
       last_last_error = last_error
       last_error = error
    END IF

  END SUBROUTINE W_ERR_COOL



  REAL(dp) FUNCTION analytical_cooling(n,T)
    USE MODULE_CONSTANTS, ONLY : kB

    IMPLICIT NONE

    REAL(dp), PARAMETER :: B  = 1.0e-11_dp
    REAL(dp), PARAMETER :: T1 = 1.0e+04_dp
    REAL(dp), PARAMETER :: T2 = 1.0e+01_dp
    REAL(dp), PARAMETER :: nc = 1.0e+00_dp
    REAL(dp)            :: n,T

    analytical_cooling = B * n * kB * ( T - T1/(1d0+(n/nc)**2) - T2 )

  END FUNCTION analytical_cooling



END MODULE MODULE_MOLECULAR_COOLING