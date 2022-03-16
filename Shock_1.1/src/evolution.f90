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

MODULE MODULE_EVOLUTION

  !*****************************************************************************
  !** The module 'MODULE_EVOLUTION' contains 3 subroutines (CHEMISTRY, SOURCE **
  !** and DIFFUN) used to calculate at each step :                            **
  !**           - the source terms (CHEMISTRY and SOURCE)                     **
  !**           - the partial derivatives for DVODE integration package.      **
  !*****************************************************************************
  ! changes made by SC since last release 16II04: search for "SC" in file
  ! 
  ! coding of the evolution.f90 file name:
  ! - kn: bug corrected in Newton solving for dVn/dz when KN=1 (2 May 2006)
  !       + minimum gradient = thermal gradient = Vsound/(Z+dmin)
  ! - rg: use only core species to calculate mean radius of grain core
  !       assuming a Matthis distribution (2 May 2006)
  ! - lay: add mantle layers to grain cross section assuming Nsites = cste
  !        using: rsquare_grain = (r_core + Nlayer * layer_thickness)^2
  ! modified variable: Nsites_grain = 5e4 instead of 1e6 (variables_mhd)
  ! new variables needed: layer_thickness = 2e-8 cm (variables_mhd)
  !                       R_gr_scale12 (variables_mhd, calculated in initialize)
  ! - G: species C60,C60+,C60- replaced by G,G+,G-
  ! Tabone 09-18: heat_exchange terms are now computed in CHEMISTRY to include any
  ! type of reaction
  USE tex
  USE MODULE_TOOLS, ONLY : MINIMUM, MAXIMUM
  IMPLICIT NONE
  INCLUDE "precision.f90"
  INTEGER, PARAMETER :: QP = selected_real_kind(P=33) ! quadruple precision

  
  !-----------------------------------------------
  ! vectors containing variables and derivatives
  ! at the last call to DIFFUN
  ! dimension = 1:d_v_var
  ! allocated and initialized in MHD
  !-----------------------------------------------
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: v_dvariab ! Y * dY at the last call of DIFFUN
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: v_l_var   ! Y at the last call of DIFFUN
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: v_l_der   ! DY at the last call to DIFFUN

  REAL(KIND=DP), PARAMETER :: tiny = 1.0D-40 ! avoids overflow in DIFFUN

  !------------------------------------------------------------------
  ! source terms for each fluid
  ! (N -> neutrals, I -> ions, NEG -> negative ions)
  ! = rates of change in number, mass, momentum and energy densities
  ! calculated in CHEMISTRY, and in SOURCE
  !------------------------------------------------------------------
  REAL(KIND=DP) :: CupN, CupI, CupA, CupNEG  ! increase in number density (cm-3.s-1)
  REAL(KIND=DP) :: CdownN, CdownI, CdownA, CdownNEG  ! decrease in number density (cm-3.s-1)
  REAL(KIND=DP) :: YNn, YNi, YNa, YNneg      ! change in number density (cm-3.s-1)
  REAL(KIND=DP) :: Sn, Si, Sa, Sneg          ! change in mass density (g.cm-3.s-1)
  REAL(KIND=DP) :: An, Ai, Aa, Aneg          ! change in momentum density (g.cm-2.s-2) !neutral, ionized and electron fluids
  REAL(KIND=DP) :: A_CH_n, YA_CH, YS_CH      ! note: A_CH_n=change in momentum density (g.cm-2.s-2) for CH
  REAL(KIND=DP) :: A_S_n,  YA_S,  YS_S       ! note: A_S_n =change in momentum density (g.cm-2.s-2) for S
  REAL(KIND=DP) :: A_SH_n, YA_SH, YS_SH      ! note: A_SH_n=change in momentum density (g.cm-2.s-2) for SH
  REAL(KIND=DP) :: Bn, Bi, Ba, Bneg          ! change in energy density (erg.cm-3.s-1)
  REAL(KIND=DP) :: B_heat_exchange_n , B_heat_exchange_i, B_heat_exchange_neg ! change in enthalpy btw fluids (Tabone 09-18)
  REAL(KIND=DP) :: H2_int_energy             ! internal energy of H2 (K cm-3)
  !------------------------------
  ! source terms for each species
  ! dimension = 0:Nspec
  ! calculated in CHEMISTRY
  !------------------------------
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: YN  ! change in number density (cm-3.s-1)
  REAL(KIND=QP), DIMENSION(:), SAVE, ALLOCATABLE :: YNQ ! change in number density (cm-3.s-1) - quadruple preicsion

  PRIVATE :: tiny!, SOURCE ! source is needed bby thermochemistry

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE CHEMISTRY
    !---------------------------------------------------------------------------
    ! called by :
    !     DIFFUN
    ! purpose :
    !     calculates source terms related to chemistry. Each type of reaction
    !     has its own reaction rate (with different meanings and units for
    !     gamma, alpha and beta ).
    !     For each species :
    !         * YN : rate at which it is formed through chemical reactions
    !                [1], eqs. (14) and (20)
    !         * YS : rate of change in mass (g.cm-3.s-1)
    !         * YA : rate of change in momentum (g.cm-2.s-2)
    !         * YB : rate of change in energy (erg.cm-3.s-1)
    !     For each fluid (x= n -> neutrals, i -> ions, neg -> negative ions and e-)
    !         * CupX   : rate of increase in number density (cm-3.s-1)
    !         * CdownX : rate of decrease in number density (cm-3.s-1)
    !         * Nx     : rate of change in number density (cm-3.s-1)
    !         * Sx     : rate of change in mass density (g.cm-3.s-1)
    !         * Ax     : rate of change in momentum density (g.cm-2.s-2)
    !         * Bx     : rate of change in energy density (erg.cm-3.s-1)
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     YN
    !     CupN,   CupI,   CupNEG
    !     CdownN, CdownI, CdownNEG
    !     YNn,    YNi,    YNa,    YNneg
    !     Sn,     Si,     Sa,     Sneg
    !     An,     Ai,     Aa,     Aneg
    !     Bn,     Bi,     Ba,     Bneg
    ! references :
    !     [1] Flower et al., 1985, MNRAS 216, 775
    !     [2] Flower et al., 1986, MNRAS 218, 729
    !     [3] Pineau des Forets et al., 1986, MNRAS 220, 801
    !     [4] Viala et al., 1988, A&A 190, 215
    !---------------------------------------------------------------------------
    USE MODULE_CHEM_REACT
    USE MODULE_GRAINS
    USE MODULE_H2
    USE MODULE_CO
    ! USE MODULE_PHYS_VAR,    ONLY : Av, RAD, Zeta, Te, Tn, Vi, Vn, nH, ABS_DeltaV, NH2_lev
    USE MODULE_PHYS_VAR
    USE MODULE_CONSTANTS,   ONLY : pi, kB, me, EVerg, Zero, qe2
    USE MODULE_DEBUG_JLB
    USE MODULE_SHIELD
    USE MODULE_RADIATION,   ONLY : fluphfree, fluph, fluph0
    USE MODULE_GAMMA,       ONLY : GAMMA1

    IMPLICIT NONE

    ! index 0 -> inexistent species, Nspec_plus -> added species (e-, GRAIN, ...)
    REAL(KIND=DP), DIMENSION(0:Nspec_plus) :: up_N, up_MV, up_MV2, up_DE
    REAL(KIND=DP), DIMENSION(0:Nspec_plus) :: down_N, down_MV, down_MV2
    REAL(KIND=DP), DIMENSION(0:Nspec)      :: YS, YA, YB
    INTEGER(KIND=LONG) :: i, IR1, IR2, IR3, IP1, IP2, IP3, IP4                               ! modif Tabone 12/2017 3Body
    INTEGER(KIND=LONG) :: idx_raw
    REAL(KIND=DP)      :: massR1, massR2, massR3, densityR1, densityR2, densityR3            ! modif Tabone 12/2017 3Body
    REAL(KIND=DP)      :: V_R1, V_R2, V_R3, Mass_prod, Nprod_m1                              ! modif Tabone 12/2017 3Body
    REAL(KIND=DP)      :: Vcm, Vcm2
    REAL(KIND=DP)      :: creation_rate, destr_R1, destr_R2, destr_R3                        !  modif Tabone 12/2017 3Body
    REAL(KIND=DP)      :: creation_Vcm, creation_Vcm2, creation_DE
    REAL(KIND=DP)      :: destr1_Vcm, destr2_Vcm, destr3_Vcm, destr1_Vcm2, destr2_Vcm2, destr3_Vcm2  !  modif Tabone 12/2017 3Body
    REAL(KIND=DP)      :: exp_factor, exp_factor1, exp_factor2, coeff, Tr, Ts, energy, Teff
    REAL(KIND=DP)      :: react_rate, r_alph, r_beta, r_gamm
    REAL(KIND=DP)      :: fracH2
    REAL(KIND=DP)      :: shielding
    REAL(KIND=DP)      :: nui
    INTEGER            :: lev
    INTEGER            :: j
    REAL(KIND=DP)      :: dum
    REAL(KIND=DP)      :: sum_enth_reac ! used to compute B_heat_exchange_*

    REAL(KIND=DP)      :: stick_H_gr = 1.0_DP

    INTEGER(KIND=LONG) :: CR1, CR2
    REAL(KIND=DP)      :: adust
    REAL(KIND=DP)      :: rap_z
    REAL(KIND=DP)      :: tauph, cterec
    REAL(KIND=DP)      :: stick_e
    REAL(KIND=DP)      :: qi_ds, tau_ds, nu_ds
    REAL(KIND=DP)      :: thet_ds, Jtld_ds

    ! Test Tabone
    INTEGER            :: indexR1exo, indexR1endo, indexR2exo, indexR2endo
    REAL(KIND=DP)      :: creation_max, creation_min

    INTEGER            :: idxf, idxd

    creation_max = 0._DP
    creation_min = 0._DP

    ! initialization
    up_N     = Zero
    up_MV    = Zero
    up_MV2   = Zero
    up_DE    = Zero
    down_N   = Zero
    down_MV  = Zero
    down_MV2 = Zero

    Sel_ch_H2(1:NH2_lev) = Zero
    Sel_ne_H2(1:NH2_lev) = Zero
    Sel_io_H2(1:NH2_lev) = Zero
    Sel_el_H2(1:NH2_lev) = Zero

    B_heat_exchange_n   = Zero
    B_heat_exchange_i   = Zero
    B_heat_exchange_neg = Zero

    !=====================================================================================
    ! 1 - calculates creation & destruction terms for all reaction (i=idx of the reaction)
    !=====================================================================================
    DO i=1,Nreact
         !  write(*,*) '****CHECK REACTION number****', i, 'over', Nspec
         !  write(*,*) '****CHECK 3BODY REACT IR3****', react(i)%R(3)
       !--- index, mass, abundance, velocity of each reactant and product ---
       IR1       = react(i)%R(1)       ; IR2       = react(i)%R(2)
       IR3       = react(i)%R(3)                                                 !Tabone 12/2017 3Body
       IP1       = react(i)%P(1)       ; IP2       = react(i)%P(2)
       IP3       = react(i)%P(3)       ; IP4       = react(i)%P(4)
       massR1    = speci(IR1)%mass     ; massR2    = speci(IR2)%mass
       massR3    = speci(IR3)%mass
       densityR1 = speci(IR1)%density  ; densityR2 = speci(IR2)%density
       densityR3 = speci(IR3)%density
       V_R1      = speci(IR1)%velocity ; V_R2      = speci(IR2)%velocity
       V_R3      = speci(IR3)%velocity                                          ! Tabone 12/2917 3Body => velocity optional 
       Mass_prod = react(i)%Mass_prod  ! sum of the masses of the products
       Nprod_m1  = react(i)%Nprod_m1   ! number of products minus one

       r_alph = react(i)%alpha
       r_beta = react(i)%beta
       r_gamm = react(i)%gamma

       !--- center of mass velocity (cm/s) and its square ---
       Vcm  = ( massR1 * V_R1 + massR2 * V_R2 ) / ( massR1 + massR2 )
       Vcm2 = Vcm * Vcm

       !--- index of raw integration ---
       idx_raw = react(i)%idx_raw_integ

       ! ------------------------------------------------------
       ! now calculate the reaction rate (cm3.s-1) according to 
       ! the type of the reaction (effective temperature, 
       ! significance of gamma alpha, beta, ...).
       ! ------------------------------------------------------
       SELECT CASE (react(i)%type)

       ! ------------------------------------------------------
       ! 3Body GAS PHASE reactions
       ! ------------------------------------------------------
       CASE ('3BODY')
          !--- Here temperatures = termerature R1 => to be updated if ion / neutral /neutral 3Body reaction Tabone 12/2017 3Body
          Tr             = speci(IR1)%temperature
          !--- reaction rate (cm-3.s-1) ---
          exp_factor    = r_beta/Tr
          react_rate    = r_gamm * TEXP(-exp_factor) * ((Tr/300._DP)**r_alph)
          creation_rate = react_rate * densityR1 * densityR2 * densityR3

       ! ------------------------------------------------------
       ! photo reactions
       ! ------------------------------------------------------
       CASE ('PHOTO')
          !--- reaction rate coefficient (s-1) ---
          ! BG 2017-08
          !   analytical expressions for all reactions.
          !   fits of photoreaction rates as functions of Av have been obtained
          !   assuming a given extinction curve and grain properties which lead
          !   to a NH/Av ratio of 1.87e21. Changing grain properties greatly
          !   change the significance of Av compared to the photon flux between
          !   911 and 2000 A.
          !   Exponential fit of the photon flux predicted by the PDR code for 
          !   0 < Av < 3 with the PDR code gives
          !   ===> fluph = fluph0 * exp(-2.8 Av)
          !   for a galactic extinction curve. This means that the rates usually 
          !   written
          !   ===> k = G0 * gamma * exp(-beta*Av)
          !   can be modified to
          !   ===> k = G0 * gamma * (fluph / fluph0)**(beta / 2.8)
          !   We adopt these rates here because they depend on a more general
          !   and physical definition of the absorption of FUV photons

          IF      (F_COUP_RAD == 1) THEN
             react_rate = RAD * r_gamm * (fluph/fluph0)**(r_beta/2.8)
          ELSE IF (F_COUP_RAD == 2) THEN
             react_rate = r_gamm * fluph / fluphfree
          ELSE
             react_rate = RAD * r_gamm * TEXP(-r_beta*Av)
          ENDIF

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! photo reactions integrated over cross sections
       ! ------------------------------------------------------
       CASE ('CROSS')
          !--- reaction rate coefficient (s-1) ---
          ! BG 2017-08
          !     integration of photodestruction cross sections for a few reactions
          !     listed in chemistry.in

          IF(F_COUP_RAD == 1 .OR. F_COUP_RAD == 2) THEN
             react_rate = phdest_rate(idx_raw)
          ELSE
             react_rate = RAD * r_gamm * TEXP(-r_beta*Av)
          ENDIF

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! self-shielded reactions (H2 & CO photodiss)
       ! ------------------------------------------------------
       CASE ('SHILD')
          ! Identify molecule and set the dissociation rate
          IF (IR1==ind_H2) THEN
             IF(F_COUP_RAD == 1 .OR. F_COUP_RAD == 2) THEN
                react_rate = probdissH2
             ELSE
                IF (react(i)%ref == 'PDR16') THEN
                   ! How can we include the width of the shock (and the variation 
                   ! in density) in the computation of the photodissociation rate 
                   ! of H2 and CO, with the sole results of the Meudon PDR code ? 
                   ! For now, use data at the surface layer, as in steady-state
                   react_rate = photoH2
                ELSE
                   shielding= shield_H2(Av,N_H2_0+coldens_h2,react(i)%ref)
                   react_rate = r_gamm * RAD * shielding
                ENDIF
             ENDIF
          ELSE IF (IR1==ind_CO) THEN
             !!! -----------------------------------------------
             !!! There is a problem in using the FGK approx for 
             !!! computing the photodissociation rate of CO. FGK
             !!! works only if the doppler broadening is larger 
             !!! than natural broadening. As a result the self-
             !!! shielding of CO obtained with FGK approx is 
             !!! overestimated. It is far better for now to 
             !!! take the values derived by Lee et al. (1996)
             !!! => We could redo the work of FGK and extend it 
             !!!    to a broader range of scenarios in order to
             !!!    remove the comments belows - TO DO LATER
             !!! -----------------------------------------------
             !!! IF(F_COUP_RAD == 1 .OR. F_COUP_RAD == 2) THEN
             !!!    react_rate = probdissCO
             !!! ELSE
             !!! -----------------------------------------------
                IF (react(i)%ref == 'PDR16') THEN
                   ! How can we include the width of the shock (and the variation 
                   ! in density) in the computation of the photodissociation rate 
                   ! of H2 and CO, with the sole results of the Meudon PDR code ? 
                   ! For now, use data at the surface layer, as in steady-state
                   react_rate = photoCO
                ELSE
                   shielding= shield_CO(Av,N_H2_0+coldens_h2, N_CO_0+coldens_co,react(i)%ref)
                   ! shielding= shield_CO(Av,N_H2_0+coldens_h2,N_CO_0,react(i)%ref)
                   ! shielding= shield_CO(Av,N_H2_0,N_CO_0,react(i)%ref)
                   react_rate = r_gamm * RAD * shielding
                ENDIF
             !!! -----------------------------------------------
             !!! ENDIF
             !!! -----------------------------------------------
          ENDIF

          ! Save the photo-dissociation rates (s-1)
          IF (IR1==ind_H2) THEN
             photodiss_H2 = react_rate * densityR2
          ELSE
             photodiss_CO = react_rate * densityR2
          ENDIF

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! photoelectric reaction
       ! ------------------------------------------------------
       CASE ('PHELE')
          !--- reaction rate coefficient (s-1) ---
          ! BG 2017-08
          ! a - use the expression of Bakes and Tielens 1994 stored in phele_reac
          ! b - analytical expressions
          !     fits of photoreaction rates as functions of Av have been obtained
          !     assuming a given extinction curve and grain properties which lead
          !     to a NH/Av ratio of 1.87e21. Changing grain properties greatly
          !     change the significance of Av compared to the photon flux between
          !     911 and 2000 A.
          !     Exponential fit of the photon flux predicted by the PDR code for 
          !     0 < Av < 3 with the PDR code gives
          !     ===> fluph = fluph0 * exp(-2.8 Av)
          !     for a galactic extinction curve. This means that the rates usually 
          !     written
          !     ===> k = G0 * gamma * exp(-beta*Av)
          !     can be modified to
          !     ===> k = G0 * gamma * (fluph / fluph0)**(beta / 2.8)
          !     We adopt these rates here because they depend on a more general
          !     and physical definition of the absorption of FUV photons

          IF(F_COUP_RAD == 1) THEN
             IF ( idx_raw /= 0 ) THEN
                ! The rate computed in COMPUTE_PHOTOELEC is the exact calculation 
                ! of the photoelectric ejection rate
                ! However, the value can be too large to be correctly treated
                ! by dvode -> if so reduce the rate to a 0.01 micron size grain
                react_rate = phele_reac(idx_raw)%rate
                ! react_rate = phele_reac(idx_raw)%rate * 1e-12 / rsq_grm
             ELSE
                ! Careful, the r_gamm in the chemical file obtained for a grain 
                ! of 0.01 microns => equivalent rsq_grm = 1e-12
                ! Should be scaled by (rsq_grm / 1e-12). However, the value can then
                ! be too large to be correctly treated by dvode -> if so forget the scaling
                react_rate = RAD * r_gamm * (fluph/fluph0)**(r_beta/2.8) * (rsq_grm / 1e-12)
                ! react_rate = RAD * r_gamm * (fluph/fluph0)**(r_beta/2.8) ! * (rsq_grm / 1e-12)
             ENDIF
          ELSE IF(F_COUP_RAD == 2) THEN
             IF ( idx_raw /= 0 ) THEN
                react_rate = phele_reac(idx_raw)%rate
                ! react_rate = phele_reac(idx_raw)%rate * 1e-12 / rsq_grm
             ELSE
                react_rate = r_gamm * fluph / fluphfree * (rsq_grm / 1e-12)
                ! react_rate = r_gamm * fluph / fluphfree ! * (rsq_grm / 1e-12)
             ENDIF
          ELSE
             ! Careful, the r_gamm in the chemical file obtained for a grain
             ! of 0.01 microns => equivalent rsq_grm = 1e-12
             ! Should be scaled by (rsq_grm / 1e-12). However, the value can then
             ! be too large to be correctly treated by dvode -> If so forget the scaling
             react_rate = r_gamm * TEXP(-r_beta*Av) * RAD * (rsq_grm / 1e-12)
             ! react_rate = r_gamm * TEXP(-r_beta*Av) * RAD ! * (rsq_grm / 1e-12)
          ENDIF

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! Recombination on grains
       ! ------------------------------------------------------
       CASE ('GRREC')
          !--- reaction rate coefficient (s-1) ---

          ! ---------------------------------------------------
          ! Weingartner & Draine (2001)
          ! ---------------------------------------------------

          ! ---------------------------------------------------
          ! Choose the size of grains - 6.4 angstroms for PAH
          ! ---------------------------------------------------
          IF ( (IR1 == ind_G).OR.(IR1 == ind_Gminus).OR.(IR1 == ind_Gplus) ) THEN
             adust = r_grm
          ELSE
             ! adust = 6.4e-8_dp
             adust = 5.0e-8_dp
          ENDIF

          ! ---------------------------------------------------
          ! effective temperature of the reaction
          ! see [1], eq. (43) and [3], eqs.(2.1)-(2.3)
          ! for reactions between members of the same fluid
          ! Ts = 0, Teff = Tr
          ! ---------------------------------------------------
          tauph  = adust * kB * speci(IR2)%temperature / qe2
          cterec = SQRT(8.0_dp * kB * speci(IR2)%temperature * pi) * adust**2.0

          ! ---------------------------------------------------
          ! compute electron attachement sticking coefficient
          ! Eqs 28, 29 & 30 of Weingartner & Draine (2001)
          ! ---------------------------------------------------
          stick_e = 0.5_dp * ( 1.0_dp - EXP( -adust*1.e8_dp/10.0_dp ) )

          CR1   = speci(IR1)%charge
          CR2   = speci(IR2)%charge
          rap_z = DBLE(CR1) / DBLE(CR2)
          ! ---------------------------------------------------
          ! relevant variables for gaz-grain collisions
          ! Eq 3.1 of Draine & Sutin (1987)
          ! ---------------------------------------------------
          qi_ds  = DBLE(CR2)
          tau_ds = tauph / qi_ds / qi_ds
          nu_ds  = DBLE(CR1) / qi_ds
          IF ( rap_z < 0 ) THEN
             ! ------------------------------------------------
             ! Eqs. 2.4a, 3.4, & 3.9 of Draine & Sutin (1987)
             ! ------------------------------------------------
             thet_ds = 0.0_dp
             Jtld_ds = ( 1.0_dp - nu_ds / tau_ds ) &
                     * ( 1.0_dp + SQRT(2.0_dp/(tau_ds - 2.0_dp * nu_ds)) )
          ELSE IF ( rap_z == 0 ) THEN
             ! ------------------------------------------------
             ! Eqs. 2.4a, 3.3, & 3.8 of Draine & Sutin (1987)
             ! ------------------------------------------------
             thet_ds = 0.0_dp
             Jtld_ds = 1.0_dp + SQRT(pi / (2.0_dp * tau_ds))
          ELSE
             ! ------------------------------------------------
             ! Eqs. 2.4a, 3.5, & 3.10 of Draine & Sutin (1987)
             ! ------------------------------------------------
             thet_ds = nu_ds / (1.0_dp + 1.0_dp/SQRT(nu_ds))
             Jtld_ds = EXP( -thet_ds / tau_ds ) &
                     * ( 1.0_dp + 1.0_dp / SQRT(4.0_dp * tau_ds + 3.0_dp * nu_ds) )**2
          ENDIF
          ! ---------------------------------------------------
          ! attach. electrons on grain
          ! ---------------------------------------------------
          IF (IR2 == ind_e) THEN
             react_rate = cterec * stick_e * Jtld_ds / SQRT(me)
          ! ---------------------------------------------------
          ! collisions grain - ions
          ! ---------------------------------------------------
          ELSE
             react_rate = cterec * Jtld_ds / SQRT(massR2)
          ENDIF

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! cosmic ray ionization or dissociation
       ! ------------------------------------------------------
       CASE ('CR_IO')
          !--- reaction rate coefficient (cm3.s-1) ---
          ! We take into account :
          !     0) direct cosmic ray (CR) induced ionization/dissociation
          !     1) the CR induced secondary photons
          !     2) the UV photons produced by collisional excitation of H2 by electrons
          !        following shock heating
          !     3) allow for the differential compression of the charged grains and the
          !        neutrals (Vi/Vn). This factor is applicable to cases (1) and
          !        (2) above, but not to case (0). However, the direct CR induced
          !        processes are completely negligible within the shock, which is
          !        where Vi/Vn differs from 1.
          ! note that react(i)%beta is preset at a very large value (in REACTION_TYPE)
          ! for direct cosmic ray ionization : react(i)%beta = 1.0D8

          ! JLB + GP : 13 juin 2002  - correction de la fraction de H2 dans le gaz
          ! JLB     exp_factor = MIN(r_beta/Te,180._DP)
          exp_factor = r_beta / Te
          ! print*,'Te',Te,Av
          coeff = SQRT(8.0_DP*kB*Te/(pi*me)) * 1.0D-16 * (4.29D-1 + 6.14D-6*Te)
          IF (r_beta < 0.99D8) THEN
            ! proportionality factor, fracH2: 0.46 is the probability of cosmic ray ionization of H,
            ! 1.0 is the probability of cosmic ray ionization of H2
            ! The secondary photons are produced by radiative cascade, following collisional 
            ! excitation of Rydberg states of H2 by secondary electrons, generated by cosmic rays
            fracH2 = (0.46*Dens_H + Dens_H2) / (Dens_H + Dens_H2) * Dens_H2 / (Dens_H + dens_H2)
          ELSE
            fracH2 = 1.0_DP
          ENDIF
          !
          ! DRF+GDF: Jan-2017
          ! For the treatment of the grains (G, G-, G+), see F & PdesF 2003, MNRAS, 343, 390, Sect. 2.2.2
          ! Rate per unit volume of detachment of electrons from grains by the H_2 fluorescence photons (SECPHO):
          ! 0.15 zeta n_H f y (cm-3 s-1)
          ! where
          ! -- 0.15 is the fraction of the secondary electrons that collisionally excite the Rydberg states of H_2, 
          !    generating the secondary photons in the subsequent radiative cascade;
          ! -- zeta (s-1) is the rate of cosmic ray ionization of H_2;
          ! -- n_H = (n_H / n_g) n_g 
          ! -- f depends on the supposed fractional PAH abundance: f = 0.87 for n_PAH / n_H = 1e-6;
          !    f = 0.99 for n_PAH / n_H = 1e-7
          ! -- y = 0.03 for ionization of neutral grains, y = 0.2 for electron photo-detachment from negatively 
          !    charged grains (both values are uncertain).
          ! In the above paper, 
          ! n_H / n_g = 1 / 7.1e-11 = 1.4e10
          ! and hence the rates per unit volume of ionization/detachment of electrons from grains 
          ! are, when n_PAH / n_H = 0.0,
          ! zeta * n_g * 6.3e7 cm-3 s-1 for neutral grains
          ! and
          ! zeta * n_g- * 4.2e8 cm-3 s-1 for negatively charged grains
          ! The values of r_gamm must be scaled to the current value of n_H / n_g and also adapted 
          ! to the value of n_PAH / n_H
          IF (IR2 == ind_SECPHO) THEN
             IF ((IR1 == ind_Gminus).or.(IR1 == ind_G)) THEN
                r_gamm = r_gamm * 7.1e-11_dp / (dens_grc / (dens_H + 2.*dens_H2))
             ! Following the analysis of Gredel et al. 1989, ApJ, 347,289, the probabilities, p_M,
             ! of photo-ionization and photo-dissociation of atoms and molecules are inversely
             ! proportional to sigma_g, the grain extinction cross section per hydrogen nucleus;
             ! sigma_g = <n_g sigma> / n_H cm^2
             ! Gredel et al. adopt a value of sigma_g = 2e-21 cm^2, Thus, r_gamm must be scaled by
             ! a factor 2.e-21 / (<n_g sigma> / n_H)
             ELSE
                !*************  Tabone 03/2018  ************
                !Modification: opactitŽ du gaz 0.01*2e-21
                !r_gamm = r_gamm * 2.e-21_dp / (dens_grc * pi * rsq_grm / (dens_H + 2.*dens_H2))
                r_gamm = r_gamm * 1._dp /(1.e-2_dp + (dens_grc * pi * rsq_grm / (dens_H + 2.*dens_H2)) / 2.e-21_dp )
             ENDIF
          ENDIF

          react_rate = Zeta * ((Tn/300._DP)**r_alph) * r_gamm &
                     + Dens_e * coeff * TEXP(-exp_factor) / 0.15_DP * r_gamm * Vi / Vn
          react_rate = react_rate * fracH2

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! cosmic ray induced desorption from grains
       ! ------------------------------------------------------
       CASE ('CR_DE')
         !--- reaction rate (s-1) ---

         react_rate = r_gamm * dens_grc / ab_ads * pi * rsq_grm

         !  values of the rates of desorption, relative to that for CO, are obtained from:
         !         exp[-(E_ads - 855)/T_g(max)], where E_ads is the adsorption energy (K)
         !  of the species (Aikawa et al. 1996, ApJ, 467, 684, Table 1; see also Bergin et al.
         !  2002, ApJ, 570, L101) and 855 K is the adsorption energy of CO on CO ice
         !  T_g(max) = 70 K (Hasegawa & Herbst 1983, MNRAS, 261, 83, equ. (15))

         ! Note : in the formula below, and adsorption energy of 855 is adopted for CO
         ! based on Hasegawa. However, we assumed an adsorption energy of 960 K in the
         ! chemical file based on Aikawa (1996)
         ! This lead to a factor of 5 difference for the desorption of CO.

         ! This formalism and those desorption rates are based on the work of Leger 
         ! et al. (1985). We assume that the CR ionization rate of 1e-17 s-1 used by 
         ! these authors was that of H2.
         
         react_rate = zeta / 1e-17_dp * react_rate * exp(-(r_beta-855._dp)/70._dp)

         creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! H2 and HD formation on grains
       ! ------------------------------------------------------
       CASE ('H2_FO')
          !--- reaction rate (cm3.s-1) ---

          ! react_rate = r_gamm * nH/Dens_H * ((Tn/300._DP)**r_alph)

          ! JLB - 21 juin 2001 - Correction du taux de formation de H2

          !--- effective temperature ---
          Teff_grain = massR1*(Vi-Vn)*(Vi-Vn)/(3.0_DP*kB) + speci(IR1)%temperature
          stick_H_gr = 1.0_DP / sqrt(1.0_DP + Teff_grain / 30.0_DP)

          !--- reaction rate (cm3.s-1) (Andersson & Wannier for H adsorption) ---
          !    or try Hollenbach & McKee (1979, ApJ Suppl, 41, 555)
          stick_H_gr = 1.0_DP / (1.D0 + 0.04_DP*(Teff_grain+Tgrain)**0.5_DP &
                                 + 2.0D-3*Teff_grain + 8.0D-6*Teff_grain*Teff_grain)
          react_rate = stick_H_gr * pi * rsq_grm * &
               SQRT(8.0_DP*kB*Teff_grain/massR1/pi)

          creation_rate = react_rate * dens_grc * densityR2

          ! Formula of Le Bourlot et al. (2012) - BG
          ! Be careful
          ! -> to be tested more thoroughly before adoption
          ! react_rate = r_gamm &
          !            * 1.0_DP / (1.0_DP + speci(IR1)%temperature / 464.0_DP)**(1.5_DP) &
          !            * (speci(IR1)%temperature / 100._DP)**r_alph &
          !            * (Dens_H + 2.0_DP*Dens_H2 + Dens_Hplus) / Dens_H

          ! creation_rate = react_rate * densityR1 * densityR2

          if (IR2 == ind_H) then
            For_gr_H2 = creation_rate
          endif

       ! ------------------------------------------------------
       ! three-body reactions on grains surface
       ! ------------------------------------------------------
       CASE ('THREE')
          !--- effective temperature ---

          Teff_grain = massR1*(Vi-Vn)*(Vi-Vn)/(3.0_DP*kB) + speci(IR1)%temperature

          !--- reaction rate (cm3.s-1) (Andersson & Wannier for H adsorption) ---

          react_rate = r_gamm * pi * rsq_grm * SQRT(8.0_DP*kB*Teff_grain/massR1/pi)
          react_rate = react_rate / (Teff_grain/r_beta+1.0_DP)**2.0_DP
          react_rate = react_rate * dens_grc / ab_ads

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! sputtering of grain mantles
       ! ------------------------------------------------------
       CASE ('SPUTT')
          ! Flower & Pineau des Forets 1994 (MNRAS 268 724)
          ! Barlow (1978) MNRAS 183 367

          !--- effective temperature of grain sputtering ---

          Ts = massR2*(V_R2-V_R1)*(V_R2-V_R1)/(3.0_DP*kB)
          Tr =  speci(IR2)%temperature
          Teff_grain = Ts + Tr

          !--- reaction rate (cm3.s-1) ---
          ! reaction  : X* + IR2 -> X + IR2 + GRAIN
          ! Nlayers = number of layers in the mantles
          ! IF Nlayers > 1 => the detachment probability is independent of Nlayers
          ! IF Nlayers < 1 => the detachment probability is proportional to Nlayers

          Nlayers = d_ads / d_site
          exp_factor1   = (4.0_DP*r_beta-1.5_DP*Ts)/Tr
          exp_factor2   = 4.0_DP*r_beta/Teff_grain
          exp_factor    = MAXIMUM(exp_factor1,exp_factor2)
          ! JLB
          ! exp_factor    = MIN(exp_factor,180._DP)

          react_rate = 4.0_DP * r_gamm * dens_grc/ab_ads &
                     * pi * rsq_grm * SQRT( 8.0_DP * kB * Teff_grain / massR2 / pi)
          react_rate = react_rate * (1.0_DP+2.0_DP*Teff_grain/(4.0_DP*r_beta)) * &
               TEXP(-exp_factor)
          ! smooth MIN function
          ! react_rate = react_rate * MIN(Nlayers,1.0_DP)
          react_rate = react_rate * 1.0_DP / sqrt( 1.0_DP / Nlayers**2.0_DP + 1.0_DP )

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! erosion of grain cores
       ! ------------------------------------------------------
       CASE ('EROSI')
          !--- reaction rate (cm3.s-1) ---
          react_rate = Zero
          !--- absolute value of the ion-neutral velocity drift ---

          energy = massR2 * ABS_DeltaV*ABS_DeltaV / 2.0_DP / EVerg ! energy (eV)

          IF (energy > r_alph) THEN
          react_rate = pi * rsq_grc * ABS_DeltaV * &
               r_gamm * TEXP(-r_beta/(Energy-r_alph))
          END IF

          ! the rate coefficient is proportional to the grain density
          ! and not to the core species density

          ! Original form:
          !         react_rate = react_rate * dens_grc / densityR1

          IF (densityR1 > dens_grc) THEN
            creation_rate = react_rate * dens_grc * densityR2
          ELSE
            creation_rate = react_rate * densityR1 * densityR2
          ENDIF

       ! ------------------------------------------------------
       ! adsorption on to grains
       ! ------------------------------------------------------
       CASE ('ADSOR')
          !--- effective temperature (K) ---

          Teff_grain    = massR1*(Vi-Vn)*(Vi-Vn)/(3.0_DP*kB) + speci(IR1)%temperature

          !--- reaction rate (cm3.s-1) ---
          ! here, gamma = sticking coefficient
          ! react_rate = STICK * SIGMA(GRAIN) * V(MOYENNE)

          react_rate = r_gamm * pi * dens_grc * rsq_grm * SQRT(8.0_DP*kB*Teff_grain/massR1/pi)
          IF(r_beta /=Zero) THEN
             ! Tielens & Hollenbach
             react_rate = react_rate / (1.0_DP + 0.04_DP*(Teff_grain+Tgrain)**0.5_DP + &
                  2.0D-3 * Teff_grain + 8.0D-6*Teff_grain*Teff_grain)
          END IF

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! photodesorption
       ! ------------------------------------------------------
       CASE ('PHDES')

          ! Hollenbach  et al., 2009, ApJ 690, 1497
          ! Equations (6) and (7)

          !--- reaction rate (cm3.s-1) ---
          ! reaction  : X* + PHOTONS -> X + GRAIN
          ! Nlayers = number of layers in the mantles
          ! IF Nlayers > 1 => the detachment probability is independent of Nlayers
          ! IF Nlayers < 1 => the detachment probability is proportional to Nlayers

          Nlayers = d_ads / d_site

          IF(F_COUP_RAD == 1 .OR. F_COUP_RAD == 2) THEN
             ! direct calculation of photon flux
             dum = fluph
          ELSE
             ! Equation (7) of Hollenbach et al. (2009)
             dum = RAD * fluphfree * exp(-1.8*Av)
          ENDIF

          react_rate = r_gamm * dum * pi * rsq_grm * dens_grc / ab_ads
          ! react_rate = r_gamm * dum * 4 * pi * rsq_grm * dens_grc / ab_ads
          ! smooth MIN function
          ! react_rate = react_rate * MIN(Nlayers,1.0_DP)
          react_rate = react_rate * 1.0_DP / sqrt( 1.0_DP / Nlayers**2.0_DP + 1.0_DP )

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! photodesorption
       ! ------------------------------------------------------
       CASE ('SECDE')

          !--- reaction rate (cm3.s-1) ---
          ! reaction  : X* + SECPHO -> X + GRAIN
          ! Nlayers = number of layers in the mantles
          ! IF Nlayers > 1 => the detachment probability is independent of Nlayers
          ! IF Nlayers < 1 => the detachment probability is proportional to Nlayers

          Nlayers = d_ads / d_site

          ! proportionality factor, fracH2: 0.46 is the probability of cosmic ray ionization of H,
          ! 1.0 is the probability of cosmic ray ionization of H2
          ! The secondary photons are produced by radiative cascade, following collisional 
          ! excitation of Rydberg states of H2 by secondary electrons, generated by cosmic rays
          fracH2 = (0.46*Dens_H + Dens_H2) / (Dens_H + Dens_H2) * Dens_H2 / (Dens_H + dens_H2)

          ! --- old treatment
          ! Hollenbach  et al., 2009, ApJ 690, 1497
          ! page 1503, last paragraph of first column
          ! IMPORTANT - to use, set r_gamma = 1e-3 in the chemical file
          ! react_rate = r_gamm * fphsec * frach2 * pi * rsq_grm * dens_grc / ab_ads
          ! --- new treatment
          ! According to DRF+GDF Jan-2017 (see above), the ionization by secpho of G- 
          ! and G occurs with a yield of 0.2 and 0.03 respectively. Taking into account
          ! the fraction of secondary electrons that collisionally excite the Rydberg
          ! states (0.15), this lead to a rate
          ! zeta * n_g * 6.3e7 cm-3 s-1 for neutral grains
          ! 6.3e7 is the r_gamma set in the input file (which include yield * 0.15)
          ! The desorption by secpho of atoms and molecules occurs with a yield of 1e-3
          ! The r_gamma set in the input file must therefore be
          ! 6.3e7 * 1e-3 / 0.03 = 2.1e6
          ! The desorption rate by secpho has therefore the same expression as ionization
          ! of grains
          r_gamm = r_gamm * 7.1e-11_dp / (dens_grc / (dens_H + 2.*dens_H2))
          react_rate = Zeta * ((Tn/300._DP)**r_alph) * r_gamm * frach2 * dens_grc / ab_ads

          ! smooth MIN function
          ! react_rate = react_rate * MIN(Nlayers,1.0_DP)
          react_rate = react_rate * 1.0_DP / sqrt( 1.0_DP / Nlayers**2.0_DP + 1.0_DP )

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! thermal desorption
       ! ------------------------------------------------------
       CASE ('THDES')

          ! Hollenbach  et al., 2009, ApJ 690, 1497
          ! Equations (2) & following paragraph
          ! Binding energy from Hasegawa & Herbst (1993)
          !                     Aikawa et al. (1996)
          !                 and Saptarsy
          ! -> see chemistry.in (r_beta)

          !--- reaction rate (cm3.s-1) ---
          ! reaction  : X* + VOISIN -> X + GRAIN
          ! the input grain temperature is assumed
          ! -> better solution would be to find the
          !    temperature that mimic the integration
          !    of the rate over the MRN distribution
          !    TO DO

          nui = 1.6e11_dp * sqrt( r_beta / (massR1 / mass_H) )
          exp_factor = r_beta / Tgrain

          react_rate = nui * TEXP(-exp_factor)

          creation_rate = react_rate * densityR1 * densityR2

       ! ------------------------------------------------------
       ! C+ + H2 for each excitation level of H2
       ! ------------------------------------------------------
       CASE ('EXCIT')
          !--- effective temperature of the reaction                               ---
          !--- see [1], eq. (43) and [3], eqs.(2.1)-(2.3)                          ---
          !--- for reactions between members of the same fluid : Ts = 0, Teff = Tr ---

          Ts = massR1*massR2*(V_R1-V_R2)*(V_R1-V_R2)/(massR1+massR2)/3.0_DP/kB
          Tr = (massR1*speci(IR2)%temperature + massR2*speci(IR1)%temperature) / &
               (massR1+massR2)
          Teff = Ts + Tr

          !--- destruction is computed state by state                             ---
          !--- reduce dissociation energy r_beta by the excitation energy         ---
          !--- of the level
          Sel_rx_H2=0_DP
          do lev=1,NH2_lev
             ! Line 2 and 3 of table 1 in Agundez et al (2010)
             ! case v=0
             if (H2_lev(lev)%v==0 .and.H2_lev(lev)%j<8) then
                ! Line 2 : v=0, j=0..7
                exp_factor1 = (r_beta - H2_lev(lev)%energy - 3.0_DP*Ts) / Tr
                exp_factor2 = (r_beta - H2_lev(lev)%energy) / Teff
                exp_factor  = MAXIMUM(exp_factor1,exp_factor2)
                exp_factor  = MAXIMUM(exp_factor,0.0_DP)
                Sel_rx_H2(lev) = TEXP(-exp_factor)
                !--- reaction rate (cm3.s-1) ---
                Sel_rx_H2(lev) = r_gamm * Sel_rx_H2(lev) * ((Teff/300._DP)**r_alph)
             else
                ! Line 3 : v=1
                !--- reaction rate (cm3.s-1) ---
                Sel_rx_H2(lev) = 1.6e-9_DP
             endif
             ! Other levels are simply discarded
          end do
 
          !--- reaction rate (cm3.s-1) ---
          react_rate = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_rx_H2))/Dens_H2

          creation_rate = react_rate * densityR2 *densityR1

       ! ------------------------------------------------------
       ! Collisional dissociation of H2
       ! ------------------------------------------------------
       CASE ('DISSO')
          !--- effective temperature of the reaction                               ---
          !--- see [1], eq. (43) and [3], eqs.(2.1)-(2.3)                          ---
          !--- for reactions between members of the same fluid : Ts = 0, Teff = Tr ---

          Ts = massR1*massR2*(V_R1-V_R2)*(V_R1-V_R2)/(massR1+massR2)/3.0_DP/kB
          Tr = (massR1*speci(IR2)%temperature + massR2*speci(IR1)%temperature) / &
               (massR1+massR2)
          Teff = Ts + Tr

          !--- destruction is computed state by state                             ---
          !--- reduce dissociation energy r_beta by the excitation energy         ---
          !--- of the level

          do lev=1,NH2_lev

            !--- exponential factor ([3], PAGE 810) ---
            exp_factor1 = (r_beta - H2_lev(lev)%energy - 3.0_DP*Ts) / Tr
            exp_factor2 = (r_beta - H2_lev(lev)%energy) / Teff
            exp_factor  = MAXIMUM(exp_factor1,exp_factor2)
            exp_factor  = MAXIMUM(exp_factor,0.0_DP)

            Sel_rx_H2(lev) = TEXP(-exp_factor)

          end do

          !--- reaction rate (cm3.s-1) ---
          Sel_rx_H2 = r_gamm * Sel_rx_H2 * ((Teff/300._DP)**r_alph)

          !--- Note: we avoid dividing by n(H2) just to be able to remultiply...

          react_rate = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_rx_H2))

          creation_rate = react_rate * densityR2

          Sel_ch_H2 = Sel_ch_H2 - Sel_rx_H2 * densityR2
          if (IR2 >= b_neu .AND. IR2 <= e_neu) then
            Sel_ne_H2 = Sel_ne_H2 - Sel_rx_H2 * densityR2
          else if (IR2 >= b_ion .AND. IR2 <= e_ion) then
            Sel_io_H2 = Sel_io_H2 - Sel_rx_H2 * densityR2
          else if (IR2 == ind_e) then
            Sel_el_H2 = Sel_el_H2 - Sel_rx_H2 * densityR2
          endif

       ! ------------------------------------------------------
       ! all other reactions (including reverse reactions)
       ! ------------------------------------------------------
       CASE ('OTHER', 'REVER')
          !--- effective temperature of the reaction                               ---
          !--- see [1], eq. (43) and [3], eqs.(2.1)-(2.3)                          ---
          !--- for reactions between members of the same fluid : Ts = 0, Teff = Tr ---

          Ts = massR1*massR2*(V_R1-V_R2)*(V_R1-V_R2)/(massR1+massR2)/3.0_DP/kB
          Tr = (massR1*speci(IR2)%temperature + massR2*speci(IR1)%temperature) / &
               (massR1+massR2)
          Teff = Ts + Tr

          !--- exponential factor ([3], PAGE 810) ---
          exp_factor1 = (r_beta-3.0_DP*Ts) / Tr
          exp_factor2 = r_beta / Teff
          exp_factor  = MAXIMUM(exp_factor1,exp_factor2)
          ! JLB 
          ! exp_factor  = MIN(exp_factor,180._DP)

          !--- reaction rate (cm3.s-1) ---
          react_rate = r_gamm * TEXP(-exp_factor) * ((Teff/300._DP)**r_alph)

          ! include temperature dependent Coulomb enhancement factor in the cases of neutralization
          ! of charged grains, G+ and G-, by electrons and positive ions, respectively
          if (IP1 == ind_G) react_rate = react_rate * (1._DP + 450._DP/Teff)

          creation_rate = react_rate * densityR1 * densityR2

       END SELECT ! end of selection for the reaction rate

       !--- Save the reaction rate in the react structure ---
       react(i)%rate = creation_rate

       !--- creation rate (cm-3.s-1), destruction rate (s-1) of each reactant ---
       !  creation computed in select case above
       !  multiplication by density done already
       !      creation_rate = react_rate * densityR1 * densityR2
       !      destr_R1      = react_rate * densityR2
       !      destr_R2      = react_rate * densityR1

       destr_R1      = creation_rate
       destr_R2      = creation_rate
       destr_R3      = creation_rate                                ! Tabone 12/2017 3Body
       
       !--- source terms for destruction (multiplied by mass later) ---
       destr1_Vcm   = destr_R1 * Vcm
       destr2_Vcm   = destr_R2 * Vcm
       destr3_Vcm   = destr_R3 * Vcm                                ! Tabone 12/2017 3Body
       
       destr1_Vcm2  = destr_R1 * Vcm2
       destr2_Vcm2  = destr_R2 * Vcm2
       destr3_Vcm2  = destr_R3 * Vcm2                               ! Tabone 12/2017 3Body
       
       !--- source terms for creation ---
       creation_Vcm  = creation_rate * Vcm
       creation_Vcm2 = creation_rate * Vcm2
       creation_DE   = creation_rate * react(i)%DE

       !--- for the reactants, increase in the rate of destruction (s-1) of : ---
       !---      down_N(IR)   -> number                                       ---
       !---      down_MV(IR)  -> momentum                                     ---
       !---      down_MV2(IR) -> energy                                       ---

       down_N(IR1)   = down_N(IR1)   + destr_R1
       down_N(IR2)   = down_N(IR2)   + destr_R2
       down_N(IR3)   = down_N(IR3)   + destr_R3                       ! Tabone 12/2017 3Body
       
       down_MV(IR1)  = down_MV(IR1)  + destr1_Vcm
       down_MV(IR2)  = down_MV(IR2)  + destr2_Vcm
       down_MV(IR3)  = down_MV(IR3)  + destr3_Vcm                     ! Tabone 12/2017 3Body      
       
       down_MV2(IR1) = down_MV2(IR1) + destr1_Vcm2
       down_MV2(IR2) = down_MV2(IR2) + destr2_Vcm2
       down_MV2(IR3) = down_MV2(IR3) + destr3_Vcm2                    ! Tabone 12/2017 3Body       

       !--- for the products, increase in the rate of creation (s-1) of :    ---
       !---      up_N(IP)   -> number                                        ---
       !---      up_MV(IP)  -> momentum                                      ---
       !---      up_MV2(IP) -> kinetic energy                                ---
       !---      up_DE(IP)  -> energy defect of the reaction                 ---
       !---                    DE is distributed over the reaction products  ---
       !---                    with a weighing factor : [1], eq. (31)        ---
       !---                        =  DE/(n-1)*[1-mass(i)/MASS] if n > 1     ---
       !---                        =  DE if n=1                              ---
       !---                    where n=number of produts, i=index of product ---
       !---                    and MASS=sum of mass(i), i=1..n               ---
       !---                    By this way, energy is conserved whatever     ---
       !---                    the number of products in the reaction.       ---

       ! first product, always present !
       up_N(IP1)   = up_N(IP1)   + creation_rate
       up_MV(IP1)  = up_MV(IP1)  + creation_Vcm
       up_MV2(IP1) = up_MV2(IP1) + creation_Vcm2


       IF (Nprod_m1 == 0) THEN
          ! allow for reactions with only one product : H2_FO, ADSOR

          ! JLB - May 2001
          ! 1/3 of H2 enthalpy formation goes To kinetic energy only
          ! Maybe should be corrected else-where ???
          !        up_DE(IP1) = up_DE(IP1) + creation_DE
          ! + (24 I 2002) account for lost of kinetic energy by H attaching to grain
          !   (proposed by David)

          if (react(i)%type == "H2_FO") then
            if (ikinH2 == 1) then
              up_DE(IP1) = up_DE(IP1) + 0.5_DP * (creation_DE - creation_rate * H2_int_E * kB / EVerg) &
              - creation_rate * 1.5_DP * kB * speci(IR2)%temperature / EVerg
            else
              up_DE(IP1) = up_DE(IP1) + MINIMUM(creation_DE / 3.0_DP, creation_DE - creation_rate * H2_int_E * kB / EVerg) &
              - creation_rate * 1.5_DP * kB * speci(IR2)%temperature / EVerg
            endif
          else
            up_DE(IP1) = up_DE(IP1) + creation_DE
          endif
       ELSE
          ! two or more products
          up_DE(IP1) = up_DE(IP1) + creation_DE / Nprod_m1 * &
                      (1.0_DP-speci(IP1)%mass/Mass_prod)
       ENDIF

       ! second product
       IF (IP2 > 0) THEN
          up_N(IP2)   = up_N(IP2)   + creation_rate
          up_MV(IP2)  = up_MV(IP2)  + creation_Vcm
          up_MV2(IP2) = up_MV2(IP2) + creation_Vcm2
          ! in this case, the reaction has at least two products !
          up_DE(IP2)  = up_DE(IP2)  + creation_DE/Nprod_m1 * &
                       (1.0_DP-speci(IP2)%mass/Mass_prod)
       END IF
       ! third product
       IF (IP3 > 0) THEN
          up_N(IP3)   = up_N(IP3)   + creation_rate
          up_MV(IP3)  = up_MV(IP3)  + creation_Vcm
          up_MV2(IP3) = up_MV2(IP3) + creation_Vcm2
          ! in this case, the reaction has at least three products !
          up_DE(IP3)  = up_DE(IP3)  + creation_DE / Nprod_m1 * &
                       (1.0_DP-speci(IP3)%mass/Mass_prod)
       END IF
       ! fourth product
       IF (IP4 > 0) THEN
          up_N(IP4)   = up_N(IP4)   + creation_rate
          up_MV(IP4)  = up_MV(IP4)  + creation_Vcm
          up_MV2(IP4) = up_MV2(IP4) + creation_Vcm2
          ! in this case, the reaction has four products !
          up_DE(IP4)  = up_DE(IP4)  + creation_DE / Nprod_m1 * &
                       (1.0_DP - speci(IP4)%mass / Mass_prod)
       END IF

       !-------Compute heat_exchange btw fluids in erg s-1 cm-3 (Tabone 09-18)---
       IF( react(i)%type /= "PHELE" .AND.&
           react(i)%type /= "GRREC" .AND.&
           react(i)%type /= "ADSOR" .AND.&
           react(i)%type /= "PHDES" .AND.&
           react(i)%type /= "CR_DE" .AND.&
           react(i)%type /= "SECDE" .AND.&
           react(i)%type /= "THDES" .AND.&
           react(i)%type /= "EROSI" .AND.&
           react(i)%type /= "SPUTT" ) THEN

           sum_enth_reac  = (react(i)%Nreac_3F(1)*Tn+react(i)%Nreac_3F(2)*Ti+react(i)%Nreac_3F(3)*Te) &
                          / (react(i)%Nprod_3F(1)+react(i)%Nprod_3F(2)+react(i)%Nprod_3F(3)) ! total enthalpy available per products
    
           B_heat_exchange_n   = B_heat_exchange_n   + GAMMA1 * kB * creation_rate &
                                * ( react(i)%Nprod_3F(1)*sum_enth_reac - react(i)%Nreac_3F(1)*Tn )
           
           B_heat_exchange_i   = B_heat_exchange_i   + GAMMA1 * kB * creation_rate &
                                * ( react(i)%Nprod_3F(2)*sum_enth_reac - react(i)%Nreac_3F(2)*Ti )
            
           B_heat_exchange_neg = B_heat_exchange_neg + GAMMA1 * kB * creation_rate &
                            * ( react(i)%Nprod_3F(3)*sum_enth_reac - react(i)%Nreac_3F(3)*Te )
       ELSE
          sum_enth_reac       = 0.0_dp
          B_heat_exchange_n   = 0.0_dp
          B_heat_exchange_i   = 0.0_dp
          B_heat_exchange_neg = 0.0_dp
       ENDIF

       !--------Tabone  display main chemical exo/endothermic reaction on terminal-------
       IF      (creation_DE >= creation_max) THEN
          indexR1exo   = IR1
          indexR2exo   = IR2
          creation_max = creation_DE
       ELSE IF (creation_DE <= creation_min) THEN
          indexR1endo  = IR1
          indexR2endo  = IR2
          creation_min = creation_DE
       ENDIF

       !-----------------------------------------------------------------------------
       ! Tabone: write all reaction with its heating/cooling term (in erg/s/cm3)
       ! only for debug purposes 
       !-----------------------------------------------------------------------------
       ! WRITE(*,*) &
       !     speci(IR1)%name, '+', speci(IR2)%name, '+', speci(IR3)%name, '->',  
       !     speci(IP1)%name, '+', speci(IP2)%name, '+', speci(IP3)%name, 'heating', creation_DE*EVerg
       !-----------------------------------------------------------------------------

    END DO ! end of the loop on chemical reactions


    !=====================================================================================
    ! 2 - Calculate for every chemical species the change per unit volume and time in :
    !       YN(i) -> number (cm-3.s-1)
    !       YS(i) -> mass (g.cm-3.s-1)
    !       YA(i) -> momentum (g.cm-2.s-2)
    !       YB(i) -> energy (erg.cm-3.s-1)
    !=====================================================================================
    ! Compute the derivatives using sorted reaction rates or not
    IF ( F_SORT == 1 ) THEN
       YNQ(0) = Zero
       DO j = 1, Nspec
          YNQ(j) = 0.0_dp
          idxd = speci(j)%nbreac_d
          idxf = speci(j)%nbreac_f
          DO WHILE (idxd > 0 .AND. idxf > 0)
             IF ( react(speci(j)%sort_reac_f(idxf))%rate > react(speci(j)%sort_reac_d(idxd))%rate ) THEN
                YNQ(j) = YNQ(j) + react(speci(j)%sort_reac_f(idxf))%rate
                idxf = idxf - 1
             ELSE
                YNQ(j) = YNQ(j) - react(speci(j)%sort_reac_d(idxd))%rate
                idxd = idxd - 1
             ENDIF
          ENDDO
          DO WHILE (idxf > 0)
             YNQ(j) = YNQ(j) + react(speci(j)%sort_reac_f(idxf))%rate
             idxf = idxf - 1
          ENDDO
          DO WHILE (idxd > 0)
             YNQ(j) = YNQ(j) - react(speci(j)%sort_reac_d(idxd))%rate
             idxd = idxd - 1
          ENDDO
       ENDDO
    ELSE
       YNQ(0) = Zero
       YNQ(1:Nspec) = up_N(1:Nspec) - down_N(1:Nspec)
    ENDIF
    YN(1:Nspec) = YNQ(1:Nspec)

    ! ---------------------------------------------------------
    ! Include conservation equations
    ! ---------------------------------------------------------
    ! 1 - in principle, this manipulation results only
    !     in a reordering of the terms and is therefore
    !     not necessary
    ! 2 - in practice, it insures conservation of hydrogen
    !     even when we include adsorption reactions where
    !     hydrogen atoms are artificially produced by the
    !     reaction network itself
    ! ==> ADDED - BG2016
    !     It would be best to write a consistant chemical 
    !     network (e.g. correct the adsorption processes) 
    !     from the start which conserve matters. In this 
    !     case, it would then be best to call CONSISTENCY 
    !     procedure in mhd_vode to try to insure conservation 
    !     of elements.
    ! ---------------------------------------------------------
    !!! DO i = 1, Nspec
    !!!    k = speci(i)%icon
    !!!    !!!! To only apply this treatment to hydrogen conservation
    !!!    !!!! uncomment the following line and comment the next one
    !!!    IF(k /= 1) CYCLE
    !!!    !!!! IF(k == 0) CYCLE
    !!!    dum = 0.0_dp
    !!!    DO j = 1, Nspec
    !!!       IF(j == i) CYCLE
    !!!       dum = dum + YN(j) * speci(j)%formula(k)
    !!!    ENDDO
    !!!    YN(i) = -dum / speci(i)%formula(k)
    !!! ENDDO
    !!! ! charge
    !!! YN(Nspec) = SUM(YN(1:Nspec-1) * speci(1:Nspec-1)%charge)

    ! MAUVAISES DERIVEES EN QUANTITE DE MOUVEMENT ET EN ENERGIE
    ! D(N * m * V)/DT = DN/DT * m * V + N * m * DV/DT (SEUL LE PREMIER TERME EST COMPTE ICI)
    ! D(N * m * V2)/DT = DN/DT * m * V2 + N * m * 2*V*DV/DT (SEUL LE PREMIER TERME EST COMPTE ICI)
    ! POUR LES GRAINS, G+, G-, G0, IL FAUT AUSSI COMPTER LE FAIT QUE LA MASSE CHANGE
    ! D(N * m * V)/DT = DN/DT * m * V + N * Dm/DT * V + N * m * DV/DT
    ! D(N * m * V2)/DT = DN/DT * m * V2 + N * Dm/DT * V2 + N * m * 2*V*DV/DT

    YS(1:Nspec) = YN(1:Nspec) * speci(1:Nspec)%mass
    YA(1:Nspec) = (up_MV(1:Nspec) - down_MV(1:Nspec)) * speci(1:Nspec)%mass
    YB(1:Nspec) = (up_MV2(1:Nspec) - down_MV2(1:Nspec)) * speci(1:Nspec)%mass
    YA_CH = YA(ind_CH)
    YS_CH = YS(ind_CH)
    YA_S  = YA(ind_S)
    YS_S  = YS(ind_S)
    YA_SH = YA(ind_SH)
    YS_SH = YS(ind_SH)


    !=====================================================================================
    ! 3 - Calculate the source terms (the changes per unit volume and time)
    !     of each fluid by summing over the source terms for the chemical species :
    !       Cup    -> increase in number (cm-3.s-1)
    !       Cdown  -> decrease in number (cm-3.s-1)
    !       YN     -> net change in number (cm-3.s-1) [1], eq. (16,17,20)
    !       S      -> mass (g.cm-3.s-1)               [1], eq. (18,19)
    !       A      -> momentum (g.cm-2.s-2)           [1], eq. (21)
    !       B      -> energy (erg.cm-3.s-1)           [1], eq. (25,26,31)
    !=====================================================================================

    ! neutral species
    CupN   = SUM(up_N(b_neu:e_neu))
    CdownN = SUM(down_N(b_neu:e_neu))
    YNn    = SUM(YN(b_neu:e_neu))
    Sn     = SUM(YS(b_neu:e_neu))
    An     = SUM(YA(b_neu:e_neu))
    Bn     = 0.5_DP * SUM(YB(b_neu:e_neu)) + SUM(up_DE(b_neu:e_neu)) * EVerg

    ! positive ions
    CupI   = SUM(up_N(b_ion:e_ion))
    CdownI = SUM(down_N(b_ion:e_ion))
    YNi    = SUM(YN(b_ion:e_ion))
    Si     = SUM(YS(b_ion:e_ion))
    Ai     = SUM(YA(b_ion:e_ion))
    Bi     = 0.5_DP * SUM(YB(b_ion:e_ion)) + SUM(up_DE(b_ion:e_ion)) * EVerg

    ! anions
    CupA   = SUM(up_N(b_ani:e_ani))
    CdownA = SUM(down_N(b_ani:e_ani))
    YNa    = SUM(YN(b_ani:e_ani))
    Sa     = SUM(YS(b_ani:e_ani))
    Aa     = SUM(YA(b_ani:e_ani))
    Ba     = 0.5_DP * SUM(YB(b_ani:e_ani)) + SUM(up_DE(b_ani:e_ani)) * EVerg

    ! negative species
    CupNEG   = SUM(up_N(b_neg:e_neg))
    CdownNEG = SUM(down_N(b_neg:e_neg))
    YNneg    = SUM(YN(b_neg:e_neg))
    Sneg     = SUM(YS(b_neg:e_neg))
    Aneg     = SUM(YA(b_neg:e_neg))
    Bneg     = 0.5_DP * SUM(YB(b_neg:e_neg)) + SUM(up_DE(b_neg:e_neg)) * EVerg

  END SUBROUTINE CHEMISTRY



  SUBROUTINE SOURCE

    !---------------------------------------------------------------------------
    ! called by :
    !     DIFFUN
    ! purpose :
    !     First calculates the cooling of the different fluids, then computes
    !     the source terms needed in DIFFUN: momentum and energy. These terms
    !     are updated from their values calculated in CHEMISTRY. Note that
    !     source terms for number and mass density are not modified.
    !
    !     Source takes into account the cooling due to H2 lines, fine structure
    !     lines and other molecular cooling.
    !     It includes also other physical processes :
    !        - heat transfer between neutral and charged fluids;
    !        - photo-electric effect on grains;
    !        - thermalisation with grains;
    !        - momentum and energy transfer between neutral and charged fluids
    !          through elastic scattering;
    !        - energy transfer between charged fluids through
    !          Coulomb scattering.
    ! subroutine/function needed :
    !     LINE_THERMAL_BALANCE
    !     COMPUTE_MOLECULAR_COOLING
    !     COMPUTE_H2
    ! input variables :
    ! ouput variables :
    ! results :
    !     An,  Ai,  Aneg  -> source terms for momentum (g.cm-2.s-2)
    !     Bn,  Bi,  Bneg  -> source terms for energy (erg.cm-3.s-1)
    ! references :
    !     [1]  Flower et al., 1985, MNRAS 216, 775
    !     [2]  Flower et al., 1986, MNRAS 218, 729
    !     [3]  Pineau des Forets et al., 1986, MNRAS 220, 801
    !     [4]  Monteiro & Flower, 1987, MNRAS 228, 101
    !     [5]  Black, 1987, in Interstellar processes, eds. Hollenbach & !
    !          Thronson (Reidel, Dordrecht), P.731
    !     [6]  Viala et al., 1988, A&A 190, 215
    !     [7]  Spitzer, 1962, Physics of fully ionised gases (Wiley, New York)
    !     [8]  Spitzer & Scott, 1969, AP.J. 158, 161
    !     [9]  Flower & Watt, 1984, MNRAS 209, 25
    !     [10] Roueff 1990, A&A 234, 567
    ! modification history :
    !   - nov 2000 : David Wilgenbus (suggestion of David Flower)
    !       correction of the source term for momentum transfer between grains
    !       and neutrals.
    !       "A_grain_n = RhoN * 1.1D-21 * nH * ABS_DeltaV * DeltaV"
    !       is replaced by
    !       "A_grain_n = RhoN*dens_grc *pi*rsq_grm * ABS_DeltaV*DeltaV"
    !   - sept 2018 : Tabone (+BG + GPdF). Generalize Eq. (27)-(28)-(29) in [1]
    !              to include any type of reaction in the heat_exchange terms
    !               B_heat_exchange_n , B_heat_exchange_i, B_heat_exchange_neg are
    !               now computed in CHEMISTRY to treat "heat exchange" consistently
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_PHYS_VAR
    USE MODULE_VAR_TH_BALANCE
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_H2
    USE MODULE_GRAINS, ONLY : Tgrain, dens_grc, rsq_grm
    USE MODULE_DEBUG_JLB
    USE MODULE_LINE_EXCIT
    USE MODULE_CHEM_REACT                                 ! Tabone 01/18 photochemical heating
    USE MODULE_RADIATION, ONLY : fluph, fluph0, fluphfree

    IMPLICIT NONE
    REAL(KIND=DP) :: cooling_n, cooling_i, cooling_neg ! cooling rate (erg/cm3/s)
    REAL(KIND=DP) :: Teff, Ve, alpha_n
    REAL(KIND=DP) :: muN_plus_muI, muN_plus_muA, RhoN_times_RhoI, RhoN_times_RhoA
    REAL(KIND=DP) :: Sigma_IN, Sigma_eN, Sigma_AN, Sigma_inelastic_eN
    REAL(KIND=DP) :: Sigma_CH_n, Sigma_S_n, Sigma_SH_n
    REAL(KIND=DP) :: A_i_n, A_e_n, A_a_n, A_grpos_n, A_grneg_n
    REAL(KIND=DP) :: B_DE_n, B_DE_i, B_DE_neg
    REAL(KIND=DP) :: B_i_n, B_a_n, B_e_n, B_inelastic_e_n, B_grpos_n, B_grneg_n
    REAL(KIND=DP) :: B_therm_grain_n, B_ioniz_RC_n
    REAL(KIND=DP) :: B_n_phchem    ! Tabone 01/18 photochemical heating vector containing heating rate for each computed photoreactions with cross sections
    REAL(KIND=DP) :: lambda, B_i_e, B_i_a
    REAL(KIND=DP) :: dum

    REAL(KIND=DP) :: Sel_net_n_H2_n
    REAL(KIND=DP) :: Sel_net_neg_H2_e
    REAL(KIND=DP) :: Sel_net_n_H2_i
    REAL(KIND=DP) :: Sel_net_i_H2_i
    REAL(KIND=DP) :: cooling_i_rvH2_i
    REAL(KIND=DP) :: cooling_n_rvH2_i
    REAL(KIND=DP) :: cooling_neg_rvH2_e
    REAL(KIND=DP) :: cooling_neg_rCO_e
    REAL(KIND=DP) :: B_inelastic_n_e

    !----------------------------------------------------------
    ! save the heating / cooling rates induced by the energy
    ! defect of any reaction (computed in CHEMISTRY)
    !----------------------------------------------------------
    B_DE_n   = Bn
    B_DE_i   = Bi
    B_DE_neg = Bneg

    !----------------------------------------------------------
    ! in FINE_STRUCTURE_COOLING, cooling rates (erg/cm3/s) for
    ! fine structure lines of C+, Si+, C, Si, O, S+, and N+
    ! are calculated.
    ! results are :
    ! cooling_Cplus, cooling_Siplus, cooling_C, cooling_Si
    ! cooling_O, cooling_Splus, cooling_Nplus, cooling_Cplus_e,
    ! cooling_Siplus_e, cooling_Splus_e, cooling_Nplus_e
    !----------------------------------------------------------
    ! CALL FINE_STRUCTURE_COOLING

    ! JLB - september 2002
    ! FINE_STRUCTURE_COOLING is obsolete and superseded by :

    CALL LINE_THERMAL_BALANCE

    !---------------------------------------------------------
    ! in COMPUTE_MOLECULAR_COOLING, cooling rates (erg/cm3/s)
    ! 13CO, OH, NH3, CO and H2O are calculated.
    ! useful results for SOURCE is :
    !    molec_cool
    !---------------------------------------------------------
    !
    ! 29 mai 2001 - JLB : compute molec_cool directly in DIFFUN
    !
    ! CALL COMPUTE_MOLECULAR_COOLING

    !-----------------------------------------------------
    ! Fills the radiative pumping matrix
    ! in COMPUTE_H2, cooling rate (erg/cm3/s) for
    ! ro-vibrational lines of H2 and source terms for
    ! H2 levels are calculated.
    ! results are :
    ! Cooling_H2, YN_rovib_H2, H2_Energy
    !-----------------------------------------------------
    CALL COMPUTE_H2

    ! analytic cooling supersedes all other cooling processes.
    IF (cool_KN==2) THEN
       RETURN
    ENDIF
    !------------------------------------------------------------------
    ! Add different contributions to the cooling of the neutral fluid
    ! cooling_n (erg/cm3/s) :
    ! (1) fine structure lines (from FINE_STRUCTURE_COOLING)
    !     remark : we assume that, in collisions with C+ or Si+, 
    !              the kinetic energy is lost by the (light) neutrals, 
    !              not by the (heavy) ions; this is why the terms  
    !              cooling_Cplus and cooling_Siplus are here and
    !              not in the cooling of the positively charged fluid.
    ! (1) - bis : Superseded by LINE_THERMAL_BALANCE
    !     For each collision, energy loss is allocated accurately to
    !     each fluid.
    ! (2) molecular (not H2) cooling (from COMPUTE_MOLECULAR_COOLING)
    ! (3) H2 rovibrational cooling (from COMPUTE_H2)
    ! (4) dissociation of H2 by collisions with H, He and H2
    !     we assume a rate 8 times smaller for H2 than for H
    !     and 10 times smaller for He
    !     remark :
    !        this process is not taken into account in CHEMISTRY as DE=0._DP
    !        for endothermic reactions (see ENERGY_DEFECT_REACTION).
    !------------------------------------------------------------------
    cooling_n = Zero
    ! (1)
    cooling_n = cooling_n - Cool_n
    ! (2)
    !
    ! 29 mai 2001 - JLB : added directly in DIFFUN
    !
    ! cooling_n = cooling_n + molec_cool
    !
    ! (3)
    cooling_n = cooling_n + cooling_H2
    ! (4) to be supressed if endothermic reaction don't have DE = Zero
    !     Use real collisional dissociation by neutrals 
    !     (56641.590 is energy of highest H2 level in K)
    ! NOTE - DE is automatically set to 0 for all collisional dissociation
    !        reactions of H2 (see REACTION_TYPE in chemical reaction.f90)
    Sel_net_n_H2_n = SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_ne_H2 &
                   * (56641.590-H2_lev(1:NH2_lev)%energy)) * kB

    cooling_n = cooling_n - Sel_net_n_H2_n

    ! cooling_n = cooling_n + &
    !      Dens_H * Dens_H2 * 7.2D-22*TEXP(-52000._DP/Tn) + &
    !      Dens_H2* Dens_H2 * 9.0D-23*TEXP(-52000._DP/Tn)

    !------------------------------------------------------------------------
    ! Add the different contributions to the rate of radiative cooling of the
    ! negatively charged fluid cooling_neg (erg/cm3/s) through excitation by
    ! electron collisions of :
    ! (1) O and C : 3P -> 1D; N 4S -> 2D (MENDOZA, 1983, PNE SYMPOSIUM)
    ! (2) C+, Si+, S+, and N+ (from FINE_STRUCTURE_COOLING)
    ! (1) and (2) - bis : this done in LINE_THERMAL_BALANCE now
    ! (3) rovibrational excitation of H2, rotational excitation of CO
    ! 
    ! (4) Add the cooling of the negatively charged fluid by dissociation
    !     of H2 through collision with e-
    !-------------------------------------------------------------------------
    ! (1) and (2)
    cooling_neg = Zero
    cooling_neg = cooling_neg - Cool_e
    ! (3)
    Ve = SQRT(8.0_DP*kB*Te/(pi*me)) ! thermal velocity of e- (cm/s)
    !-------------------------------------------------------------------------
    !--- WARNING: NO, at high density, needs to take into account collisional
    !---          desexcitation of H2 by e => should be done in H2.f90 Tabone
    !-------------------------------------------------------------------------
    cooling_neg_rvH2_e = Dens_e *  Dens_H2 * EVerg * Ve * (4.0D-18*TEXP(-510._DP/Te)*0.044_DP &
                       + 4.0D-17*TEXP(-6330._DP/Te)*0.544_DP)
    cooling_neg_rCO_e  = Dens_e * Dens_CO * EVerg * Ve * 1.0D-15*TEXP(-5.0_DP/Te)*0.00043_DP
    ! (4)
    Sel_net_neg_H2_e   = SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_el_H2 &
                       * (116300.0_DP-H2_lev(1:NH2_lev)%energy)) * kB

    cooling_neg = cooling_neg + cooling_neg_rvH2_e + cooling_neg_rCO_e - Sel_net_neg_H2_e

    !-------------------------------------------------------------------------
    ! Add the different contributions to the rate of radiative cooling of the
    ! positively charged fluid, cooling_i (erg/cm3/s) :
    ! (1) ro-vibrational excitation of H2 by collisions with ions
    ! (2) dissociation of H2 by collisions with ions
    !     remark :
    !        this is not taken into account in CHEMISTRY as DE=Zero
    !        for endothermic reactions (see ENERGY_DEFECT_REACTION).
    ! (3) ionization of H2 by collisions with ions
    ! The energy loss by the positively charged fluid is proportional to
    !     mass_H2/(mass_H2 + muI); the corresponding loss by the neutral
    ! fluid is proportional to
    !     muI/(mass_H2 + muI)
    ! (4) ionization of H by collisions with ions
    ! The energy loss by the positively charged fluid is proportional to
    !     mass_H/(mass_H + muI); the corresponding loss by the neutral
    ! fluid is proportional to
    !     muI/(mass_H + muI)
    !-------------------------------------------------------------------------
    cooling_i = Zero
    ! effective temperature, used in (1), (2) and (3)
    Teff = ( mass_H2 * muI * ABS_DeltaV * ABS_DeltaV / 3.0_DP / kB &
           + (mass_H2*Ti + muI*Tn) ) / (mass_H2 + muI)
    !-------------------------------------------------------------------------
    ! Collisions with ions not taken into account in H2.f90 - adopt an
    ! analytical formula to describe those.
    !-------------------------------------------------------------------------
    ! (1)
    cooling_i_rvH2_i = ABS_DeltaV * EVerg * Dens_H2 * DensityI * &
                      ( 4.0D-16*TEXP( -510._DP/Teff)*0.044_DP &
                      + 4.0D-17*TEXP(-6330._DP/Teff)*0.544_DP )
    cooling_i_rvH2_i = cooling_i_rvH2_i * mass_H2 / (mass_H2 + muI)
    cooling_n_rvH2_i = cooling_i_rvH2_i * muI / mass_H2

    cooling_i = cooling_i + cooling_i_rvH2_i
    cooling_n = cooling_n + cooling_n_rvH2_i

    ! (2) to be supressed if endothermic reaction don't have DE = Zero
    Sel_net_i_H2_i = SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_io_H2 &
                   * (56641.590-H2_lev(1:NH2_lev)%energy)) * kB
    Sel_net_i_H2_i = Sel_net_i_H2_i * mass_H2 / (mass_H2 + muI)
    Sel_net_n_H2_i = Sel_net_i_H2_i * muI / mass_H2
    
    cooling_i = cooling_i - Sel_net_i_H2_i
    cooling_n = cooling_n - Sel_net_n_H2_i

    ! ----------------------------------------------------------------
    ! (3) to be supressed if endothermic reaction don't have DE = Zero
    ! ----------------------------------------------------------------
    ! removed. We now take into account endothermic reaction - BG2017
    ! ----------------------------------------------------------------
    ! rate = DensityI * Dens_H2 * 1.10D-13*(Teff/300._DP)**0.5_DP* &
    !        TEXP(-179160.0_DP/Teff)*15.43_DP * EVerg
    ! cooling_i = cooling_i + rate * mass_H2 / (mass_H2 + muI)
    ! cooling_n = cooling_n + rate * muI / (mass_H2 + muI)
    ! ----------------------------------------------------------------
    ! (4) to be supressed if endothermic reaction don't have DE = Zero
    ! ----------------------------------------------------------------
    ! removed. We now take into account endothermic reaction - BG2017
    ! ----------------------------------------------------------------
    ! Teff = & ! effective temperature used in (4)
    !      (mass_H*muI*ABS_DeltaV*ABS_DeltaV/3.0_DP/kB + (mass_H*Ti + muI*Tn)) / &
    !      (mass_H + muI)
    ! rate = DensityI * Dens_H * 1.30D-13*(Teff/300._DP)**0.5_DP* &
    !        TEXP(-157890.0_DP/Teff)*13.60_DP * EVerg
    ! cooling_i = cooling_i + rate * mass_H/(mass_H + muI)
    ! cooling_n = cooling_n + rate * muI/(mass_H + muI)

    ! (0) add contribution from LINE_THERMAL_BALANCE
    ! Added after (1)-(4) so as not to add it also to cooling_n
    cooling_i = cooling_i - Cool_i

    !--- some useful variables ---
    muN_plus_muI    = muN  + muI  ! sum of neutral and positive ion mass (g)
    muN_plus_muA    = muN  + muA  ! sum of neutral and negative ion mass (g)
    RhoN_times_RhoI = RhoN * RhoI ! product of neutral and positive ion mass densisities (g2/cm6)
    RhoN_times_RhoA = RhoN * RhoA ! product of neutral and negative ion mass densisities (g2/cm6)

    !--- polarisability of the neutrals                    ---
    !--- = weighted average of polarisability of H and H2  ---
    alpha_n = ( Dens_H*alpha_H + Dens_H2*alpha_H2 ) / ( Dens_H + Dens_H2 )

    !------------------------------------------------------------------
    ! calculate the cross-sections (cm-2)
    ! (1) Sigma_IN           -> elastic ions >0  - neutral scattering
    !                           ([1], eq.(23))
    ! (2) Sigma_eN           -> elastic electron - neutral scattering
    !                           ([1], eq.(34))
    ! (3) Sigma_AN           -> elastic ions <0  - neutral scattering
    !                           ([1], eq.(23))
    ! (4) Sigma_inelastic_eN -> inelastic electron-neutral scattering
    ! (5) Sigma_CH_n -> elastic coupling between CH and neutrals
    ! (6) Sigma_S_n  -> elastic coupling between S  and neutrals
    ! (7) Sigma_SH_n -> elastic coupling between SH and neutrals
    !------------------------------------------------------------------
    ! (1)
    Sigma_IN = 2.41_DP * pi * qe * SQRT ( muN_plus_muI * alpha_n / muN / muI)
    Sigma_IN = MAXIMUM(Sigma_IN, 1.0D-15*ABS_DeltaV)
    ! (2)
    Sigma_eN = 1.0D-15 * SQRT( 8.0_DP * kB * Te / ( pi * me ) )
    ! (3)
    Sigma_AN = 2.41_DP * pi * qe * SQRT ( muN_plus_muA * alpha_n / muN / muA )
    Sigma_AN = MAXIMUM(Sigma_AN, 1.0D-14 * ABS_DeltaV)
    ! (4)
    Sigma_inelastic_eN = 1.0D-16 * ( 4.29D-1 + 6.14D-6*Te ) * SQRT( 8.0_DP * kB * Te / ( pi * me ) )
    ! (5) !CH! Flower Pineau des Forets (1998) MNRAS, 297, 1182
    ! Be careful: The reference given by GPDF (1998) for the collision rate of CH and H2
    ! is Zabarnick et al. (1986). However the value given in this paper is that for CH2
    ! + H and not CH + H2.
    Teff = mass_CH * muN * (Vn-v_CH)**2 / (muN + mass_CH)
    Teff = Tn + Teff / (3.*kB)
    Sigma_CH_n = 5.00D-10 * (Teff/300.)**0.25
    ! (6) !S!  same as CH
    Teff = mass_S * muN * (Vn-v_S)**2 / (muN + mass_S)
    Teff = Tn + Teff / (3.*kB)
    Sigma_S_n = 5.00D-10 * (Teff/300.)**0.25
    ! (7) !SH! same as CH
    Teff = mass_SH * muN * (Vn-v_SH)**2 / (muN + mass_SH)
    Teff = Tn + Teff / (3.*kB)
    Sigma_SH_n = 5.00D-10 * (Teff/300.)**0.25

    !--------------------------------------------------------------------
    ! Rates of production of momentum (g.cm-2.s-2) in the neutral fluid
    !   A_i_n     -> elastic ions >0  - neutral scattering
    !   A_e_n     -> elastic electron - neutral scattering
    !   A_a_n     -> elastic ions <0  - neutral scattering
    !   A_grpos_n -> elastic positive charged grains - neutral scattering
    !   A_grneg_n -> elastic negative charged grains - neutral scattering
    !--------------------------------------------------------------------
    A_i_n     = Sigma_IN * DeltaV * RhoN_times_RhoI / muN_plus_muI !/fudge
    A_e_n     = DensityN * Dens_e * me * Sigma_eN * DeltaV
    A_a_n     = Sigma_AN * DeltaV * RhoN_times_RhoA / muN_plus_muA
    ! grain - neutral scatt -> G0 and G+ are assumed to belong to the positive fluid
    !                       -> G- to the negative fluid
    A_grpos_n = RhoN * dens_grc * pi * rsq_grm * ABS_DeltaV * DeltaV &
              * (Dens_G + Dens_Gplus)/ (Dens_Gminus + Dens_Gplus + Dens_G)
    A_grneg_n = RhoN * dens_grc * pi * rsq_grm * ABS_DeltaV * DeltaV &
              * (Dens_Gminus)/ (Dens_Gminus + Dens_Gplus + Dens_G)
    !CH!
    A_CH_n = RhoN * Dens_CH * mass_CH/(muN+mass_CH)*Sigma_CH_n*(Vn-v_CH) 
    !S!
    A_S_n  = RhoN * Dens_S  * mass_S /(muN+mass_S )*Sigma_S_n *(Vn-v_S ) 
    !SH!
    A_SH_n = RhoN * Dens_SH * mass_SH/(muN+mass_SH)*Sigma_SH_n*(Vn-v_SH) 

    !------------------------------------------------------------------
    ! Source terms of number and mass density are calculated in
    ! CHEMISTRY and not in SOURCE.
    ! These terms are (for neutrals, positive fluid, negative fluid) :
    ! YNn, YNi, YNneg -> change in number density (cm-3.s-1)
    ! Sn,  Si,  Sneg  -> change in mass density (g.cm-3.s-1)
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Source terms of momentum (g.cm-2.s-2) of each fluid (update of
    ! the calculation in CHEMISTRY).
    !   An   -> neutral fluid
    !   Ai   -> positive fluid
    !   Aneg -> negative fluid
    !----------------------------------------------------------------
    An   = An   + A_i_n + A_e_n + A_a_n + A_grpos_n + A_grneg_n
    Ai   = Ai   - A_i_n                 - A_grpos_n
    Aneg = Aneg         - A_e_n - A_a_n             - A_grneg_n

    !--------------------------------------------------------------------------
    ! Source term of energy Bn (erg.cm-3.s-1) in the neutral fluid
    ! obtained by adding up :
    ! +Bn                 -> exo/endothermicity of the reactions
    !                        (from CHEMISTRY)
    ! -cooling_n          -> energy loss through cooling
    ! +B_i_n              -> elastic diffusion of the neutrals by the
    !                        positive ions ([1], eq.(32))
    ! +B_a_n              -> elastic diffusion of the neutrals by the
    !                        anions ([1], eq.(32))
    ! +B_e_n              -> elastic diffusion of the neutrals by the
    !                        electrons ([1], eq.(33))
    ! +B_inelastic_e_n    -> inelastic scattering of the electrons on H2;
    !                        mean excitation energy 12 eV = 140 000 K for the
    !                        singlet excited states, 10 eV = 116 300 K for the
    !                        triplet, of which 10 - 4.48 = 5.52 eV are recovered
    !                        by the neutrals through the dissociation products;
    !                        collisional excitation of the N=2 state of H.
    ! +B_grpos_n          -> elastic diffusion of the neutrals by the
    !                        positively charged grains
    ! +B_grneg_n          -> elastic diffusion of the neutrals by the
    !                        negatively charged grains
    ! +B_heat_exchange_n  -> heat exchange
    !                        [1],eq.(27) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1) + computed in CHEMISTRY (BT 09-18)
    ! +B_photo_grain_n    -> heating through the photo-electric effect on grains
    !                        [5]
    ! +B_therm_grain_n    -> energy loss/gain through thermalisation with grains
    !                        TIELENS  & HOLLENBACH (1985, AP.J. 291, 722)
    ! +B_ioniz_RC_n       -> heating through cosmic ray ionization [8]
    !                        (not included in CHEMISTRY as DE = Zero for
    !                        reactions of the type 'CR_IO')
    !--------------------------------------------------------------------------
    
    B_i_n = ( 3.0_DP * kB * (Ti-Tn) + DeltaV * (muI * Vi + muN * Vn) ) &
          * RhoN_times_RhoI / muN_plus_muI * Sigma_IN / muN_plus_muI

    B_a_n = ( 3.0_DP * kB * (Te-Tn) + DeltaV * (muA * Vi + muN * Vn) ) & 
          * RhoN_times_RhoA / muN_plus_muA * Sigma_AN / muN_plus_muA

    B_e_n = ( 4.0_DP * kB * (Te-Tn) + muN * Vn * DeltaV ) &
          * DensityN * Dens_e * me / muN * Sigma_eN

    B_inelastic_e_n = - SUM(DBLE(H2_lev(1:NH2_lev)%density) * Sel_el_H2) * 5.52_DP * EVerg

    B_grpos_n = A_grpos_n * Vi
    B_grneg_n = A_grneg_n * Vi

    ! B_photo_grain_n = 4.0D-26*nH*RAD*TEXP(-2.5_DP*Av)
    ! Correction : from Bakes and Tielens (1994) - Eqs. 42 and 43
    ! modif Tabone 01/2018 + correction 02/2018: dens_grc
    IF      (F_COUP_RAD == 1) THEN 
       dum = 4.87e-2_dp / (1.0_dp + 4e-3_dp * (RAD * sqrt(Tn) / Dens_e )**0.73_dp ) &
           + 3.65e-2_dp / (1.0_dp + 2e-4_dp * (RAD * sqrt(Tn) / Dens_e) ) * (Tn / 1e4_dp)**0.7_dp
       B_photo_grain_n = 1e-24_dp * RAD * nH * dum * (fluph/fluph0)**(2.5/2.8) * (dens_grc/nH) / 7.0e-11_dp
    ELSE IF (F_COUP_RAD == 2) THEN
       dum = 4.87e-2_dp / (1.0_dp + 4e-3_dp * (fluph / fluphfree * sqrt(Tn) / Dens_e )**0.73_dp ) &
           + 3.65e-2_dp / (1.0_dp + 2e-4_dp * (fluph / fluphfree * sqrt(Tn) / Dens_e) ) * (Tn / 1e4_dp)**0.7_dp
       B_photo_grain_n = 1e-24_dp * fluph / fluphfree * nH * dum * (dens_grc/nH)/7.0D-11
    ELSE
       dum = 4.87e-2_dp / (1.0_dp + 4e-3_dp * (RAD * sqrt(Tn) / Dens_e )**0.73_dp ) &
           + 3.65e-2_dp / (1.0_dp + 2e-4_dp * (RAD * sqrt(Tn) / Dens_e) ) * (Tn / 1e4_dp)**0.7_dp
       B_photo_grain_n = 1e-24_dp * RAD * nH * dum * TEXP(-2.5_DP*Av) * (dens_grc/nH) / 7.0e-11_dp
    ENDIF

    ! modif Tabone 01/2018 + correction 02/2018: dens_grc
    B_therm_grain_n = 3.5D-34 * SQRT(Tn) * (Tgrain - Tn) * nH**2.0_DP * (dens_grc/nH) / 7.0e-11_dp

    B_ioniz_RC_n = Zeta * EVerg * DensityN * &
         MAXIMUM(5.7_DP,32.0_DP-7.0_DP*LOG10(DensityN/Dens_e))
    !------ Tabone 01/18 photochemical heating
    !-------phdest is computed in initalise.f90 and in the difun if Av is integrated
    IF(F_COUP_RAD == 1 .OR. F_COUP_RAD == 2) THEN
        B_n_phchem = SUM(speci(phdest(1:Ncross)%ispe)%density*phheating_rate(1:Ncross)) ! photochemical heating in erg.cm-3.s-1
        ! write(*,*) "Photoheating=", heat_i_phchem, "erg cm-3 s-1"
     ELSE
        B_n_phchem = 0.0_DP
    ENDIF

    !--- sum of the different contributions to Bn ---
    Bn = B_DE_n            &
       - cooling_n         &
       + B_i_n             &
       + B_a_n             &
       + B_e_n             &
       + B_inelastic_e_n   &
       + B_grpos_n         &
       + B_grneg_n         &
       + B_photo_grain_n   &
       + B_therm_grain_n   &
       + B_ioniz_RC_n      &
       + B_n_phchem        & ! Tabone 01/18 photochemical heating
       + B_heat_exchange_n   ! Tabone add B_heat_exchange_n previously discaded for neutrals

    !--------------------------------------------------------------------------
    ! Source term of energy Bi (erg.cm-3.s-1) for the positive ions fluid
    ! obtained by adding up :
    !  +Bi                -> exo/endothermicity of the reactions
    !                        (from CHEMISTRY)
    !  -cooling_i         -> energy loss through cooling
    !  -B_i_e             -> diffusion of electrons by positive ions
    !                        [7], see also [1], eq.(35)
    !  -B_i_a             -> diffusion of negative ions by positive ions
    !                        [7]
    !  -B_i_n             -> elastic diffusion of the neutrals by the
    !                        positive ions. See the calculation for Bn above
    !  -B_grpos_n          -> elastic diffusion of the neutrals by the
    !                         positively charged grains, see calculation for Bn above.
    !  +B_heat_exchange_i -> heat exchange
    !                        [1],eq.(28) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1)
    !--------------------------------------------------------------------------
    lambda = 1.5_DP / (qe**3.0_DP) * SQRT( (kB * Te)**3.0_DP / (Dens_e * pi) )
    B_i_e = 4.0_DP * qe**4._DP / muI / Te * SQRT( 2.0_DP * me * pi / (kB * Te) ) &
          * DensityI * Dens_e * LOG(lambda)*(Ti-Te)

    lambda = 1.5_DP / (qe**3.0_DP) * SQRT( (kB * Te)**3.0_DP / (DensityI * pi) )
    B_i_a = qe**4.0_DP / (muI * muA) * (Ti / muI + Te / muA)**(-1.5_DP) &
          * (32.0_DP*pi/kB)**0.5_DP * DensityI * DensityA * LOG(lambda) * (Ti-Te)

    !--- sum of the different contributions to Bi ---
    Bi = B_DE_i             &
       - cooling_i          &
       - B_i_e              &
       - B_i_a              &
       - B_i_n              &
       - B_grpos_n          &
       + B_heat_exchange_i

    !--------------------------------------------------------------------------
    ! CALCULATE THE SOURCE TERM OF ENERGY Bneg (ERG.CM-3.S-1) FOR THE
    ! NEGATIVE FLUID BY ADDING UP
    !  +Bneg               -> exo/endothermicity of the reactions
    !                        (from CHEMISTRY)
    !   B_i_e              -> diffusion of electrons by positive ions
    !                         [7], see also [1], eq.(35).
    !                         see calculation for Bi above.
    !   B_i_a              -> diffusion of negative ions by positive ions
    !                         [7], see calculation for Bi above.
    !  -B_e_n              -> elastic diffusion of the neutrals by the
    !                         electrons ([1], eq.(33)).
    !                         see calculation for Bn above.
    !  -B_inelastic_n_e    -> inelastic scattering of the electrons on H2 and H.
    !                         For H, excitation of the N=2 state :
    !                         Aggarwal et al. 1991, J. PHYS. B, 24, 1385
    !                                ionization:
    !                         Hummer, 1963, MNRAS, 125, 461.
    !                         For H2, excitation of the repulsive B triplet state
    !                         and of the B and C singlet states
    !                                ionization:
    !                         Rapp & Englander-Golden, 1965, JCP, 43, 1464.
    !  -B_grneg_n          -> elastic diffusion of the neutrals by the
    !                         negatively charged grains, see calculation for Bn above.
    !  -B_a_n              -> elastic diffusion of the neutrals by the
    !                         negative ions ([1], eq.(32)).
    !                         see calculation for Bn above.
    !  cooling_neg         -> energy loss through cooling
    !  B_heat_exchange_neg -> heat exchange
    !                        [1],eq.(28) with the 5/2 factors corrected to 3/2
    !                        (GAMMA -> GAMMA1) + computed in CHEMISTRY for any reaction
    !--------------------------------------------------------------------------
    B_inelastic_n_e = Dens_H2 * Dens_e * Sigma_inelastic_eN * TEXP(-1.4D5/Te) &
                    * 12.0_DP * EVerg

    ! remove the following term -> H excitation now taken into account in line_excit.f90 cooling_Hat
    ! B_inelastic_n_e = B_inelastic_n_e &
    !                 + Dens_H*Dens_e*EVerg*3.8D-06*TEXP(-1.1842D5/Te)*(Te/300._DP)**(-0.5_DP)
    ! remove the following term -> collisional ionization already included in chemical cooling
    ! B_inelastic_n_e = B_inelastic_n_e &
    !                 + Dens_H*Dens_e*EVerg*13.6_DP*9.2D-10*TEXP(-1.5789D5/Te)*(Te/300.0D0)**0.5_DP &
    !                 + Dens_H2*Dens_e*EVerg*15.43_DP*1.40D-09*TEXP(-1.7916D5/Te)*(Te/300._DP)**0.5_DP

    ! I don't understand the expression below (B_heat_exchange_neg). Why adding a term 
    ! of kinetic energy while this term seems to be neglected for all other fluids ? 
    ! Beside, why not treat the heat exchange term similarly to that computed for the ions ?
    ! B_heat_exchange_neg = B_heat_exchange_neg + 0.5_DP * me * Vi * Vi * (YNi-YNneg)

    !--- sum of the different contributions to Bneg ---
    Bneg = B_DE_neg             &
         - cooling_neg          &
         + B_i_e                &
         + B_i_a                &
         - B_e_n                &
         - B_inelastic_n_e      &
         - B_grneg_n            &
         - B_a_n                &
         + B_heat_exchange_neg


    !--------------------------------------------------------------------------
    ! BT/BG 17-05-2019: save cooling/heating terms in output variables for 
    !                   ascii and HFD5
    !--------------------------------------------------------------------------

    !----Atomic cooling---- 
    ! Already computed in LINE_THERMAL_BALANCE (line_excit.f90)
    ! line_rad_n_atoms  = Cool_n ! Cool_n is negative
    ! line_rad_i_atoms  = Cool_i
    ! line_rad_e_atoms  = Cool_e

    !----Molecular cooling---
    ! => cooling of the neutrals by CO/H20/OH/NH3: 
    !    defined in COMPUTE_MOLECULAR_COOLING (molecular_cooling.f90)
    ! => there, line_rad_n_molecules is also summed
    line_rad_n_H2     = - cooling_H2       &      ! cooling_H2 is positive
                        - cooling_n_rvH2_i        ! cooling_n_rvH2_i is positive

    line_rad_i_H2     = - cooling_i_rvH2_i        ! cooling_i_rvH2_i is positive

    line_rad_e_H2     = - cooling_neg_rvH2_e &    ! cooling_neg_rvH2_e is positive
                        - B_inelastic_n_e         ! B_inelastic_n_e    is positive
    line_rad_e_CO     = - cooling_neg_rCO_e       ! cooling_neg_rCO_e  is positive

    ! To study in detail the H2 cooling, uncomment those lines and output the results
    ! h2excit_n_coh = - cooling_H2
    ! h2excit_n_add = - cooling_n_rvH2_i
    ! h2excit_i_add = - cooling_i_rvH2_i
    ! h2excit_e_add = - cooling_neg_rvH2_e
    ! h2excit_e_LW  = - B_inelastic_n_e

    !----defect energy---
    chem_DE_n_other   =  B_DE_n
    chem_DE_i_other   =  B_DE_i
    chem_DE_e_other   =  B_DE_neg

    chem_DE_n_diss_H2 =  Sel_net_n_H2_n   &       ! Sel_net_n_H2_n   is negative
                      +  Sel_net_n_H2_i   &       ! Sel_net_n_H2_i   is negative
                      +  B_inelastic_e_n          ! B_inelastic_e_n  is positive (heating by dissociation)
    chem_DE_i_diss_H2 =  Sel_net_i_H2_i           ! Sel_net_i_H2_i   is negative
    chem_DE_e_diss_H2 =  Sel_net_neg_H2_e         ! Sel_net_neg_H2_e is negative

    chem_DE_n_phgas   = B_n_phchem
    chem_DE_i_phgas   = 0._DP
    chem_DE_e_phgas   = 0._DP
    
    chem_DE_n_phgrn   = B_photo_grain_n
    chem_DE_i_phgrn   = 0._DP
    chem_DE_e_phgrn   = 0._DP

    chem_DE_n_cosmic  = B_ioniz_RC_n
    chem_DE_i_cosmic  = 0._DP
    chem_DE_e_cosmic  = 0._DP

    chem_DE_n         = chem_DE_n_other   &
                      + chem_DE_n_diss_H2 &
                      + chem_DE_n_phgas   &
                      + chem_DE_n_phgrn   &
                      + chem_DE_n_cosmic
    chem_DE_i         = chem_DE_i_other   &
                      + chem_DE_i_diss_H2 &
                      + chem_DE_i_phgas   &
                      + chem_DE_i_phgrn   &
                      + chem_DE_i_cosmic
    chem_DE_e         = chem_DE_e_other   &
                      + chem_DE_e_diss_H2 &
                      + chem_DE_e_phgas   &
                      + chem_DE_e_phgrn   &
                      + chem_DE_e_cosmic
    chem_DE_tot       = chem_DE_n         &
                      + chem_DE_i         &
                      + chem_DE_e

    elast_scat_n      =   B_i_n + B_e_n + B_a_n                    + B_grneg_n + B_grpos_n
    elast_scat_i      = - B_i_n                 - B_i_a - B_i_e                - B_grpos_n
    elast_scat_e      =         - B_e_n - B_a_n + B_i_a + B_i_e    - B_grneg_n
    elast_scat_tot    = elast_scat_n + elast_scat_i + elast_scat_e

    exch_eint_n       = B_heat_exchange_n
    exch_eint_i       = B_heat_exchange_i
    exch_eint_e       = B_heat_exchange_neg
    exch_eint_tot     = B_heat_exchange_n + B_heat_exchange_i + B_heat_exchange_neg

    therm_grain_tot   = B_therm_grain_n

  END SUBROUTINE SOURCE



  SUBROUTINE DIFFUN(N, Z, Y, DY)

    !---------------------------------------------------------------------------
    ! called by :
    !    STIFF
    !    PSET
    ! purpose :
    !    calculates DY, the partial derivatives of the Y vector.
    !               DY(i)/DZ= F(Y,Z), 1<= i <= N
    !    see [1], equations (3)-(14)
    !    The source terms (change per unit volume and time) are calculated in
    !    CHEMISTRY and SOURCE.
    !    DIFFUN also extracts all physical variables from the Y vector
    !    (this variables are used in SOURCE).
    ! subroutine/function needed :
    !    COMPUTE_OP_H2
    !    CHEMISTRY
    !    SOURCE
    ! input variables :
    !    N -> dimension of Y and DY
    !    Z -> distance from the cloud (useless here)
    !    Y -> vector containing MHD variables + species + para-H2 + H2 levels
    ! output variables :
    !    DY -> first order derivatives DY/dZ
    ! results :
    !    v_l_var, v_l_der ...
    ! references :
    !   [1] Flower et al., 1985, MNRAS 216, 775
    !---------------------------------------------------------------------------
    USE MODULE_VAR_VODE,       ONLY : integ_type
    USE MODULE_SHIELD
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_GAMMA
    USE MODULE_GRAINS
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_H2
    USE MODULE_CO
    USE MODULE_CONSTANTS,      ONLY : kB, pi, Zero
    USE MODULE_DEBUG_JLB
    USE MODULE_LINE_EXCIT,     ONLY : total_n_cool
    USE MODULE_DUST_TREATMENT, ONLY : COMPUTE_QCOEF, COMPUTE_DERO, COMPUTE_DADS
    USE MODULE_RADIATION
    USE MODULE_CHEM_REACT,     ONLY : Ncross, phdest, phdest_rate, phheating_rate ! Tabone photoheating
    
    USE tex
    USE MODULE_VAR_VODE,       ONLY : Tout_V

    IMPLICIT NONE
    INTEGER(KIND=LONG),          INTENT(in)    :: N           ! dimension of Y and DY
    REAL(KIND=DP),               INTENT(in)    :: Z           ! useless here : the dependence in Z is implicit
    REAL(KIND=DP), DIMENSION(N), INTENT(inout) :: Y           ! vector containing the log of the variables
    REAL(KIND=DP), DIMENSION(N), INTENT(inout) :: DY          ! vector containing the derivatives
    REAL(KIND=DP)                              :: Bfield2
    REAL(KIND=DP)                              :: kT_muN
    REAL(KIND=DP)                              :: Vi2
    REAL(KIND=DP)                              :: Vn2
    REAL(KIND=DP)                              :: aux1_dvn
    REAL(KIND=DP)                              :: aux2_dvn
    REAL(KIND=DP)                              :: dvn1
    REAL(KIND=DP)                              :: dvn2
    REAL(KIND=DP)                              :: gvn1
    REAL(KIND=DP)                              :: gvn2
    REAL(KIND=DP)                              :: divu_n
    REAL(KIND=DP)                              :: divu_i
    REAL(KIND=DP)                              :: corrdv
    REAL(KIND=DP)                              :: dmin
    INTEGER                                    :: ii, i
    REAL(KIND=DP)                              :: XLL2
    REAL(KIND=DP)                              :: Bne
    REAL(KIND=DP)                              :: Ane
    REAL(KIND=DP)                              :: Sne
    REAL(KIND=DP)                              :: YNne
    REAL(KIND=DP)                              :: DensityNe
    REAL(KIND=DP)                              :: RhoNe

    REAL (KIND=dp)                             :: diff_dist   ! dumy variable for a linear or logarithmic integration of z
    REAL (KIND=dp)                             :: bt_s_bd1    ! dumy variable for calculating the doppler weighted H2 column density
    REAL (KIND=dp)                             :: bt_s_bd2    ! dumy variable for calculating the doppler weighted H2 column density


    ! ===========================================================================
    ! STEP 1 - EXTRACT VARIABLES AND COMPUTE DERIVATIVES FROM THE Y TABLE
    ! ===========================================================================

    !-------------------------------------------------------------------
    ! choose linear or logarithmic integration for Y
    !-------------------------------------------------------------------
    ! v_variab(1:N) = Y(1:N)
    !-------------------------------------------------------------------
    v_variab(1:N) = EXP(Y(1:N))

    ! if v_variab <= Zero (unphysical for these variables) for this interpolation of DVODE
    ! use the last value v_l_var
    WHERE (v_variab <= Zero)
       v_variab = v_l_var
    END WHERE

    IF      (integ_type == 0) THEN
       diff_dist      = (Z-distance_old)
       v_dvariab(1:N) = v_variab(1:N) * DY(1:N)
    ELSE IF (integ_type == 1) THEN
       diff_dist      = (EXP(Z)-distance_old)
       v_dvariab(1:N) = v_variab(1:N) * DY(1:N) / EXP(Z)
    ENDIF

    !-------------------------------------------------------------------
    !  JLB test - 19 fev 2002
    !  ...meme genre
    !-------------------------------------------------------------------
    ! Removed - dvode doesn't like such treatments BG2017
    !-------------------------------------------------------------------
    ! IF (Force_I_C == 1) THEN
    !    DensityI = SUM(v_variab(bv_ion:ev_ion))
    !    IF (Z > 1.0e6) THEN
    !       Y(iv_DensityI) = log(DensityI)
    !    ENDIF
    !    v_variab(iv_DensityI) = DensityI
    ! ENDIF
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! Saveguard against round-off errors in the derivation of Vi, Ti, Te
    ! for one fluid model, we have :
    !    Vn = Vi (velocities), Tn = Ti = Te (temperatures)
    ! for two fluids model, we have :
    !    Ti = Te (temperatures)
    ! This can be checked in printouts periodically 
    ! performed in MHD
    !-------------------------------------------------------------------
    SELECT CASE (Nfluids)
    CASE (1)
       v_variab(iv_Vi) = v_variab(iv_Vn)
       v_variab(iv_Ti) = v_variab(iv_Tn)
       v_variab(iv_Te) = v_variab(iv_Tn)
    CASE (2)
       v_variab(iv_Te) = v_variab(iv_Ti)
    CASE DEFAULT
    END SELECT

    !-------------------------------------------------------------------
    ! extract physical variables (MODULE_PHYS_VAR)
    !-------------------------------------------------------------------

    !-----------------------------------------------------
    ! a -  temperatures (K)
    !-----------------------------------------------------
    Tn = v_variab(iv_Tn)
    Ti = v_variab(iv_Ti)
    Te = v_variab(iv_Te)

    !-----------------------------------------------------
    ! b - velocities (cm/s)
    !-----------------------------------------------------
    Vn = v_variab(iv_Vn)
    Vi = v_variab(iv_Vi)
    DeltaV     = Vi - Vn
    ABS_DeltaV = ABS(DeltaV)
    !-----------------------------------------------------
    ! Treat CH velocity independently if F_CH=1
    !-----------------------------------------------------
    IF(F_CH == 1) THEN
       v_CH = v_variab(iv_vCH)
    ELSE
       v_CH = v_variab(iv_Vn)
    ENDIF
    !-----------------------------------------------------
    ! Treat S  velocity independently if F_S=1
    !-----------------------------------------------------
    IF(F_S == 1) THEN
       v_S = v_variab(iv_vS)
    ELSE
       v_S = v_variab(iv_Vn)
    ENDIF
    !-----------------------------------------------------
    ! Treat SH velocity independently if F_SH=1
    !-----------------------------------------------------
    IF(F_SH == 1) THEN
       v_SH = v_variab(iv_vSH)
    ELSE
       v_SH = v_variab(iv_Vn)
    ENDIF

    grad_V = v_variab(iv_gv)

    !-----------------------------------------------------
    ! c - column densities and extinction
    !-----------------------------------------------------
    ! H  column-density
    coldens_h =v_variab(iv_nh)
    ! H2 column-density
    coldens_h2=v_variab(iv_nh2)
    ! CO column-density
    coldens_co=v_variab(iv_nco)

    ! Refresh Av according to column-density.
    ! Note : Always use N(H) + 2*N(H2) to compute Av, and not N(H2)
    !        alone as it was permitted in old PL version - BG 2017
    IF(F_AV /= 0) THEN
       Av = Av0 + (coldens_h2*2_DP + coldens_h) * inv_Av_fac
    ELSE
       Av = Av0
    ENDIF

    !-----------------------------------------------------
    ! d - number densities (cm-3)
    !   - mass densities (g/cm3)
    !   - average masses (g)
    !-----------------------------------------------------
    DensityN   = v_variab(iv_DensityN)
    DensityI   = v_variab(iv_DensityI)
    DensityA   = v_variab(iv_DensityA)
    DensityNEG = v_variab(iv_DensityNeg)

    RhoN   = v_variab(iv_RhoN)
    RhoI   = v_variab(iv_RhoI)
    RhoA   = v_variab(iv_RhoA)
    RhoNEG = v_variab(iv_RhoNEG)

    muN   = RhoN   / DensityN
    muI   = RhoI   / DensityI
    muA   = RhoA   / DensityA
    muNEG = RhoNEG / DensityNEG

    !-----------------------------------------------------
    ! e - characteristic speeds
    !   - sound speed and ion magnetosonic speed (cm.s-1)
    !-----------------------------------------------------
    Vsound  = SQRT ( Gamma * kB * Tn / muN )
    Valfven = SQRT ( (Bfield * Vs_cm / Vi)**2._DP / (4._DP * pi * RhoN) )

    !-----------------------------------------------------
    ! f - variables useful for calculation
    !-----------------------------------------------------
    kT_muN = kB * Tn / muN
    Vn2 = Vn * Vn
    Vi2 = Vi * Vi
    ! square of the magnetic field
    Bfield2 = (Vs_cm*Vs_cm) * (Bfield*Bfield) / (4.0_DP*pi*Vi2)
    XLL2 = XLL * XLL

    !-------------------------------------------------------------------
    ! density of some species (MODULE_CHEMICAL_SPECIES)
    !    used in SOURCE and in RSF
    ! NOTE : if the specy is not present, density = Zero
    !-------------------------------------------------------------------
    Dens_H = v_variab(bv_specy+ind_H)
    IF (ind_H==0) Dens_H = Zero
    Dens_H2 = v_variab(bv_specy+ind_H2)
    IF (ind_H2==0) Dens_H2 = Zero
    ! JLB test - 25 jan 2002
    IF (NH2_lev_var>1) THEN
       sum_H2 = SUM(v_variab(bv_H2_lev:ev_H2_lev))
    ELSE
       sum_H2 = Dens_H2
    ENDIF
    Dens_He = v_variab(bv_specy+ind_He)
    IF (ind_He==0) Dens_He = Zero
    Dens_O = v_variab(bv_specy+ind_O)
    IF (ind_O==0) Dens_O = Zero
    Dens_Oplus = v_variab(bv_specy+ind_Oplus)
    IF (ind_Oplus==0) Dens_Oplus = Zero
    Dens_N = v_variab(bv_specy+ind_N)
    IF (ind_N==0) Dens_N = Zero
    Dens_C = v_variab(bv_specy+ind_C)
    IF (ind_C==0) Dens_C = Zero
    Dens_S = v_variab(bv_specy+ind_S)
    IF (ind_S==0) Dens_S = Zero
    Dens_Si = v_variab(bv_specy+ind_Si)
    IF (ind_Si==0) Dens_Si = Zero
    Dens_H2O = v_variab(bv_specy+ind_H2O)
    IF (ind_H2O==0) Dens_H2O = Zero
    Dens_OH = v_variab(bv_specy+ind_OH)
    IF (ind_OH==0) Dens_OH = Zero
    Dens_CO = v_variab(bv_specy+ind_CO)
    IF (ind_CO==0) Dens_CO = Zero
    Dens_CH = v_variab(bv_specy+ind_CH)
    IF (ind_CH==0) Dens_CH = Zero
    Dens_SH = v_variab(bv_specy+ind_SH)
    IF (ind_SH==0) Dens_SH = Zero
    Dens_NH3 = v_variab(bv_specy+ind_NH3)
    IF (ind_NH3==0) Dens_NH3 = Zero
    Dens_G = v_variab(bv_specy+ind_G)
    IF (ind_G==0) Dens_G = Zero
    Dens_Cplus = v_variab(bv_specy+ind_Cplus)
    IF (ind_Cplus==0) Dens_Cplus = Zero
    Dens_Siplus = v_variab(bv_specy+ind_Siplus)
    IF (ind_Siplus==0) Dens_Siplus = Zero
    Dens_Hplus = v_variab(bv_specy+ind_Hplus)
    IF (ind_Hplus==0) Dens_Hplus = Zero
    Dens_Splus = v_variab(bv_specy+ind_Splus)
    IF (ind_Splus==0) Dens_Splus = Zero
    Dens_Nplus = v_variab(bv_specy+ind_Nplus)
    IF (ind_Nplus==0) Dens_Nplus = Zero
    Dens_Feplus = v_variab(bv_specy+ind_Feplus)
    IF (ind_Feplus==0) Dens_Feplus = Zero
    Dens_Gplus = v_variab(bv_specy+ind_Gplus)
    IF (ind_Gplus==0) Dens_Gplus = Zero
    Dens_Gminus = v_variab(bv_specy+ind_Gminus)
    IF (ind_Gminus==0) Dens_Gminus = Zero
    Dens_e = v_variab(bv_specy+ind_e)
    IF (ind_e==0) Dens_e = Zero

    !--- 'proton' density (cm-3) nH=n(H)+2n(H2)+n(H+) ---
    nH = Dens_H + 2._DP * Dens_H2 + Dens_Hplus

    !-------------------------------------------------------------------
    ! chemical species (used in CHEMISTRY)
    ! density, velocity, temperature
    !-------------------------------------------------------------------
    !--- ordinary species ---
    speci(1:Nspec)%density = v_variab(bv_speci:ev_speci)
    ! neutrals : V = Vn, T = Tn
    speci(b_neu:e_neu)%velocity    = Vn
    speci(b_neu:e_neu)%temperature = Tn
    ! species on grain mantles : V = Vi, T = Tgrain (not used)
    speci(b_gra:e_gra)%velocity    = Vi
    speci(b_gra:e_gra)%temperature = Tgrain
    ! species in grain cores : V = Vi, T = Tgrain (not used)
    speci(b_cor:e_cor)%velocity    = Vi
    speci(b_cor:e_cor)%temperature = Tgrain
    ! positive ions : V = Vi, T = Ti
    speci(b_ion:e_ion)%velocity    = Vi
    speci(b_ion:e_ion)%temperature = Ti
    ! negative ions : V = Vi, T = Te
    speci(b_neg:e_neg)%velocity    = Vi
    speci(b_neg:e_neg)%temperature = Te
    ! Species with their own velocities:
    speci(ind_CH)%velocity = v_CH
    speci(ind_S )%velocity = v_S
    speci(ind_SH)%velocity = v_SH
    !--- Neutral grains ---
    speci(ind_G)%velocity    = Vi
    speci(ind_G)%temperature = Ti

    !-------------------------------------------------------------------
    ! compression factors
    !-------------------------------------------------------------------
    compr_n = v_variab(iv_compr_n)
    compr_i = v_variab(iv_compr_i)

    !-------------------------------------------------------------------
    ! dust mass, adsorption, erosion, and effective radius
    ! Note : A_grc, dens_grc (which should be constant)
    !        can increase because of compression
    !-------------------------------------------------------------------
    A_grc    = A_grc0 * compr_i
    dens_grc = dens_grc0 * compr_i
    Mgrc     = DOT_PRODUCT(v_variab(bv_cor:ev_cor),speci(b_cor:e_cor)%mass)
    Mgrm     = DOT_PRODUCT(v_variab(bv_gra:ev_gra),speci(b_gra:e_gra)%mass)
    Mgrain   = Mgrc + Mgrm
    ab_cor   = SUM( v_variab(bv_cor:ev_cor) )
    ab_ads   = SUM( v_variab(bv_gra:ev_gra) )
    IF ( (shock_type(1:1) == "S") .OR. &
         (shock_type(1:1) == "P") .OR. &
         (shock_type(1:1) == "W") ) THEN
       d_ero = 0.0_DP
    ELSE
       CALL COMPUTE_DERO
    ENDIF
    CALL COMPUTE_DADS
    ! d_ero   = v_variab(iv_d_ero)
    ! d_ads   = v_variab(iv_d_ads)
    d_site  = ( Mgrm / ab_ads / rho_grm )**(1.0_DP/3.0_DP)
    Nlayers = d_ads / d_site
    rsq_grc = ( f3_mrn - 2.0_DP * f2_mrn * d_ero + f1_mrn * d_ero**2.0_DP ) / f1_mrn
    ! 1st formalism
    rsq_grm = ( f3_mrn + 2.0_DP * f2_mrn * d_ads + f1_mrn * d_ads**2.0_DP ) / f1_mrn
    !!!!! WRITE(*,*) d_ero, d_ads,  ab_cor, ab_ads, Mgrm
    !!!!! DO i = bv_gra, ev_gra
    !!!!!    WRITE(*,*) speci(i-bv_speci+1)%name, v_variab(i), v_dvariab(i)
    !!!!! ENDDO
    r_grm   = sqrt(rsq_grm)
    ! 2nd formalism
    ! rsq_grm  = ( f3_mrn &
    !            + 2.0_DP * f2_mrn * ( d_ads - d_ero ) &
    !            + f1_mrn * ( d_ero**2.0_DP + d_ads**2.0_DP - 2.0_DP * d_ero * d_ads ) ) / f1_mrn
    ! r_grm    = sqrt(rsq_grm)

    !-------------------------------------------------------------------
    ! Add grains in charges density and mass for the computation of 
    ! Vmagnet and Vi. This is not done in RhoI, RhoNEG, DensityI or 
    ! DensityNEG because momentum and energy exchanges are based on 
    ! ions-neutral collision cross sections - BG2017
    !-------------------------------------------------------------------
    DensityCharges = DensityI &
                   + DensityNEG &
                   + v_variab(bv_specy+ind_G)
    RhoCharges     = RhoI + RhoNEG + Mgrain
    muCharges      = RhoCharges / DensityCharges

    ! recompute magnetosonic speed (cm.s-1) allowing for charged grains
    ! Attention - neutral grains also contribute to inertia
    Vmagnet = Gamma * kB * DensityCharges * (Ti+Te) / RhoCharges &
            + ( Bfield * Vs_cm / Vi )**2._DP / ( 4._DP * pi * RhoCharges )
    Vmagnet = SQRT(Vmagnet)

    !-------------------------------------------------------------------
    ! density of each H2 level
    !-------------------------------------------------------------------
    IF (NH2_lev_var>1) THEN
       H2_lev(1:NH2_lev)%density = v_variab(bv_H2_lev:ev_H2_lev)
    ELSE
       H2_lev(1)%density = Dens_H2 * op_H2_in / (1.d0 + op_H2_in)
       H2_lev(2)%density = Dens_H2 * 1.d0 / (1.d0 + op_H2_in)
    ENDIF

    ! JLB + GF (25 VI 2002) - Update H2_lev%Form_gr if required
    IF (iforH2 == 4) THEN
      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%density / Dens_H2
      H2_int_E = SUM(H2_lev(1:NH2_lev)%Form_gr * H2_lev(1:NH2_lev)%Energy)
    ENDIF

    !-------------------------------------------------------------------
    ! In subroutine COMPUTE_OP are calculated the 
    ! ortho:para H2 ratio, and the densities of ortho-H2 
    ! and para-H2 (cm-3).
    ! results : op_H2, Dens_orthoH2, Dens_paraH2
    ! (used in FINE_STRUCTURE_LINE, called in SOURCE)
    !-------------------------------------------------------------------
    CALL COMPUTE_OP_H2

    !-------------------------------------------------------------------
    ! Treatment of the radiation field depending quantities
    ! 1 - compute the dust absorption coefficients (change with r_grm)
    !     compute the rovibrational level populations of CO
    ! 2 - compute the radiative transfer if Av is integrated and
    !     therefore subsequently update the
    !     - FGK coefficients
    !     - phodestruction and photoheating rates
    ! 3 - necessarily update the
    !     - UV pumping of H2 and its photodestruction rate
    !     - photoejection rate from grains
    !     - grain temperature
    !-------------------------------------------------------------------

    ! 1 --- !
    ! Comment the following lines to have constant dust absorption coef
    CALL COMPUTE_QCOEF
    CALL INTERP_QABSO(nwlg,qabso_D_rf,wlg)

    CALL BOLTZMANN_ROVIB_CO

    ! 2 --- !
    IF ( ( F_AV /= 0 ) .AND. ( F_COUP_RAD == 1 ) .AND. ( Z /= Tout_V ) ) THEN
       !--------------------------------------------------
       ! set radiation field intensity to previous value
       ! compute the increment in optical depth
       ! and solve the radiative transfer
       !--------------------------------------------------
       DO i = 1, nangle
          irf(i,1:nwlg) = irf_old(i,1:nwlg)
       ENDDO
       doptdpth(1:nwlg) = pi * rsq_grm * qabso_D_rf(1:nwlg) * dens_grc
       DO i = 1, Ncross
          doptdpth(1:nwlg) = doptdpth(1:nwlg) + sigma_G_rf(i,1:nwlg) * speci(phdest(i)%ispe)%density
       ENDDO
       doptdpth(1:nwlg) = doptdpth(1:nwlg) * diff_dist
       ! optdpth(1:nwlg)  = optdpth_old(1:nwlg) + doptdpth(1:nwlg)
       CALL SIMPLE_TRANSFER

       !--------------------------------------------------
       ! compute Doppler width weighted H2 and CO column
       ! densities and update fgk coefficients
       !--------------------------------------------------
       bt_s_bd1 = vturb / SQRT(vturb*vturb + 2.0_dp * kB * Tn_old / speci(ind_H2)%mass)
       bt_s_bd2 = vturb / SQRT(vturb*vturb + 2.0_dp * kB * Tn     / speci(ind_H2)%mass)
       DO i = 1, NH2_lev
          H2_lev(i)%cd_l = H2_lev(i)%cd_l_old + 0.5_DP * diff_dist &
                         * (bt_s_bd1 * H2_lev(i)%Dens_old + bt_s_bd2 * H2_lev(i)%density)
       ENDDO
       bt_s_bd1 = vturb / SQRT(vturb*vturb + 2.0_dp * kB * Tn_old / speci(ind_CO)%mass)
       bt_s_bd2 = vturb / SQRT(vturb*vturb + 2.0_dp * kB * Tn     / speci(ind_CO)%mass)
       DO i = 1, NCO_lev
          CO_lev(i)%cd_l = CO_lev(i)%cd_l_old + 0.5_DP  * diff_dist &
                         * (bt_s_bd1 * CO_lev(i)%Dens_old + bt_s_bd2 * CO_lev(i)%density)
       ENDDO
       CALL FGKCOEF

       !--------------------------------------------------
       ! compute photodestruction and heating rates
       ! => used in CHEMISTRY and SOURCE
       !--------------------------------------------------
       DO i = 1, Ncross
          phdest_rate(i)    = SECT_INT(phdest(i))
          phheating_rate(i) = HEATING_PHOTOCHEM_SECT_INT(phdest(i))
       ENDDO
    ENDIF

    ! 3 --- !
    !-----------------------------------------------------
    ! compute H2 pumping rates and the dissociation rates 
    ! of H2 and CO
    ! => used in CHEMISTRY
    ! Note : the dissociation rates of H2 and CO are 
    !        computed within the LVG approximation. While
    !        this is correct for H2 and has been checked
    !        on numerous occasion, this is not appropriate
    !        for CO where the natural linewidths are too
    !        large to use FGK. Those rates are therefore
    !        not used in the chemistry where we apply the
    !        table of Lee et al. (1996) to compute CO
    !        dissociation
    !-----------------------------------------------------
    CALL FGKDATA

    !-----------------------------------------------------
    ! compute photoejection rate from grains
    ! => used in CHEMISTRY and SOURCE
    !-----------------------------------------------------
    CALL COMPUTE_PHOTOELEC

    !-----------------------------------------------------
    ! compute grain temperature
    ! => used in CHEMISTRY and SOURCE
    !-----------------------------------------------------
    IF( F_TGR == 1 ) THEN
       CALL COMPUTE_GRAIN_TEMP
       speci(b_gra:e_gra)%temperature = Tgrain
       speci(b_cor:e_cor)%temperature = Tgrain
    ENDIF


    ! ===========================================================================
    ! STEP 2 - COMPUTE THE SOURCES TERMS (CHEMISTRY AND DYNAMICS)
    ! ===========================================================================

    !-------------------------------------------------------------------
    ! In subroutine CHEMISTRY the source terms of ([1], equations(3)-(14))
    ! are calculated as functions of the chemistry and the abundances of
    ! the species and the global properties of the fluids.              
    ! results in :                                                      
    !    YN                                                             
    !    CupN,   CupI,   CupNEG                                         
    !    CdownN, CdownI, CdownNEG                                       
    !    YNn,    YNi,    YNneg                                          
    !    Sn,     Si,     Sneg                                           
    !    An,     Ai,     Aneg                                           
    !    Bn,     Bi,     Bneg                                           
    !-------------------------------------------------------------------
    CALL CHEMISTRY

    !-------------------------------------------------------------------
    ! In subroutine SOURCE, the source terms for MHD equations are      
    ! computed depending on the results of subroutine CHEMISTRY and     
    ! the function values v_variab(i)                                   
    ! SOURCE also computes H2 cooling and internal energy.              
    ! results in :                                                      
    !   An, Ai, Aneg                                                    
    !   Bn, Bi, Bneg                                                    
    !   H2_energy, YN_rovib_H2                                          
    !-------------------------------------------------------------------
    CALL SOURCE

    !-------------------------------------------------------------------
    ! add chemical contribution YN (from CHEMISTRY) to the source of    
    ! internal energy of H2 H2_energy (from COMPUTE_H2 via SOURCE)      
    !-------------------------------------------------------------------

    !--- Sel_tot_H2 is the rate of destruction of H2 by collisional dissociation
    Sel_tot_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_ch_H2))
    Sel_tne_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_ne_H2))
    Sel_tio_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_io_H2))
    Sel_tel_H2 = SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_el_H2))

    !-------------------------------------------------------------------
    ! H2_energy was used in the old version in the evolution equations  
    ! of Tn, Vn, gradVn. This term appears in the equations if we put   
    ! the internal energy of H2 (·niEi) in the energy conservation      
    ! equation (see BG notes - deriving shock equations).               
    ! However this term is useless if we compute the energy conservation
    ! based on the kinetic energy and the thermal energy only. Careful, 
    ! in this case, we need to modify the way the H2 cooling is         
    ! computed in Bn                                                    
    !-------------------------------------------------------------------
    if (NH2_lev_var>1) then
       H2_int_energy = SUM(DBLE(H2_lev(1:NH2_lev)%density * H2_lev(1:NH2_lev)%energy))
       H2_energy = H2_energy + kB * (YN(ind_H2) - Sel_tot_H2) / sum_H2 * &
            H2_int_energy &
            + kB * SUM(DBLE(H2_lev(1:NH2_lev)%density * Sel_ch_H2 * H2_lev(1:NH2_lev)%energy))
    else
       H2_int_energy=Zero
       H2_energy=Zero
    endif

    !-------------------------------------------------------------------
    ! Concatenate source terms depending on Nfluids                     
    !-------------------------------------------------------------------
    !
    !-------------------------------------------------------------------
    ! 29 mai 2001 - JLB : molecular cooling not included in Bn any more 
    !                     but computed directly here. Old expression was
    !                     with molec_cool in Bn, via cooling_n          
    ! one fluid model : Bne = Bn + Bi + Bneg                            
    !-------------------------------------------------------------------

    SELECT CASE (Nfluids)
    CASE (1)
       ! one-fluid model : Bne = Bn + Bi + Bneg; Ane = 0; Sne = 0
       Bne = Bn + Bi + Bneg
       Ane = 0.0_DP
       Sne = 0.0_DP
       YNne = YNn + YNi + YNneg
       DensityNe = DensityN + DensityI + DensityNeg
       RhoNe = RhoN + RhoI + RhoNEG
    CASE (2)
       ! two-fluids model : Bne = Bn
       Bne = Bn
       Ane = An
       Sne = Sn
       YNne = YNn
       DensityNe = DensityN
       RhoNe = RhoN
    CASE (3)
       ! three-fluids model : Bne = Bn
       Bne = Bn
       Ane = An
       Sne = Sn
       YNne = YNn
       DensityNe = DensityN
       RhoNe = RhoN
    END SELECT


    ! ===========================================================================
    ! STEP 3 - CALCULATION OF THE DERIVATIVES OF ALL VARIABLES WITH RESPECT
    !          TO Z (DISTANCE INTO THE CLOUD)
    !          FOR WIND MODELS, A SPHERICAL GEOMETRY IS ASSUMED
    !          Z+Z0=R (DISTANCE TO THE ORIGIN OF THE MASS FLOW)
    ! ===========================================================================

    !-------------------------------------------------------------------
    ! A - velocities Vn, Vi, Ve                                         
    !   - in the case of static models (S1, S2, P1, P2), Vn is          
    !     considered as constant.                                       
    !   - in the case of Wind models, Vn is either considered as        
    !     constant (W1) or is computed self consistently (W2) combining 
    !     the mass, momentum, and energy conservation equations (see    
    !     notes BG 2019)                                                
    !   - It should be noted that a constant velocity wind necessarily  
    !     violate the momentum conservation equation, as an unknown     
    !     external force is necessary applied to sustain the gas        
    !     velocity. The same remarks holds of course for static         
    !     isochoric or isobaric models with constant velocity           
    !   - for J- or C-type shocks, the molecular cooling is a function  
    !     of DY(iv_Vn), which leads to an implicit equation             
    !     But in case shock_type='C', grad_V is not a Vode variable     
    !     (DY(iv_gv)=0), so we use a Newton method to implicit this     
    !     equation at the level of diffun.                              
    !     (Ideally, we should use an Algebraic Differential Equations   
    !     solver rather than Vode)                                      
    !     The same problem also holds for self-consistent spherical     
    !     wind models where grad_V is not a Vode variable. For this     
    !     case we should apply the same implicit method than that       
    !     derived for C-type shocks or else derive an equation for      
    !     grad_V. This is not done yet and the molecular cooling is     
    !     computed in this case assuming grad_V = 0                     
    !-------------------------------------------------------------------
    !
    !-------------------------------------------------------------------
    ! Note : new integrations of the evolution equations (BG 2017)      
    !        We dont include the internal energy of H2 in the energy    
    !        conservation equation anymore. H2 cooling is treated as    
    !        atomic and chemical cooling in the source terms directly   
    !-------------------------------------------------------------------

    !-----------------------------------------------------
    ! Compute Vgrad and molecular cooling
    ! SC Mar 2004 : Vgrad should be larger than 1 km s-1 pc-1
    !               minimum gradient = thermal gradient
    !               over scale z=(distance+dmin)
    !-----------------------------------------------------
    dmin = 1.0D13
    Vgrad = SQRT(grad_V*grad_V + vsound*vsound/(distance+dmin)**2) * 1.0D-5 ! SC
    CALL COMPUTE_MOLECULAR_COOLING

    !-----------------------------------------------------
    ! A1 - Neutral velocity
    !-----------------------------------------------------

    SELECT CASE (shock_type(1:1))

    CASE ('P','S')
       DY(iv_Vn) = 0_DP

    CASE ('W')
       SELECT CASE (shock_type(2:2))
       !--------------------------------------------------
       ! type 1 - constant velocity spherical wind
       !--------------------------------------------------
       CASE ('1')
          DY(iv_Vn) = 0_DP
       !--------------------------------------------------
       ! type 2 - self-consistent spherical wind
       !        - be careful, the wind might cross a sonic
       !          point
       !--------------------------------------------------
       CASE ('2')
          aux1_dvn = Bne - 2.0_DP * GAMMA2 * Vn * DensityNe * kB * Tn / (Z+Z0)
          aux2_dvn = (DensityNe*GAMMA2*kB*Tn-RhoNe*GAMMA1*Vn2)
          DY(iv_Vn) = (aux1_dvn - molec_cool) / aux2_dvn
       END SELECT

    CASE ('J')
       aux1_dvn = GAMMA3*Sne*Vn2 - GAMMA2*Ane*Vn + Bne
       aux2_dvn = (DensityNe*GAMMA2*kB*Tn-RhoNe*GAMMA1*Vn2)
       !--------------------------------------------------
       ! Add Magnetic term for single fluid J shocks
       !--------------------------------------------------
       IF (  Nfluids == 1 ) THEN
         aux2_dvn = aux2_dvn + GAMMA1 * Bfield2
       ENDIF
       !--------------------------------------------------
       ! 1 - Compute Vn derivative without viscosity
       ! 2 - Save this value for comparison with the 2nd 
       !     order term (including viscosity) in mhd_vode
       ! 3 - If viscosity, include second order term
       !     If not, use the derivative without viscosity
       !--------------------------------------------------
       dvn2    = (aux1_dvn - molec_cool) / aux2_dvn
       save_dv = dvn2
       IF ( viscosity ) THEN
          DY(iv_Vn) = - grad_V
       ELSE
          DY(iv_Vn) = dvn2
       ENDIF

    CASE ('C')
       aux1_dvn = GAMMA3*Sne*Vn2 - GAMMA2*Ane*Vn + Bne
       aux2_dvn = (DensityNe*GAMMA2*kB*Tn-RhoNe*GAMMA1*Vn2)
       !--------------------------------------------------
       ! Trick to implicit grad_V (SC May 06)
       !      dvn1 = (aux1_dvn - molec_cool) / aux2_dvn
       !      Vgrad = SQRT(dvn1*dvn1 + 1.0D10 / (parsec*parsec)) * 1.0D-5
       !--------------------------------------------------
       dvn1 = v_l_der(iv_Vn)  
       Vgrad = SQRT(dvn1*dvn1 + vsound*vsound/(distance+dmin)**2) * 1.0D-5
 
       CALL COMPUTE_MOLECULAR_COOLING
       dvn2 = (aux1_dvn - molec_cool) / aux2_dvn
       gvn1 = dvn1 - dvn2
 
       Vgrad = SQRT(dvn2*dvn2 + vsound*vsound/(distance+dmin)**2) * 1.0D-5
       CALL COMPUTE_MOLECULAR_COOLING
       gvn2 = dvn2 - (aux1_dvn - molec_cool) / aux2_dvn
 
       !--------------------------------------------------
       ! SC: proposed loop for Newton's method of solving 
       !     for dVn/dz
       !--------------------------------------------------
       ii = 0
       DO WHILE ( abs(gvn1-gvn2) > 1.0e-6*(abs(gvn1)).AND. ii < 1000)
         corrdv = dvn2 - gvn2 * (dvn1 - dvn2) / (gvn1 - gvn2)
         dvn1 = dvn2
         gvn1 = gvn2
         dvn2 = corrdv
         Vgrad = SQRT(dvn2*dvn2 + vsound*vsound/(distance+dmin)**2) * 1.0D-5
         CALL COMPUTE_MOLECULAR_COOLING
         gvn2 = dvn2 - (aux1_dvn - molec_cool) / aux2_dvn
         ii = ii + 1
       ENDDO
 
       DY(iv_Vn) = dvn2
       Vgrad = SQRT(dvn2*dvn2 + vsound*vsound/(distance+dmin)**2) * 1.0D-5

    END SELECT

    !-----------------------------------------------------
    ! A2 - Ions / electrons velocity
    !      according to the number of fluids
    !-----------------------------------------------------

    SELECT CASE (Nfluids)
    !-----------------------------------------------------
    ! one-fluid model : Vi = Vn
    !-----------------------------------------------------
    CASE (1)
       DY(iv_Vi) = DY(iv_Vn)
    !-----------------------------------------------------
    ! two or three fluids model :
    ! ions and electrons have the same velocity
    !-----------------------------------------------------
    CASE DEFAULT
       DY(iv_Vi) = (-GAMMA3*Sn*Vi2 + GAMMA2*An*Vi + (Bi + Bneg)) &
                 / ( GAMMA2 * kB * DensityI * (Ti + Te) + GAMMA1*Bfield2 &
                   - GAMMA1 * Vi2 * RhoCharges )
    END SELECT


    !-------------------------------------------------------------------
    ! B - Compression factors                                           
    !   - in the case of static models, there is either no compression  
    !     (S1, P1), or an isobaric compression (S2, P2)                 
    !   - in the case of wind models, the compression is either computed
    !     in the framework of a constant velocity wind (W1) or computed 
    !     self consistently (W2), taking into account the variation of  
    !     the velocity along the Z (i.e. radial) direction              
    !   - for J- or C-type shocks, compression is computed              
    !     self-consistently                                             
    !   Note : for all the cases, we compute effective velocity         
    !          divergence terms, divu_n and divu_i, which are used to   
    !          take into account compression factors even if the        
    !          velocity along Z is considered as constant               
    !          very useful for isobaric models (to simply account for   
    !          the work of external forces) and wind models (to simply  
    !          account for spherical dilution)                          
    !-------------------------------------------------------------------

    SELECT CASE (shock_type(1:1))

    CASE ('P','S')
       SELECT CASE (shock_type(2:2))
       !--------------------------------------------------
       ! type 1 - isochoric model
       !--------------------------------------------------
       CASE ('1')
          divu_n = 0.0_DP
          divu_i = divu_n
       !--------------------------------------------------
       ! type 2 - isobaric model
       !--------------------------------------------------
       CASE ('2')
          divu_n = ( Bne - molec_cool ) / ( GAMMA2 * DensityNe * kB * Tn )
          divu_i = divu_n
       END SELECT

    CASE ('W')
       SELECT CASE (shock_type(2:2))
       !--------------------------------------------------
       ! type 1 - constant velocity spherical wind
       !        - geometric dilution only
       !--------------------------------------------------
       CASE ('1')
          divu_n = 2.0_DP / (Z+Z0) * Vn
          divu_i = divu_n
       !--------------------------------------------------
       ! type 2 - self-consistent spherical wind
       !--------------------------------------------------
       CASE ('2')
          divu_n = 2.0_DP / (Z+Z0) * Vn + DY(iv_Vn)
          divu_i = divu_n
       END SELECT

    CASE ('J', 'C')
       divu_n = DY(iv_Vn)
       divu_i = DY(iv_Vi)

    END SELECT

    DY(iv_compr_n) = - v_variab(iv_compr_n) * divu_n / Vn
    DY(iv_compr_i) = - v_variab(iv_compr_i) * divu_i / Vi


    !-------------------------------------------------------------------
    ! C - densities (mass and volume)
    !     Once the effective velocity divergence divu of a fluid is
    !     computed, all densities simply write
    !     Dt(n)   = Y - n   * divu   (Y = volume density source term)
    !     Dt(rho) = S - rho * divu   (S = mass   density source term)
    !     time derivative translate into spatial derivative via
    !     Dt = Uz * Dz
    !-------------------------------------------------------------------

    !-----------------------------------------------------
    ! mass densities of the fluids
    !-----------------------------------------------------
    DY(iv_RhoN)       = ( Sn    - RhoN       * divu_n ) / Vn
    DY(iv_RhoI)       = ( Si    - RhoI       * divu_i ) / Vi
    DY(iv_RhoA)       = ( Sa    - RhoA       * divu_i ) / Vi
    DY(iv_RhoNEG)     = ( Sneg  - RhoNEG     * divu_i ) / Vi

    !-----------------------------------------------------
    ! volume densities of the fluids
    !-----------------------------------------------------
    DY(iv_DensityN)   = ( YNn   - DensityN   * divu_n ) / Vn
    DY(iv_DensityI)   = ( YNi   - DensityI   * divu_i ) / Vi
    DY(iv_DensityA)   = ( YNa   - DensityA   * divu_i ) / Vi
    DY(iv_DensityNeg) = ( YNNeg - DensityNeg * divu_i ) / Vi

    !-----------------------------------------------------
    ! volume densities of the chemical species
    !-----------------------------------------------------
    DY(bv_neu:ev_neu) = ( YN(b_neu:e_neu) - v_variab(bv_neu:ev_neu) * divu_n ) / Vn
    DY(bv_gra:ev_gra) = ( YN(b_gra:e_gra) - v_variab(bv_gra:ev_gra) * divu_i ) / Vi
    DY(bv_cor:ev_cor) = ( YN(b_cor:e_cor) - v_variab(bv_cor:ev_cor) * divu_i ) / Vi
    DY(bv_ion:ev_ion) = ( YN(b_ion:e_ion) - v_variab(bv_ion:ev_ion) * divu_i ) / Vi
    DY(bv_neg:ev_neg) = ( YN(b_neg:e_neg) - v_variab(bv_neg:ev_neg) * divu_i ) / Vi

    !-----------------------------------------------------
    ! volume densities of the ro-vibrational levels of H2 
    !-----------------------------------------------------
    ! Note : The contribution of chemical reactions       
    !        YN(ind_H2) is included in proportion to the  
    !        fractional population of the level, namely   
    !        v_variab(J) / sum_H2, except for a few       
    !        reactions where specific ponderations are    
    !        known                                        
    !-----------------------------------------------------
    IF ( NH2_lev_var > 1) THEN
       DY(bv_H2_lev:ev_H2_lev) = &
            ( YN_rovib_H2(1:NH2_lev) &
            + (YN(ind_H2) - Sel_tot_H2 - For_gr_H2) * v_variab(bv_H2_lev:ev_H2_lev) / sum_H2 &
            + H2_lev(1:NH2_lev)%density * Sel_ch_H2 &
            + For_gr_H2 * H2_lev(1:NH2_lev)%Form_gr &
            - v_variab(bv_H2_lev:ev_H2_lev) * divu_n ) / Vn
    ENDIF


    !-------------------------------------------------------------------
    ! D - temperatures                                                  
    !   - in the case of isobaric models (S2 and P2), the temperature   
    !     evolution equation is deduced from that of the density in     
    !     order to keep the pressure constant                           
    !   - in all other cases, the temperature is computed from the      
    !     internal energy evolution equation taking into account the    
    !     production of momentum and kinetic energy within each fluid   
    !     if a multi-fluid evolution is computed                        
    !   - viscous heating is including for J-type shocks                
    !-------------------------------------------------------------------

    !-----------------------------------------------------
    ! D1 - Neutral temperature
    !-----------------------------------------------------
    SELECT CASE (shock_type(1:1))

    CASE ('P','S')
       SELECT CASE (shock_type(2:2))
       !--------------------------------------------------
       ! type 1 - isochoric model
       !--------------------------------------------------
       CASE ('1')
          DY(iv_Tn) = Bne - molec_cool &
                    - GAMMA1 * kB * Tn * YNne
          DY(iv_Tn) = DY(iv_Tn) / (GAMMA1 * DensityNe * kB * Vn )
       !--------------------------------------------------
       ! type 2 - isobaric model
       !        - Dz(T) deduced from Dz(N) using P=cst
       !--------------------------------------------------
       CASE ('2')
          DY(iv_Tn) = - Tn / DensityNe * YNne + Tn * divu_n
          DY(iv_Tn) = DY(iv_Tn) / Vn
       END SELECT

    CASE ('W')
       !--------------------------------------------------
       ! Same equation for type 1 or type 2 winds
       ! Temperature evolution deduced from internal
       ! energy equation:
       !    => Dt(e) = - P / rho * divu + (G-L) / rho
       !--------------------------------------------------
       DY(iv_Tn) = Bne - molec_cool &
                 - GAMMA1 * kB * Tn * YNne &
                 - DensityNe * kB * Tn * divu_n
       DY(iv_Tn) = DY(iv_Tn) / (GAMMA1 * DensityNe * kB * Vn )

    CASE ('J', 'C')
       !--------------------------------------------------
       ! Same equation for J-type shocks or C-type shocks
       ! (modulo viscosity effects only taken into account
       ! in J-type shocks)
       !--------------------------------------------------
       DY(iv_Tn) = Bne - molec_cool - Ane * Vn + 0.5_dp * Sne * Vn2 &
                 - GAMMA1 * kB * Tn * YNne &
                 - DensityNe * kB * Tn * divu_n
       IF (viscosity) THEN
          DY(iv_Tn) = DY(iv_Tn) &
                    - RhoNe * XLL2 * DY(iv_Vn) * DY(iv_Vn) * DY(iv_Vn)
       ENDIF
       DY(iv_Tn) = DY(iv_Tn) / (GAMMA1 * DensityNe * kB * Vn )

    END SELECT

    !-----------------------------------------------------
    ! D2 - Ion and electron temperatures
    !      according to the number of fluids
    !-----------------------------------------------------

    SELECT CASE (Nfluids)
    !-----------------------------------------------------
    ! one-fluid model : Ti = Te = Tn
    !-----------------------------------------------------
    CASE (1)
       DY(iv_Ti) = DY(iv_Tn)
       DY(iv_Te) = DY(iv_Tn)
    !-----------------------------------------------------
    ! two-fluid model : Ti = Te
    ! note : Ai + Aneg = - An
    !        Si + Sneg = - sn
    !        ... hence the expressions given below
    !-----------------------------------------------------
    CASE (2)
       DY(iv_Ti) = 0.5_DP * (Bi + Bneg) + 0.5_DP * An * Vi -0.25_DP * Sn * Vi2 &
                 - GAMMA1 * kB * Ti * YNi &
                 - DensityI * kB * Ti * divu_i
       DY(iv_Ti) = DY(iv_Ti) / ( GAMMA1 * DensityI * kB * Vi )
       DY(iv_Te) = DY(iv_Ti)
    !-----------------------------------------------------
    ! three-fluid model
    !-----------------------------------------------------
    CASE (3)
       DY(iv_Ti) = Bi - Ai * Vi + 0.5_DP * Si * Vi2 &
                 - GAMMA1 * kB * Ti * YNi &
                 - DensityI * kB * Ti * divu_i
       DY(iv_Ti) = DY(iv_Ti) / ( GAMMA1 * DensityI * kB * Vi )
       DY(iv_Te) = Bneg - Aneg * Vi + 0.5_DP * Sneg * Vi2 &
                 - GAMMA1 * kB * Te * YNi &
                 - DensityI * kB * Te * divu_i
       DY(iv_Te) = DY(iv_Te) / ( GAMMA1 * DensityI * kB * Vi )
    END SELECT

    IF (ieqth == 0) THEN
       DY(iv_Tn) = 0.0
       DY(iv_Ti) = 0.0
       DY(iv_Te) = 0.0
    ENDIF

    !-----------------------------------------------------
    ! save few heating terms for output files
    !-----------------------------------------------------
    mech_trsf_n_chem  = - An   * Vn + 0.5_dp * Sn   * Vn2 - GAMMA1 * kB * Tn * YNn  
    mech_trsf_i_chem  = - Ai   * Vi + 0.5_dp * Si   * Vi2 - GAMMA1 * kB * Ti * YNi  
    mech_trsf_e_chem  = - Aneg * Vi + 0.5_dp * Sneg * Vi2 - GAMMA1 * kB * Te * YNneg
    mech_trsf_n_compr = - DensityN   * kB * Tn * divu_n
    mech_trsf_i_compr = - DensityI   * kB * Ti * divu_i
    mech_trsf_e_compr = - DensityNeg * kB * Te * divu_i
    IF ( viscosity ) THEN
       mech_trsf_n_visc = - RhoNe * XLL2 * DY(iv_Vn) * DY(iv_Vn) * DY(iv_Vn)
    ELSE
       mech_trsf_n_visc = 0.0_DP
    ENDIF
    mech_trsf_i_visc = 0.0_DP
    mech_trsf_e_visc = 0.0_DP
    mech_trsf_n = mech_trsf_n_chem + mech_trsf_n_compr + mech_trsf_n_visc
    mech_trsf_i = mech_trsf_i_chem + mech_trsf_i_compr + mech_trsf_i_visc
    mech_trsf_e = mech_trsf_e_chem + mech_trsf_e_compr + mech_trsf_e_visc
    mech_trsf_tot = mech_trsf_n + mech_trsf_i + mech_trsf_e


    !-------------------------------------------------------------------
    ! E - velocity gradient                                             
    !     only computed for J-type shocks with viscosity                
    !     This equation is otherwise not needed as the derivative of    
    !     the velocity is a first order (and not second order) equation 
    ! Note : the evolution of the velocity gradient is computed here    
    !        because it is expressed as a function of DY(iv_Tn), the    
    !        derivative of the neutral temperature                      
    !-------------------------------------------------------------------

    IF (viscosity) THEN
       IF (Nfluids == 1) THEN
          DY(iv_gv) = Bne - molec_cool - 0.5_dp * Sne * Vn2 - XLL2 * DY(iv_Vn) * DY(iv_Vn) * Sne &
                    - GAMMA2 * kB * Tn * YNne - RhoNe * Vn2 * DY(iv_Vn) &
                    - GAMMA2 * DensityNe * Vn * kB * DY(iv_Tn) + Bfield2 * DY(iv_Vn)
          DY(iv_gv) = DY(iv_gv) / (2.0_dp * RhoNe * XLL2 * DY(iv_Vn) * Vn)
       ELSE
          DY(iv_gv) = Bne - molec_cool - 0.5_dp * Sne * Vn2 - XLL2 * DY(iv_Vn) * DY(iv_Vn) * Sne &
                    - GAMMA2 * kB * Tn * YNne - RhoNe * Vn2 * DY(iv_Vn) &
                    - GAMMA2 * DensityNe * Vn * kB * DY(iv_Tn)
          DY(iv_gv) = DY(iv_gv) / (2.0_dp * RhoNe * XLL2 * DY(iv_Vn) * Vn)
       ENDIF
       !--------------------------------------------------
       ! iv_gv is for -dVn/dz
       !--------------------------------------------------
       DY(iv_gv) = - DY(iv_gv)
    ELSE
       DY(iv_gv) = 0.0_dp
    ENDIF


    !-------------------------------------------------------------------
    ! F - CH, S, and SH velocity AND densities                          
    !   - in the case of static models and wind models (S1, S2, P1, P2, 
    !     W1, W2), those velocities are set to that of the neutrals     
    !   - in the case of J- or C-type shocks, the velocities are either 
    !     set to that of the neutrals or computed as a separate fluid   
    !     according to the corresponding flags (F_CH, F_S, F_SH)        
    !   - in the last case, the evolution equations of CH, S, and SH    
    !     must be modified accordingly                                  
    ! Note : the evolutions of those velocities are computed here       
    !        because they are expressed as functions of DY(iv_Tn), the  
    !        derivative of the neutral temperature                      
    !-------------------------------------------------------------------

    SELECT CASE (shock_type(1:1))

    CASE ('P','S','W')
       DY(iv_vCH) = DY(iv_Vn)
       DY(iv_vS)  = DY(iv_Vn)
       DY(iv_vSH) = DY(iv_Vn)

    CASE DEFAULT
       !CH!
       IF ( F_CH == 1 ) THEN
          DY(iv_vCH) = (YA_CH+A_CH_n)*v_CH-YS_CH*v_CH**2
          DY(iv_vCH) = DY(iv_vCH)-YS_CH*kB*Tn/mass_CH
          DY(iv_vCH) = DY(iv_vCH)-Dens_CH*v_CH*kB*DY(iv_Tn)
          DY(iv_vCH) = DY(iv_vCH)/(mass_CH*v_CH*v_CH-kB*Tn)/Dens_CH
          DY(bv_specy+ind_CH) = ( YN(ind_CH) - v_variab(bv_specy+ind_CH) * DY(iv_vCH) ) / v_CH
       ELSE
          DY(iv_vCH) = DY(iv_Vn)
       ENDIF
       !S!
       IF ( F_S  == 1 ) THEN
          DY(iv_vS)  = (YA_S +A_S_n )*v_S -YS_S *v_S**2
          DY(iv_vS)  = DY(iv_vS )-YS_S *kB*Tn/mass_S
          DY(iv_vS)  = DY(iv_vS )-Dens_S*v_S *kB*DY(iv_Tn)
          DY(iv_vS)  = DY(iv_vS)/(mass_S*v_S*v_S-kB*Tn)/Dens_S
          DY(bv_specy+ind_S ) = ( YN(ind_S ) - v_variab(bv_specy+ind_S ) * DY(iv_vS ) ) / v_S
       ELSE
          DY(iv_vS ) = DY(iv_Vn)
       ENDIF
       !SH!
       IF ( F_SH == 1 ) THEN
          DY(iv_vSH) = (YA_SH+A_SH_n)*v_SH-YS_SH*v_SH**2
          DY(iv_vSH) = DY(iv_vSH)-YS_SH*kB*Tn/mass_SH
          DY(iv_vSH) = DY(iv_vSH)-Dens_SH*v_SH*kB*DY(iv_Tn)
          DY(iv_vSH) = DY(iv_vSH)/(mass_SH*v_SH*v_SH-kB*Tn)/Dens_SH
          DY(bv_specy+ind_SH) = ( YN(ind_SH) - v_variab(bv_specy+ind_SH) * DY(iv_vSH) ) / v_SH
       ELSE
          DY(iv_vSH) = DY(iv_Vn)
       ENDIF

    END SELECT


    !-------------------------------------------------------------------
    ! G - H2 and H column-densities                                     
    !-------------------------------------------------------------------

    IF (F_AV == 0) THEN
       DY(iv_nh ) = 0_DP
       DY(iv_nh2) = 0_DP
       DY(iv_nco) = 0_DP
    ELSE
       IF (inv_Av_fac /= 0d0) THEN
          DY(iv_nh ) = sign(Dens_H, inv_Av_fac)
          DY(iv_nh2) = sign(Dens_H2,inv_Av_fac)
          DY(iv_nco) = sign(Dens_CO,inv_Av_fac)
       ELSE
          DY(iv_nh ) = 0_DP
          DY(iv_nh2) = 0_DP
          DY(iv_nco) = 0_DP
       ENDIF
    ENDIF


    !-------------------------------------------------------------------
    ! H - Grain related quantities                                      
    !     erosion and adsorption widths                                 
    !   - erosion and adsorption widths are now computed analytically   
    !     using a N-R procedure (see above)                             
    !     To switch formalism, uncomment the following lines and use    
    !     the variables iv_d_ero and iv_d_ads                           
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! DY(iv_d_ero)   = ( 3.0 * f3_mrn - 6.0 * f2_mrn * v_variab(iv_d_ero) &
    !                  + 3.0 * f1_mrn * v_variab(iv_d_ero)**2.0 )**(-1.0) &
    !                * f4_mrn / (Mgrc0 * compr_i) &
    !                * ( Mgrc / compr_i * DY(iv_compr_i) &
    !                  - SUM( DY(bv_cor:ev_cor) * speci(b_cor:e_cor)%mass ) )
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! 1st formalism
    !-------------------------------------------------------------------
    ! DY(iv_d_ads)   = ( 3.0 * f3_mrn + 6.0 * f2_mrn * v_variab(iv_d_ads) &
    !                +   3.0 * f1_mrn * v_variab(iv_d_ads)**2.0 )**(-1.0) &
    !                * f4_mrn * rho_grc * d_site**3.0 / (Mgrc0 * compr_i) &
    !                * ( SUM( DY(bv_gra:ev_gra) ) &
    !                  - ab_ads / compr_i * DY(iv_compr_i) )
    !                ! * f4_mrn * rho_grc / rho_grm / (Mgrc0 * compr_i) &
    !                ! * ( SUM( DY(bv_gra:ev_gra) * speci(b_gra:e_gra)%mass ) &
    !                !   - Mgrm / compr_i * DY(iv_compr_i) )
    !-------------------------------------------------------------------


    !-------------------------------------------------------------------
    ! save values of v_variab and their derivative DY
    ! to the arrays v_l_var and v_l_der.
    !-------------------------------------------------------------------
    v_l_var(1:N)=v_variab(1:N)
    v_l_der(1:N)=DY(1:N)

    !-------------------------------------------------------------------
    ! Write variables & derivatives on screen - for debug purposes
    !-------------------------------------------------------------------
    ! DO i = 1, Nv_MHD
    !    WRITE(*,*) v_variab(i), DY(i)
    ! ENDDO
    ! DO i = bv_speci, ev_speci
    !    WRITE(*,*) speci(i-bv_speci+1)%name, v_variab(i), DY(i)
    ! ENDDO
    ! DO i = bv_H2_lev, ev_H2_lev
    !    WRITE(*,*) " lev ", i-bv_H2_lev, v_variab(i), DY(i)
    ! ENDDO
    !-------------------------------------------------------------------


    ! ===========================================================================
    ! STEP 4 - SWITCH TO LIN-LIN, LOG-LIN OR LOG-LOG DERIVATIVES
    ! ===========================================================================

    !-------------------------------------------------------------------
    ! switch to logarithmic derivative : DY = DY / Y
    ! - choose if we perform linear or logarithmic integration of z
    !   DY = DY   or   DY = DY * EXP(Z)
    !-------------------------------------------------------------------
    IF      (integ_type == 0) THEN
       DY(1:Nv_MHD) = DY(1:Nv_MHD) / (v_variab(1:Nv_MHD) + tiny)
       DY(bv_speci:ev_speci) = DY(bv_speci:ev_speci) / (v_variab(bv_speci:ev_speci) + tiny)
       IF (NH2_lev_var>1) THEN
          DY(bv_H2_lev:ev_H2_lev) = DY(bv_H2_lev:ev_H2_lev) / (v_variab(bv_H2_lev:ev_H2_lev) + tiny)
       ENDIF
    ELSE IF (integ_type == 1) THEN
       DY(1:Nv_MHD) = DY(1:Nv_MHD) / (v_variab(1:Nv_MHD) + tiny) * EXP(Z)
       DY(bv_speci:ev_speci) = DY(bv_speci:ev_speci) / (v_variab(bv_speci:ev_speci)+tiny) * EXP(Z)
       IF (NH2_lev_var>1) THEN
          DY(bv_H2_lev:ev_H2_lev) = DY(bv_H2_lev:ev_H2_lev) / (v_variab(bv_H2_lev:ev_H2_lev) + tiny) * EXP(Z)
       ENDIF
    ENDIF

    !-------------------------------------------------------------------
    ! Add molec_cool to Bn as before, so it may be used without         
    ! modifications in other parts of the code                          
    ! Note : New integrations of the evolution equations. We dont       
    !        include internal energy of H2 in the energy conservation   
    !        but treat the cooling by H2 as a source term (like any     
    !        other atomic or molecular radiative process)               
    !-------------------------------------------------------------------
    Bn = Bn - molec_cool
    IF ( (shock_type == "J") .AND. (viscosity) ) THEN
      Bn = Bn - RhoNe * XLL2 * DY(iv_Vn) * DY(iv_Vn) * DY(iv_Vn)
    ENDIF
    total_n_cool = total_n_cool + molec_cool

  END SUBROUTINE DIFFUN


END MODULE MODULE_EVOLUTION
