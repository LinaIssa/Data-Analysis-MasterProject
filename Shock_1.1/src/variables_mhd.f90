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

MODULE MODULE_GAMMA

  !*****************************************************************************
  !** The module 'MODULE_GAMMA' contains gamma (ratio of specific heats)      **
  !** and useful related variables.                                           **
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! gamma=5/3 for 3 degrees of freedom
  REAL(KIND=DP), PRIVATE, PARAMETER :: freedom_deg=3._DP ! degrees of freedom
  REAL(KIND=DP), PARAMETER :: GAMMA=(freedom_deg+2._DP)/freedom_deg
  REAL(KIND=DP), PARAMETER :: GAMMA1=1._DP/(GAMMA-1._DP)
  REAL(KIND=DP), PARAMETER :: GAMMA2=GAMMA/(GAMMA-1._DP)
  REAL(KIND=DP), PARAMETER :: GAMMA3=0.5_DP*(GAMMA+1._DP)/(GAMMA-1._DP)

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_GAMMA


MODULE MODULE_DEBUG_JLB

  !*****************************************************************************
  !** The module 'DEBUG_JLB' contains variables used to check                 **
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  REAL(KIND=DP) :: tr_1, tr_2, tr_3, tr_4, tr_5, tr_6
  REAL(KIND=DP) :: tr_7, tr_8, tr_9, tr_10, tr_11, tr_12
  REAL(KIND=DP) :: tr_13, tr_14, tr_15, tr_16, tr_17, tr_18, tr_19
  REAL(KIND=DP) :: tr_20, tr_21, tr_22, tr_23, tr_24, tr_25, tr_26
  REAL(KIND=DP) :: tr_27, tr_28, tr_29, tr_30

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_DEBUG_JLB


MODULE MODULE_GRAINS

  !*****************************************************************************
  !** The module 'GRAINS' contains variables related to the dust grains       **
  !*****************************************************************************
  ! SC May 2006 changes made since last release 16II04
  ! - lay: added variable R_gr_scale12 = <R> / <R**2>**1/2 
  !        added variable layer_thickness 
  ! (needed in evolution.f90 to compute grain cross section including mantles)
  !        modified Nsites_grain = 5e4 (to have initial Nlayer compatible
  !                                     with mantle density of 1 g cm-3)

  IMPLICIT NONE
  INCLUDE "precision.f90"

  !--- physical characteristics - 1 ---
  REAL(KIND=DP)            :: amin_mrn                 ! grains MRN minimum radius (in cm)
  REAL(KIND=DP)            :: amax_mrn                 ! grains MRN maximum radius (in cm)
  REAL(KIND=DP)            :: alph_mrn                 ! grain MRN index
  REAL(KIND=DP)            :: rho_grc                  ! grain core volumic mass (g/cm3)
  REAL(KIND=DP)            :: rho_grm                  ! grain mantle volumic mass (g/cm3)
  REAL(KIND=DP)            :: f1_mrn                   ! 1st MRN integral functions
  REAL(KIND=DP)            :: f2_mrn                   ! 2nd MRN integral functions
  REAL(KIND=DP)            :: f3_mrn                   ! 3rd MRN integral functions
  REAL(KIND=DP)            :: f4_mrn                   ! 4th MRN integral functions

  !--- physical characteristics - 2 ---
  REAL(KIND=DP)            :: A_grc0                   ! Initial MRN normalisation factor
  REAL(KIND=DP)            :: A_grc                    ! MRN normalisation factor
  REAL(KIND=DP)            :: dens_grc0                ! Initial grain density (in cm-3)
  REAL(KIND=DP)            :: dens_grc                 ! Grain density (in cm-3)
  REAL(KIND=DP)            :: Mgrc0                    ! Initial grain core mass (per unit of gas volume in g cm-3)
  REAL(KIND=DP)            :: Mgrc                     ! Grain core mass (per unit of gas volume in g cm-3)
  REAL(KIND=DP)            :: Mgrm0                    ! Initial grain mantle mass (per unit of gas volume in g cm-3)
  REAL(KIND=DP)            :: Mgrm                     ! Grain mantle mass (per unit of gas volume in g cm-3)
  REAL(KIND=DP)            :: Mgrain                   ! Grain mass (per unit of gas volume in g cm-3)
  REAL(KIND=DP)            :: d_ero                    ! Erosion zone size (in cm)
  REAL(KIND=DP)            :: d_ads                    ! Adsorption zone size (in cm)
  REAL(KIND=DP)            :: d_site                   ! Distance between sites
  REAL(KIND=DP)            :: rsq_grc                  ! Mean square of core radius (in cm2)
  REAL(KIND=DP)            :: rsq_grm                  ! Mean square of mantle radius (in cm2)
  REAL(KIND=DP)            :: r_grm                    ! Root of the mean square of mantle radius (in cm)
  REAL(KIND=DP)            :: ab_cor                   ! Total abundance of core species (in cm-3)
  REAL(KIND=DP)            :: ab_ads                   ! Total abundance of adsorbed species (in cm-3)
  REAL(KIND=DP)            :: Nlayers                  ! number of layers in grain mantle

  REAL(KIND=DP)            :: ratio_GRAIN_gas          ! mass ratio (grains/gas)

  REAL(KIND=DP)            :: Tgrain                   ! grain temperature (K) read in READ_PARAMETERS
  REAL(KIND=DP)            :: Teff_grain               ! effective temperature for sputtering reactions

  REAL (KIND=dp)           :: rv                       ! RV = Av/E(B-V)
  REAL (KIND=dp)           :: cdunit                   ! NH/E(B-V)

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


END MODULE MODULE_GRAINS


MODULE MODULE_PHYS_VAR

  !*****************************************************************************
  !** The module 'MODULE_PHYS_VAR' contains all                               **
  !** the physical variables of the model, stored in 2 vectors :              **
  !**    * v_variab -> contains all physical variables in the order :         **
  !**                      (number in parenthesis is the number of variables) **
  !**                  1) MHD variables (Nv_MHD)                              **
  !**                  2) density of each chemical specy (Nspec)              **
  !**                  3) density of H2 levels (Nlevels_H2)                   **
  !**                                                                         **
  !**    * v_lvariab  = LOG(v_variab)                                         **
  !**                  this vector is sent into subroutine DRIVE              **
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  INTEGER(KIND=LONG) :: Nv_MHD  ! number of MHD variables (initialized in MHD)
  INTEGER(KIND=LONG) :: d_v_var ! dimension of the 2 vectors (calculated in MHD)

  !------------------------------------------------------------------
  ! this 2 vectors are allocated and initialized in MHD
  ! their dimension is
  ! d_v_var=(Nv_MHD + Nspec + Nlevels_H2)
  !------------------------------------------------------------------
  REAL(KIND=DP),DIMENSION(:), SAVE, ALLOCATABLE :: v_variab
  REAL(KIND=DP),DIMENSION(:), SAVE, ALLOCATABLE :: v_lvariab,v_dz_lvariab

  !-------------------------------------------------------
  ! physical variables of the model :
  ! initilised in INITIALIZE and calculated in DIFFUN
  ! (exept when precised)
  !-------------------------------------------------------

  !--- neutrals ---
  REAL(KIND=DP)      :: timeN = 0.0_DP     ! flow time (year) calculated in MHD
  REAL(KIND=DP)      :: RhoN               ! mass density (g.cm-3)
  REAL(KIND=DP)      :: DensityN           ! numeric density (cm-3)
  REAL(KIND=DP)      :: muN                ! mean mass (g)
  REAL(KIND=DP)      :: Tn                 ! temperature (K) read in READ_PARAMETERS
  REAL(KIND=DP)      :: Tn_old             ! temperature (K) at the previous call to drive
  REAL(KIND=DP)      :: nH_init            ! Initial proton density (cm-3), used in READ_SPECIES
  REAL(KIND=DP)      :: Vn                 ! velocity (cm/s)
  REAL(KIND=DP)      :: v_CH               ! velocity of CH (cm/s)
  REAL(KIND=DP)      :: v_S                ! velocity of S  (cm/s)
  REAL(KIND=DP)      :: v_SH               ! velocity of SH (cm/s)
  REAL(KIND=DP)      :: dVn                ! neutral velocity gradient
  REAL(KIND=DP)      :: compr_n            ! neutral compression factor
  REAL(KIND=DP)      :: coldens_h          ! H  column-density as a VODE variable
  REAL(KIND=DP)      :: coldens_h2         ! H2 column-density as a VODE variable
  REAL(KIND=DP)      :: coldens_co         ! CO column-density as a VODE variable

  !--- Used only for J shocks ---
  REAL(KIND=DP)      :: XLL                ! characteristic viscous length (cm)
  REAL(KIND=DP)      :: grad_V             ! neutral velocity gradient too (with viscosity...)
  REAL(KIND=DP)      :: save_dv            ! neutral velocity gradient too (without viscosity...)

  !--- positive ions ---
  REAL(KIND=DP)      :: timeI = 0.0_DP     ! flow time (year) calculated in MHD
  REAL(KIND=DP)      :: RhoI               ! mass density (g.cm-3)
  REAL(KIND=DP)      :: DensityI           ! numeric density (cm-3)
  REAL(KIND=DP)      :: muI                ! mean mass (g)
  REAL(KIND=DP)      :: Ti                 ! temperature (K)
  REAL(KIND=DP)      :: Vi                 ! velocity (cm/s)
  REAL(KIND=DP)      :: dVi                ! ions velocity gradient
  REAL(KIND=DP)      :: compr_i            ! ion     compression factor

  !--- negative ions ---
  REAL(KIND=DP)      :: RhoA               ! mass density (g.cm-3)
  REAL(KIND=DP)      :: DensityA           ! numeric density (cm-3)
  REAL(KIND=DP)      :: muA                ! mean mass (g)

  !--- negative ions and electrons ---
  REAL(KIND=DP)      :: RhoNEG             ! mass density (g.cm-3)
  REAL(KIND=DP)      :: DensityNEG         ! numeric density (cm-3)
  REAL(KIND=DP)      :: muNEG              ! mean mass (g)
  REAL(KIND=DP)      :: Te                 ! temperature (K) (for e- and negative ions)

  !--- charges including grains ---
  REAL(KIND=DP)      :: RhoCharges         ! mass density (g.cm-3)
  REAL(KIND=DP)      :: DensityCharges     ! numeric density (cm-3)
  REAL(KIND=DP)      :: muCharges          ! mean mass (g)

  !--- common ---
  REAL(KIND=DP)      :: nH                 ! density of protons nH=n(H)+2n(H2)+n(H+)
  REAL(KIND=DP)      :: Col_Dens_nH        ! column-density of protons col_Dens_nH=Col_Dens(H)+2Col_Dens(H2)+Col_Dens(H+)
  REAL(KIND=DP)      :: DeltaV             ! Vi-Vn (cm/s)
  REAL(KIND=DP)      :: ABS_DeltaV         ! ABS(Vi-Vn) (cm/s)
  REAL(KIND=DP)      :: Vgrad              ! abs. value of the gradient of Vn (km.s-1.cm-1) calc. in MHD
  REAL(KIND=DP)      :: DeltaVmin = 1.D3   ! (cm/s) if ABS_DeltaV < DeltaVmin : stop
  REAL(KIND=DP)      :: Vsound             ! sound velocity (cm/s)
  REAL(KIND=DP)      :: Vmagnet            ! magnetosonic velocity (cm/s)
  REAL(KIND=DP)      :: Valfven            ! alfven speed (cm/s)

  !--- shock ---
  CHARACTER(len=2)   :: shock_type_init    ! initial shock type ('S1', 'S2', 'P1', 'P2', 'C', or 'J')
  CHARACTER(len=2)   :: shock_type_lock    ! initial shock type ('S1', 'S2', 'P1', 'P2', 'C', or 'J')
  CHARACTER(len=2)   :: shock_type         ! current shock type ('S1', 'S2', 'P1', 'P2', 'C', or 'J')
  CHARACTER(len=2)   :: shock_type_fin     ! final   shock type ('S1', 'S2', 'P1', 'P2', 'C', 'C*', 'CJ', or 'J')
  INTEGER(KIND=LONG) :: Nfluids            ! number of fluids (neutrals, ions, e-)
  INTEGER(KIND=LONG) :: F_CH               ! Flag: compute v_CH or not
  INTEGER(KIND=LONG) :: F_S                ! Flag: compute v_S  or not
  INTEGER(KIND=LONG) :: F_SH               ! Flag: compute v_SH or not
  REAL(KIND=DP)      :: Vs_km              ! shock velocity (km/s) read in READ_PARAMETERS
  REAL(KIND=DP)      :: Vs_cm              ! shock velocity (cm/s) read in READ_PARAMETERS
  REAL(KIND=DP)      :: timeJ              ! shock age (years) read in READ_PARAMETERS
  LOGICAL            :: viscosity          ! .TRUE. if we need artificial viscosity
  LOGICAL            :: viscosity_lock     ! .TRUE. if we need artificial viscosity
  INTEGER            :: ieqth = 0          ! Flag : if 1, thermal balanced solved
  INTEGER            :: Cool_KN = 0        ! Flag : if 1, Kaufman & Neufeld cooling (H20 & CO)
  REAL(KIND=DP)      :: size_shock         ! size of the shock (cm)
  REAL(KIND=DP)      :: time_shock         ! time of the shock (yr)
  INTEGER            :: isize_shock        ! index at which the shock is cut

  INTEGER(KIND=LONG) :: F_SORT             ! Flag: sort reaction rates
  INTEGER(KIND=LONG) :: F_CONS             ! Flag: test conservation equations
  INTEGER            :: F_TGR              ! compute grain temperature or not

  !--- H2 levels and co ---
  REAL(KIND=DP)      :: op_H2              ! H2-ortho/para ratio
  REAL(KIND=DP)      :: op_H2_in           ! initial H2-ortho/para ratio read in READ_PARAMETERS
  INTEGER(KIND=LONG) :: NH2_lev_var        ! Number of levels of H2 included
  INTEGER(KIND=LONG) :: NH2_lev            ! Number of levels of H2 included
  INTEGER(KIND=LONG) :: NH2_lines_out      ! Max number of H2 lines in output file
  CHARACTER(LEN=4)   :: H_H2_flag          ! Which collision data do we use? Choices:
                                           !       DRF : Flower et al.
                                           !        MM : Martin & Mandy
                                           !      BOTH : DRF if possible and MM if not
  INTEGER            :: iforH2 = 1         ! Flag : H2 formation on grains 
                                           !   0: 1/3 of 4.4781 eV in internal energy (=> 17249 K) (Allen, 1999)
                                           !   1: Proportional to Boltzman distribution at 17249 K
                                           !   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)
                                           !   3: v = 6, J = 0,1
                                           !   4: fraction = relative populations at t, initialised as H2_lev%density
                                           !                 and changed during integration
  INTEGER            :: ikinH2 = 1         ! Flag : H2 formation energy released as kinetic energy
                                           !   1: 0.5 * (4.4781 - internal)
                                           !   2: Inf(1.4927 eV, 4.4781 - internal)
  INTEGER            :: pumpH2 = 1         ! Flag : H2 pumping by UV photons
                                           !   0: not included
                                           !   1: included
  INTEGER(KIND=LONG) :: NCO_lev            ! Number of levels of CO included (only to compute the photodissociation rate)

  !--- photodissociation rates ---
  REAL(KIND=DP)      :: photodiss_H2       ! H2 photodissociation rates (used in outputs)
  REAL(KIND=DP)      :: photodiss_CO       ! CO photodissociation rates (used in outputs)
  REAL(KIND=DP)      :: photoH2            ! H2 photodissociation rates
  REAL(KIND=DP)      :: photoCO            ! CO photodissociation rates
  LOGICAL            :: pdr_file_fnd       ! boolean to search for a PDR input file

  !--- magnetic field (GAUSS), read in READ_PARAMETERS ---
  REAL(KIND=DP)      :: Bfield
  REAL(KIND=DP)      :: Bbeta

  !--- wind, read in READ_PARAMETERS ---
  REAL(KIND=DP)      :: Z0                 ! initial launch point

  !--- environment, read in READ_PARAMETERS ---
  INTEGER            :: F_ISRF             ! radiation field spectrum - 1 = Mathis, 2 = Draine
  REAL(KIND=DP)      :: RAD                ! multipicative factor for the flux radiation
  REAL(KIND=dp)      :: fphsec             ! Secondary photons flux
  REAL(KIND=DP)      :: Av,Av0             ! extinction (magnitudes), local and initial
  REAL(KIND=DP)      :: N_H2_0,N_CO_0      ! initial column-densities 
  INTEGER            :: F_COUP_RAD         ! perform a full coupling with radiation field transfer 
                                           ! (to compute dissociation rates, desorption, ...)
  INTEGER            :: F_AV               ! integrate Av or not
  INTEGER            :: F_invAv            ! use grain coefficients to compute AV/NH (0) or scale 
                                           ! grain coefficient to reproduce inv_Av_fac (1)
  REAL(KIND=DP)      :: inv_Av_fac         ! Av conversion factor.
  REAL(KIND=DP)      :: conv_coldens       ! NH conversion factor
  REAL(KIND=DP)      :: Zeta               ! cosmic ray ionization rate (s-1)
  REAL(KIND=DP)      :: vturb              ! turbulent velocity (cm s-1, used for Doppler broadening in FGK)

  !--- distance (calculated in MHD) ---
  REAL(KIND=DP)      :: distance     = 0.0_DP ! distance from beginning of the shock (cm)
  REAL(KIND=DP)      :: distance_old = 0.0_DP
  REAL(KIND=DP)      :: dist_step    = 0.0_DP ! step (cm) between 2 consecutive calls to DRIVE

  !--- slabing ---
  INTEGER            :: slab

  !--- calculation step (calc. in MHD) ---
  INTEGER(KIND=LONG) :: counter = 0          ! counts the number of calls to DRIVE
  INTEGER(KIND=LONG) :: counter_down = 0     ! counts the number of calls to DRIVE for the down trajectory
  INTEGER(KIND=LONG) :: counter_high = 0     ! counts the number of calls to DRIVE for the high trajectory
  INTEGER(KIND=LONG) :: counter_main = 0     ! counts the number of calls to DRIVE for the main trajectory
  INTEGER(KIND=LONG) :: counter_main_old = 0 ! counts the number of calls to DRIVE for the main trajectory (save previous value)

  !--- save test quantities ---
  ! REAL(KIND=DP)      :: aux1_dvn_save_tot
  ! REAL(KIND=DP)      :: aux1_dvn_save_Sn
  ! REAL(KIND=DP)      :: aux1_dvn_save_Vn
  ! REAL(KIND=DP)      :: aux1_dvn_save_Bn
  ! REAL(KIND=DP)      :: aux1_dvn_save_mol
  ! REAL(KIND=DP)      :: aux2_dvn_save_tot
  ! REAL(KIND=DP)      :: aux2_dvn_save_Cs
  ! REAL(KIND=DP)      :: aux2_dvn_save_Vn


  !----------------------------------------------------------
  ! index allowing to extract variables from the 2 vectors
  ! initialized in MHD
  !----------------------------------------------------------
  INTEGER(KIND=LONG) :: iv_Vn, iv_Vi
  INTEGER(KIND=LONG) :: iv_vCH,iv_vS,iv_vSH
  INTEGER(KIND=LONG) :: iv_RhoN, iv_RhoI, iv_RhoA, iv_RhoNEG
  INTEGER(KIND=LONG) :: iv_Tn, iv_Ti, iv_Te
  INTEGER(KIND=LONG) :: iv_DensityN, iv_DensityI, iv_DensityA, iv_DensityNeg
  INTEGER(KIND=LONG) :: iv_gv
  INTEGER(KIND=LONG) :: iv_nh,iv_nh2,iv_nco
  INTEGER(KIND=LONG) :: iv_compr_n
  INTEGER(KIND=LONG) :: iv_compr_i
  INTEGER(KIND=LONG) :: iv_d_ero,iv_d_ads
  INTEGER(KIND=LONG) :: bv_speci,ev_speci,bv_specy
  ! These are obsolete in thermochemistry but still used in evolution.
  INTEGER(KIND=LONG) :: bv_neu,ev_neu
  INTEGER(KIND=LONG) :: bv_ion,ev_ion
  INTEGER(KIND=LONG) :: bv_ani,ev_ani
  INTEGER(KIND=LONG) :: bv_neg,ev_neg
  INTEGER(KIND=LONG) :: bv_gra,ev_gra
  INTEGER(KIND=LONG) :: bv_cor,ev_cor
  INTEGER(KIND=LONG) :: bv_H2_lev,ev_H2_lev

  !--- useful for outputs (used in WRITE_OUTPUTS) ---
  CHARACTER(len=2)  :: speci_out
  CHARACTER(len=7)  :: H2_out
  CHARACTER(len=10) :: line_out
  INTEGER           :: Npthdf5     ! maximal number of points in hdf5 file
  INTEGER           :: Nstep_w     ! number of steps between 2 outputs in ascii files
  INTEGER           :: Nstep_max   ! max. number of integration steps
  REAL(KIND=DP)     :: sum_H2
  CHARACTER(len=1)  :: flag_analysis ! (Y/N) whether you want the chemical analysis

  INTEGER           :: stop_code
  INTEGER           :: count_istat4 ! Failed integration counter

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE READ_PARAMETERS

    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     Read the shock parameters from input file.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     shock_type_init, Vs_km, Vs_cm, Tn     (MODULE_PHYS_VAR)
    !     op_H2, Zeta, RAD, Av, Bfield          (MODULE_PHYS_VAR)
    !     Nstep_max, Nstep_w                    (MODULE_PHYS_VAR)
    !     Tgrain                                (MODULE_GRAINS)
    !     XLL, Eps_V                            (MODULE_VAR_VODE)
    !     duration_max, length_max              (MODULE_VAR_VODE)
    !---------------------------------------------------------------------------

    USE MODULE_VAR_VODE,  ONLY : integ_type, duration_max, length_max, Eps_V
    USE MODULE_GRAINS,    ONLY : Tgrain, amin_mrn, amax_mrn, alph_mrn, rho_grc, rho_grm, &
                                 f1_mrn, f2_mrn, f3_mrn, f4_mrn
    USE MODULE_CONSTANTS, ONLY : YEARsec, screen, parsec
    USE MODULE_TOOLS,     ONLY : GET_FILE_NUMBER
    USE MODULE_TECHCONFIG

    IMPLICIT NONE
    INTEGER                     :: file_input

    !----------------------------------------------------
    ! Read input parameters
    !----------------------------------------------------

    !--- opening file ---
    file_input = GET_FILE_NUMBER()
    OPEN(file_input,file = inpfile,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    !--- input files ---
    READ (file_input, *)
    READ (file_input, "(A40)") modele
    READ (file_input, "(A40)") specfile
    READ (file_input, "(A40)") chemfile
    READ (file_input, "(A40)") h2exfile
    READ (file_input, "(A40)") gridfile

    !--- shock parameters ---
    READ (file_input, *)
    READ (file_input, "(A2)") shock_type_init
    READ (file_input, *) Nfluids
    READ (file_input, *) Bbeta
    READ (file_input, *) Vs_km
    READ (file_input, *) DeltaVmin
    READ (file_input, *) nH_init
    READ (file_input, *) Tn
    READ (file_input, *) op_H2_in

    !--- environment ---
    READ (file_input, *)
    READ (file_input, *) Zeta
    READ (file_input, *) F_ISRF
    READ (file_input, *) RAD
    READ (file_input, *) Av0
    READ (file_input, *) F_COUP_RAD
    READ (file_input, *) F_AV
    READ (file_input, *) F_invAv
    READ (file_input, *) inv_Av_fac
    READ (file_input, *) N_H2_0
    READ (file_input, *) N_CO_0
    READ (file_input, *) vturb

    !--- grain properties ---
    READ (file_input, *)
    READ (file_input, *) F_TGR
    READ (file_input, *) Tgrain
    READ (file_input, *) amin_mrn
    READ (file_input, *) amax_mrn
    READ (file_input, *) alph_mrn
    READ (file_input, *) rho_grc
    READ (file_input, *) rho_grm

    !--- excitation & cooling ---
    READ (file_input, *)
    READ (file_input, *) ieqth
    READ (file_input, *) Cool_KN
    READ (file_input, *) NH2_lev_var
    READ (file_input, *) NH2_lines_out
    READ (file_input, "(A4)") H_H2_flag
    READ (file_input, *) iforH2
    READ (file_input, *) ikinH2
    READ (file_input, *) pumpH2
    READ (file_input, *) NCO_lev

    !--- numerical parameters ---
    READ (file_input, *)
    READ (file_input, *) integ_type
    READ (file_input, *) Nstep_max
    READ (file_input, *) timeJ
    READ (file_input, *) duration_max
    READ (file_input, *) Eps_V
    READ (file_input, *) XLL

    !--- output specifications ---
    READ (file_input, *)
    READ (file_input, *) F_W_HDF5_STD
    READ (file_input, *) F_W_HDF5_CHE
    READ (file_input, *) F_W_ASCII
    READ (file_input, *) Npthdf5
    READ (file_input, *) Nstep_w
    READ (file_input, "(A2)") speci_out
    READ (file_input, "(A7)") H2_out
    READ (file_input, "(A10)") line_out
    READ (file_input, "(A1)") flag_analysis 

    !--- developer options ---
    READ (file_input, *)
    READ (file_input, *) F_SORT
    READ (file_input, *) F_CONS
    READ (file_input, *) F_CH
    READ (file_input, *) F_S
    READ (file_input, *) F_SH

    !--- Wind ---
    READ (file_input, *) Z0

    !--- closing file ---
    CLOSE(file_input)


    shock_type_lock = shock_type_init
    shock_type_fin  = shock_type_init
    viscosity_lock  = .TRUE.

    !----------------------------------------------------
    ! Adjust filename length and apply conversion factors
    !----------------------------------------------------
    modele     = TRIM(ADJUSTL(modele))
    specfile   = TRIM(ADJUSTL(specfile))
    chemfile   = TRIM(ADJUSTL(chemfile))
    Vs_cm      = 1.0e5_dp * Vs_km
    op_H2      = op_H2_in
    Bfield     = Bbeta * sqrt(nH_init) * 1.0e-6_dp
    length_max = duration_max * abs(Vs_cm) * YEARsec
    vturb      = vturb * 1.0e+5_dp 
    ! Conversion CR ionization rate -> secondary photon flux
    ! Hollenbach 2009 - equivalent to G0 = 1e-3 for zeta = 1e-17 (Gredel 1989)
    ! CAREFUL - The value G0 = 1e-3 seems to be wrong. Comparison of the 
    !           photodissociation / ionization rates by the ISRF and by 
    !           secondary photons given by Haeys show an equivalent G0
    !           of about 1e-5. This value is in agreement with the
    !           computation of secondary photons ionization of grains
    !           based on the yield => adopt 1e-5 instead of 1e-3
    ! Note : conversion Mathis field -> number of photons between 911 and 2400 A : 1.44e8_dp
    ! Note2 : not used anymore
    fphsec = 1.0e-5_dp * zeta / 1.0e-17_dp * 1.44e8_dp

    !----------------------------------------------------
    ! Compute MRN integral functions
    ! Le Bourlot et al. (1995) - Eq. A2
    !----------------------------------------------------
    f1_mrn = (amin_mrn**(1.0_dp-alph_mrn) - amax_mrn**(1.0_dp-alph_mrn)) / (1.0_dp-alph_mrn)
    f2_mrn = (amin_mrn**(2.0_dp-alph_mrn) - amax_mrn**(2.0_dp-alph_mrn)) / (2.0_dp-alph_mrn)
    f3_mrn = (amin_mrn**(3.0_dp-alph_mrn) - amax_mrn**(3.0_dp-alph_mrn)) / (3.0_dp-alph_mrn)
    f4_mrn = (amin_mrn**(4.0_dp-alph_mrn) - amax_mrn**(4.0_dp-alph_mrn)) / (4.0_dp-alph_mrn)

    !----------------------------------------------------
    ! Read Tn in species.in if shock_type ≠ 'S' or 'P'
    !----------------------------------------------------
    PRINT *, 'shock type=<'//shock_type_init//'>', 'Tn=',Tn
    IF ( (shock_type_init(1:1) /= 'S') ) THEN
       file_input = GET_FILE_NUMBER()
       OPEN (file_input, file = TRIM(data_dir)//TRIM(specfile), status = 'OLD', action = 'READ')
       READ (file_input,'(54X,F9.2)') Tn
       CLOSE (file_input)
       PRINT *, 'read Tn in species.in file: Tn =',Tn
    ENDIF

    !----------------------------------------------------
    ! Check the coherence of the parameters
    ! 1 - Number of fluids
    ! 2 - Radiation field coupling
    ! 3 - Av integration
    ! 4 - Velocity gradient
    !     - computed in J-type (dedicated equation)
    !     - computed in C-type shock with dVn/dz
    !     - fixed in S-type models
    !       we assume 1 km s-1 / pc
    ! 5 - Cooling mode
    !     If S-type model, use KN = 0 - no opacity 
    !     dependent cooling (e.g. CO,) can be used 
    !     in this case
    ! 6 - H2 pumping
    ! 7 - timescales
    !     If S- or P-type models, set timeJ to duration_max
    !----------------------------------------------------

    ! -----------------------
    ! 1 - Number of fluids
    ! -----------------------
    IF ( (shock_type_init(1:1) == "S") .OR. &
         (shock_type_init(1:1) == "P") .OR. &
         (shock_type_init(1:1) == "W") .OR. &
         (shock_type_init(1:1) == "J") ) THEN
       IF ( Nfluids /= 1 ) THEN
          PRINT *,' ERROR in the input parameters'
          PRINT *,' -> conflict between Nfluid and shock_type_init'
          PRINT *,'    Nfluid must = 1 for shock_type_init = S, P, W, or J'
          STOP
       ENDIF
       ! --------------------
       ! old
       ! --------------------
       ! Nfluid = 1
       ! PRINT *,' WARNING: S or P are single fluid structures'
       ! PRINT *,'          -> Nfluid ...... set to   1'
       ! --------------------
    ENDIF
    ! -----------------------
    ! 2 - Rad field coupling
    ! -----------------------
    IF (RAD == 0.0_dp .AND. gridfile == 'none') THEN
       IF ( F_COUP_RAD /= 0 ) THEN
          PRINT *,' ERROR in the input parameters'
          PRINT *,' -> conflict between radiation and F_COUP_RAD'
          PRINT *,'    F_COUP_RAD must = 0'
          PRINT *,'    if no radiation field is used (RAD = 0 or gridfile = none)'
          STOP
       ENDIF
       ! --------------------
       ! old - automatic
       ! --------------------
       ! F_COUP_RAD = 0
       ! PRINT *,' WARNING: no coupling possible with the radiation field'
       ! PRINT *,'          -> F_COUP_RAD .. set to   0'
       ! --------------------
    ENDIF
    IF (RAD /= 0.0_dp .AND. gridfile == 'none') THEN
       IF ( F_COUP_RAD == 2 ) THEN
          PRINT *,' ERROR in the input parameters'
          PRINT *,' -> conflict between radiation and F_COUP_RAD and gridfile'
          PRINT *,'    F_COUP_RAD must = 0 or 1'
          PRINT *,'    if no grid radiation field is used (gridfile = none)'
          STOP
       ENDIF
       ! --------------------
       ! old - automatic
       ! --------------------
       ! F_COUP_RAD = 1
       ! PRINT *,' WARNING: bad choice of coupling with the radiation field'
       ! PRINT *,'          -> F_COUP_RAD .. set to   1'
       ! --------------------
    ENDIF
    IF (gridfile /= 'none') THEN
       IF ( RAD /= 0.0_dp .OR. F_COUP_RAD /= 2 .OR. F_AV /= 0 ) THEN
          PRINT *,' ERROR in the input parameters'
          PRINT *,' -> conflict between gridfile and RAD or F_COUP_RAD or F_AV'
          PRINT *,'    RAD and F_AV must = 0 and F_COUP_RAD must = 2'
          PRINT *,'    if a grid radiation field is used (gridfile not set to none)'
          STOP
       ENDIF
       ! --------------------
       ! old - automatic
       ! --------------------
       ! RAD        = 0.0_dp
       ! F_COUP_RAD = 2
       ! F_AV       = 0
       ! PRINT *,' WARNING: a grid of radiation field is provided'
       ! PRINT *,'          -> RAD ......... set to 0.0'
       ! PRINT *,'          -> F_COUP_RAD .. set to   2'
       ! PRINT *,'          -> F_AV ........ set to   0'
       ! --------------------
    ENDIF
    ! -----------------------
    ! 3 - Av integration
    ! -----------------------
    IF (shock_type_init(1:1) =='S') THEN
       IF ( F_AV /= 0 ) THEN
          PRINT *,' ERROR in the input parameters'
          PRINT *,' -> conflict between F_AV and shock_type_init'
          PRINT *,'    F_AV must = 0 for shock_type_init = S'
          STOP
       ENDIF
       ! --------------------
       ! old - automatic
       ! --------------------
       ! F_AV = 0
       ! PRINT *,' WARNING: enforce F_AV = 0 for shock_type_init = S1 or S2'
       ! --------------------
    ENDIF
    IF(inv_Av_fac /= 0.0_dp) THEN
       conv_coldens = 1.0_dp / inv_Av_fac
    ELSE
       conv_coldens = 0.0_dp
    ENDIF
    ! -----------------------
    ! 4 - Velocity gradient
    ! -----------------------
    IF      ( shock_type_init(1:1) == "J") THEN
       grad_V    = 1.0e-2_dp * Vs_cm / XLL
    ELSE IF ( shock_type_init(1:1) == "C") THEN
       grad_V    = 1e-30_dp ! non-zero value
    ELSE IF ( shock_type_init(1:1) == "S" .OR. &
              shock_type_init(1:1) == "W" .OR. &
              shock_type_init(1:1) == "P" ) THEN
       grad_V = SQRT(1.0e+10_dp / (parsec*parsec))
    ENDIF
    PRINT *, "    grad_V =", grad_V
    ! -----------------------
    ! 5 - Cooling mode
    ! -----------------------
    IF ( (shock_type_init(1:1) == "S") .OR. &
         (shock_type_init(1:1) == "P") .OR. &
         (shock_type_init(1:1) == "W") ) THEN
       IF ( Cool_KN /= 0 ) THEN
          PRINT *,' ERROR in the input parameters'
          PRINT *,' -> conflict between Cool_KN and shock_type_init'
          PRINT *,'    Cool_KN must = 0 for shock_type_init = S or P'
          STOP
       ENDIF
       ! --------------------
       ! old - automatic
       ! --------------------
       ! Cool_KN = 0
       ! --------------------
    ENDIF
    ! -----------------------
    ! 6 - H2 pumping
    ! -----------------------
    IF (F_COUP_RAD == 0) THEN
       IF ( pumpH2 /= 0 ) THEN
          PRINT *,' ERROR in the input parameters'
          PRINT *,' -> conflict between F_COUP_RAD and pumpH2'
          PRINT *,'    pumpH2 must = 0 if F_COUP_RAD = 0'
          STOP
       ENDIF
       ! --------------------
       ! old - automatic
       ! --------------------
       ! pumpH2 = 0
       ! PRINT *,' WARNING: pumpH2 only if full coupling with the radiation field'
       ! PRINT *,'          -> pumpH2     .. set to   0'
       ! --------------------
    ENDIF
    !------------------------
    ! 7 - Timescale
    !------------------------
    IF ( shock_type_init(1:1) == 'S' ) THEN
       timeJ = duration_max
       PRINT *,'nH=', nH, 'duration_max=', duration_max
    ENDIF
    ! -----------------------
    ! 8 - Wind
    ! -----------------------
    IF ( (shock_type_init(1:1) /= 'W') .AND. (Z0 /= 0.0_dp) ) THEN
       PRINT *,' ERROR in the input parameters'
       PRINT *,' -> conflict between Shock type and Z0'
       PRINT *,'    Z0 must = 0 if Shock type = S, P, J, or C'
       STOP
    ENDIF

    !----------------------------------------------------
    ! Check H2 parameters
    !   - number of H2 levels computed is always >= 2
    !   - It does not make sense to treat only one 
    !     independent H2 level
    !----------------------------------------------------
    NH2_lev = max(2, NH2_lev_var)
    IF (NH2_lev_var == 1) THEN
       NH2_lev_var = 0 
    ENDIF
    IF (iforH2 < -1 .OR. iforH2 > 4) THEN
       PRINT *, "  iforH2 =", iforH2, "  is not allowed"
       STOP
    ENDIF
    IF (ikinH2 /= 1 .AND. ikinH2 /= 2) THEN
       PRINT *, "  ikinH2 =", ikinH2, "  is not allowed"
       STOP
    ENDIF
    IF (H_H2_flag /= "DRF" .AND. H_H2_flag /= "MM" .AND. H_H2_flag /= "BOTH") THEN
       WRITE(screen,'(" ")')
       WRITE(screen,'("  Wrong flag for H-H2 collisions: ",A4)') H_H2_flag
       STOP
    ENDIF

    !----------------------------------------------------
    ! Prepare output files
    !----------------------------------------------------
    out_dir = TRIM(out_dir)//TRIM(ADJUSTL(modele))//'/'
    sh_cmd = 'sh -c "if [ ! -d '//TRIM(out_dir)//" ]; then mkdir "
    sh_cmd = TRIM(sh_cmd)//" "//TRIM(out_dir)//'; fi "'
    CALL SYSTEM(sh_cmd)

  END SUBROUTINE READ_PARAMETERS


  SUBROUTINE PDR_PHOTO_DAT

    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER

    IMPLICIT none
    CHARACTER(LEN=28),           PARAMETER   :: prefpdr = 'input/pdr_grid_photo/2S_a1p1'
    CHARACTER(LEN=23),           PARAMETER   :: suffpdr = 'x0p0m1p0xray-pdrdat.res'
    CHARACTER(LEN=200)                       :: pdrname
    INTEGER                                  :: file_input
    CHARACTER(LEN=1)                         :: a,b
    INTEGER                                  :: bb
    INTEGER                                  :: nb,ios
    INTEGER                                  :: i
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_av
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_temp
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_cd_h2
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_cd_co
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_n_h2
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_n_co
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_ph_h2
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pdr_ph_co
    REAL(KIND=DP)                            :: slope

    pdr_file_fnd = .true.

    !--- opening file ---
    file_input=GET_FILE_NUMBER()
    ! ----
    IF(nH_init == 0) bb=0
    IF(nH_init /= 0) bb = FLOOR(LOG10(nH_init))
    WRITE(b,'(I1)') abs(bb)
    WRITE(a,'(I1)') NINT(nH_init/10**DBLE(bb))
    IF (bb >= 0) pdrname = 'n'//a//'p'//b
    IF (bb <  0) pdrname = 'n'//a//'m'//b
    IF (RAD > 0.0_dp) THEN 
       bb = FLOOR(LOG10(RAD))
    ENDIF
    WRITE(b,'(I1)') ABS(bb)
    WRITE(a,'(I1)') NINT(RAD/10**DBLE(bb))
    IF (bb >= 0) pdrname = TRIM(ADJUSTL(pdrname))//'r'//a//'p'//b
    IF (bb <  0) pdrname = TRIM(ADJUSTL(pdrname))//'r'//a//'m'//b
    bb = FLOOR(LOG10(Zeta/1e-17))
    WRITE(b,'(I1)') bb
    WRITE(a,'(I1)') NINT(Zeta/1e-17/10**DBLE(bb))
    IF (bb >= 0) pdrname = TRIM(ADJUSTL(pdrname))//'z'//a//'p'//b
    IF (bb <  0) pdrname = TRIM(ADJUSTL(pdrname))//'z'//a//'m'//b
    ! ----
    pdrname = TRIM(ADJUSTL(prefpdr))//TRIM(ADJUSTL(pdrname))//TRIM(ADJUSTL(suffpdr))
    INQUIRE(FILE=TRIM(ADJUSTL(pdrname)), EXIST=pdr_file_fnd) 
    IF( .NOT. pdr_file_fnd ) RETURN
    OPEN(file_input,file=TRIM(ADJUSTL(pdrname)),status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    !--- read number of data in the pdr output file ---
    nb = 0
    READ(file_input,*)
    READ(file_input,*)
    READ(file_input,*)
    DO 
        READ(file_input,*,iostat=ios)
        if ( ios /= 0 ) exit
        nb = nb + 1
    END DO

    !--- allocate local tables ---
    ALLOCATE(pdr_av   (nb))
    ALLOCATE(pdr_temp (nb))
    ALLOCATE(pdr_cd_h2(nb))
    ALLOCATE(pdr_cd_co(nb))
    ALLOCATE(pdr_n_h2 (nb))
    ALLOCATE(pdr_n_co (nb))
    ALLOCATE(pdr_ph_h2(nb))
    ALLOCATE(pdr_ph_co(nb))

    !--- read pdr output data ---
    REWIND(file_input)
    READ(file_input,*)
    READ(file_input,*)
    READ(file_input,*)
    DO i = 1, nb
       READ (file_input,*) pdr_av(i), pdr_temp(i), pdr_cd_h2(i), pdr_cd_co(i), &
                           pdr_n_h2 (i), pdr_n_co (nb), pdr_ph_h2(i), pdr_ph_co(i)
    ENDDO
    CLOSE(file_input)

    !--- find the closest Av point in the pdr output file ---
    i = 1
    DO WHILE( pdr_av(i) <= Av0)
       i = i + 1
       IF(i == nb + 1) EXIT
    ENDDO
    IF (i > nb + 1) THEN
       photoH2 = pdr_ph_h2(nb)
       photoCO = pdr_ph_co(nb)
    ELSE
       slope   = ( pdr_ph_h2(i) - pdr_ph_h2(i-1) ) / ( pdr_av(i) - pdr_av(i-1) )
       photoH2 = pdr_ph_h2(i-1) + slope * ( Av0 - pdr_av(i-1) )
       slope   = ( pdr_ph_co(i) - pdr_ph_co(i-1) ) / ( pdr_av(i) - pdr_av(i-1) )
       photoCO = pdr_ph_co(i-1) + slope * ( Av0 - pdr_av(i-1) )
    ENDIF

    !--- deallocate all (hopefully) local tables ---
    DEALLOCATE(pdr_av   )
    DEALLOCATE(pdr_temp )
    DEALLOCATE(pdr_cd_h2)
    DEALLOCATE(pdr_cd_co)
    DEALLOCATE(pdr_n_h2 )
    DEALLOCATE(pdr_n_co )
    DEALLOCATE(pdr_ph_h2)
    DEALLOCATE(pdr_ph_co)

  END SUBROUTINE PDR_PHOTO_DAT

END MODULE MODULE_PHYS_VAR


