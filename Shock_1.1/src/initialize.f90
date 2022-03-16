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

MODULE MODULE_INITIALIZE
  !*****************************************************************************
  !** The module 'MODULE_INITIALIZE' contains subroutine for reading          **
  !** and initialize :                                                        **
  !**      * the TechConfig files (properties of HDF5 output files)           **
  !**      * the shock parameters                                             **
  !**      * the chemical species                                             **
  !**      * the chemical reactions                                           **
  !**      * the H2 molecule (levels, collision rates, lines)                 **
  !**      * the SiO molecule (collision rates, Einstein coefficients)        **
  !*****************************************************************************
  ! SC May 06: 
  ! - rg: Dens_GRAIN_init now calculated correctly from mass in cores only
  ! - ng: check that Dens(G,G+,G-) is consistent with Dens_GRAIN, else stop
  ! -     Automatic initialization of Dens(G,G+,G-) in the STEADY run case
  ! -     do not include mantles in grain/gas(H) mass ratio
  ! - lay: R_gr_scale12 = <rgrain>/<rsquare_grain> for current size distrib
  ! - G: C60 species replaced by G 

  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE READ_METADATA_NBL(filename, nhdf5)
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    read the number of metadata written in the files stored in 
    !    input/TechConfig directory in order to allocate the metatables.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    number of lines in the metadata files
    !---------------------------------------------------------------------------

    USE MODULE_TECHCONFIG
    USE MODULE_CONSTANTS, ONLY : iwrtmp
    USE MODULE_TOOLS, ONLY     : GET_FILE_NUMBER
    IMPLICIT NONE

    CHARACTER(len=lenfilename), INTENT(in)  :: filename  ! directory of the ASCII files
    INTEGER,                    INTENT(out) :: nhdf5     ! dimension of the metatab table
    CHARACTER(LEN=5)                        :: templ

    ! ------------------------------------------------------------------
    ! Open the inpt metadata file and count the number of lines in it
    ! ------------------------------------------------------------------
    iwrtmp = GET_FILE_NUMBER()
    OPEN(iwrtmp,file=TRIM(filename),status='old')
    READ(iwrtmp,'(12/)')
    templ = ""
    nhdf5 = 0
    DO WHILE (templ /= '#----' )
       READ(iwrtmp,'(A5)') templ
       nhdf5 = nhdf5 + 1
    ENDDO
    nhdf5 = nhdf5 -1
    CLOSE(iwrtmp)

  END SUBROUTINE READ_METADATA_NBL



  SUBROUTINE READ_METADATA_DAT(filename, metatab, nhdf5)
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    read the list of data to write in the standard hdf5 file along with
    !    their metadata. All information are stored in input/TechConfig
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    number of lines in the metadata files
    !---------------------------------------------------------------------------
  
    USE MODULE_TECHCONFIG
    USE MODULE_CONSTANTS, ONLY : iwrtmp
    USE MODULE_TOOLS, ONLY     : GET_FILE_NUMBER
    IMPLICIT NONE
  
    CHARACTER(len=lenfilename),                   INTENT(in)  :: filename  ! directory of the ASCII files
    INTEGER,                                      INTENT(in)  :: nhdf5     ! dimension of the metatab table
    TYPE (METADATA),            DIMENSION(nhdf5), INTENT(out) :: metatab   ! table of metadata
    INTEGER, PARAMETER                                        :: nchar=1700
    CHARACTER(LEN=nchar)                                      :: line_tmp
    INTEGER                                                   :: i, j, k
    INTEGER                                                   :: ki, kf
  
    ! ------------------------------------------------------------------
    ! Read the inpt metadata file and store information in metatab
    ! ------------------------------------------------------------------
    iwrtmp = GET_FILE_NUMBER()
    OPEN(iwrtmp,file=TRIM(filename),status='old')
    READ(iwrtmp,'(12/)')
    DO i = 1, nhdf5
       READ(iwrtmp,'(A1700)') line_tmp ! double check that the number of char we 
                                       ! read is consistent with the value of 
                                       ! *nchar* above...
       j  = 1
       ki = 1
       DO k = 1, nchar
          IF( line_tmp(k:k) /= ';' ) CYCLE
          kf = k - 1
          IF(j== 1) metatab(i)%path   = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 2) metatab(i)%file   = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 3) metatab(i)%datset = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 4) metatab(i)%Hname  = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 5) metatab(i)%IDname = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 6) metatab(i)%dtype  = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 7) metatab(i)%unit   = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 8) metatab(i)%SKOS   = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j== 9) metatab(i)%ucd    = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j==10) metatab(i)%utype  = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j==11) metatab(i)%group  = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j==12) metatab(i)%parent = TRIM(ADJUSTL(line_tmp(ki:kf)))
          IF(j==12) metatab(i)%descrp = TRIM(ADJUSTL(line_tmp(kf+2:nchar)))
          ki = kf + 2
          j = j + 1
       ENDDO
    ENDDO
    CLOSE(iwrtmp)
  
  END SUBROUTINE READ_METADATA_DAT



  SUBROUTINE INITIALIZE
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    reads shock parameters, chemical species and reactions,
    !    initializes all physical variables
    ! subroutine/function needed :
    !    READ_PARAMETERS
    !    INITIALIZE_ELEMENTS
    !    READ_SPECIES
    !    CHECK_SPECIES
    !    READ_REACTIONS
    !    REACTION_TYPE
    !    CHECK_REACTIONS
    !    ADD_REVERSE_REACTIONS
    !    READ_H2_LEVELS
    !    INITIALIZE_ROVIB_H2
    !    READ_H2_RATES
    !    READ_H2_LINES
    !    READ_H2_ELECTRONIC
    !    READ_SiO_RATES
    !    EINTEIN_COEFF_SiO
    !    READ_FE_DATA
    !    RF_INIT
    !    ALLOCATE_TABLES
    ! input variables :
    ! ouput variables :
    ! results :
    !   almost all physical variables are initialized here
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_CHEM_REACT
    USE MODULE_CHEMICAL_SPECIES, ONLY : density_limit, READ_SPECIES
    USE MODULE_GAMMA
    USE MODULE_PHYS_VAR
    USE MODULE_PROFIL_TABLES
    USE MODULE_VAR_VODE,       ONLY : Tout_V
    USE MODULE_GRAINS
    USE MODULE_CONSTANTS,      ONLY : pi, mP, me, amu
    USE MODULE_H2,             ONLY : READ_H2_LEVELS, LLEV_H2, INITIALIZE_ROVIB_H2, &
                                      READ_H2_RATES, READ_H2_LINES, READ_H2_ELECTRONIC,&
                                      TYPE_H2_LEVEL, H2_lev
    USE MODULE_SiO,            ONLY : READ_SiO_RATES, EINSTEIN_COEFF_SiO
    USE MODULE_CO,             ONLY : READ_CO_LEVELS, BOLTZMANN_ROVIB_CO, READ_CO_ELECTRONIC, &
                                      CO_lev
    USE MODULE_DEBUG_JLB
    USE MODULE_LINE_EXCIT
    USE MODULE_READ_FE_DATA
    USE MODULE_DUST_TREATMENT, ONLY : READ_DRAINE_FILES, COMPUTE_QCOEF, COMPUTE_DERO, COMPUTE_DADS
    USE MODULE_RADIATION

    IMPLICIT NONE

    INTEGER(KIND=LONG) :: i
    REAL(KIND=DP)      :: dum         ! dummy variable

    REAL (KIND=dp)     :: dumy        ! dumy variable for calculating the doppler weighted H2 column density
    REAL (KIND=dp)     :: bt_s_bd     ! dumy variable for calculating the doppler weighted H2 column density
    REAL (KIND=dp)     :: cdinf       ! dumy variable for calculating the doppler weighted H2 column density

    ! mkdir output
    CALL system('mkdir -p output')
    !------------------------
    !--- shock parameters ---
    !------------------------
    WRITE(*,'("parameters : ")',ADVANCE='NO')
    CALL READ_PARAMETERS
    WRITE(*,'(" read. ")')

    !-----------------------------------------------------
    ! initialize compression factors
    !-----------------------------------------------------
    compr_n   = 1.0_DP
    compr_i   = 1.0_DP

    !---------------------------------------------------------------
    !--- H2 and CO photodissociation rates from PDR code outputs ---
    !--- used if label PDR16 is used in the chemical input file  ---
    !---------------------------------------------------------------
    CALL PDR_PHOTO_DAT

    !=====================================================================================
    ! 1 - initialize chemical species
    !=====================================================================================
    ! --- elemental abundances from Ander & Grevesse 1989 ---
    ! -> names, masses (g), abundances (cm-3)
    CALL INITIALIZE_ELEMENTS

    ! read in file_input + add 5 species + initialize indexes ...
    WRITE(*,'("species : ")',ADVANCE='NO')
    CALL READ_SPECIES
    WRITE(*,'(" ",I3," + ",I1," added ... ")',ADVANCE='NO') Nspec, Nspec_plus - Nspec

    ! nH= proton density (cm-3)
    nH = SUM(speci(:Nspec)%density*dble(speci(:Nspec)%formula(ind_elem_H)))
    ! check the set of species
    WRITE(*,'("check.")')
    CALL CHECK_SPECIES

    ! DO i = 1, Nspec
    !    WRITE(*,'(A8,1X,12(I2),3X,I2)') speci(i)%name, speci(i)%formula(1:Nelements), speci(i)%charge
    ! ENDDO
    ! WRITE(*,*)
    ! DO i = 1, Nelements
    !    WRITE(*,'(A8,1X,3(ES9.2,1X))') elements(i)%name, elements(i)%ab, elements(i)%ab_init, elements(i)%ab_ref
    ! ENDDO
    ! STOP

    !-----------------------------------------------------
    ! Check the indices of the different families
    ! BG2017
    !-----------------------------------------------------
    ! WRITE(*,*) Nneutrals, b_neu, e_neu
    ! WRITE(*,*) Nions, b_ion, e_ion
    ! WRITE(*,*) Nani, b_ani, e_ani
    ! WRITE(*,*) Nneg, b_neg, e_neg
    ! WRITE(*,*) Nongrains, b_gra, e_gra
    ! WRITE(*,*) Noncores, b_cor, e_cor
    ! STOP

    !-----------------------------------------------------
    ! Initialize Av and self-shielding column densities
    ! Modif - set to 1e1_DP initially - BG2017 (test)
    !-----------------------------------------------------
    Av = Av0
    coldens_h  = 1e1_DP
    coldens_h2 = 1e1_DP
    coldens_co = 1e1_DP

    !-----------------------------------------------------
    ! Grains properties - 1
    !-----------------------------------------------------
    ! - Compute the initial grain core mass (per unit of
    !   gas volume) Mgrc0 based on the ** species initial 
    !   abundances. This mass is supposed to be spread on 
    !   grain cores following a MRN distribution from 
    !   amin_mrn to amax_mrn
    ! - Use Mgrc0 to compute A_grc0, the MRN normalisation
    !   factor, and dens_grc0, the density of grain core
    !   (or incidentally of grains) per unit of gas volume
    ! - Initialize d_ero (erosion zone size) to 1E-30_DP
    !   (something small but not 0)
    ! - Initialize rsq_grc
    !-----------------------------------------------------
    Mgrc0     = DOT_PRODUCT( speci(b_cor:e_cor)%density, speci(b_cor:e_cor)%mass )! / amu
    Mgrc      = Mgrc0
    A_grc0    = Mgrc0 * 3.0_DP / ( 4 * pi * rho_grc * f4_mrn )! * amu
    A_grc     = A_grc0 * compr_i
    dens_grc0 = A_grc0 * f1_mrn
    dens_grc  = dens_grc0 * compr_i
    ab_cor    = SUM( speci(b_cor:e_cor)%density )
    IF ( (shock_type_init(1:1) == "S") .OR. &
         (shock_type_init(1:1) == "P") .OR. &
         (shock_type_init(1:1) == "W") ) THEN
       d_ero = 0.0_DP
    ELSE
       CALL COMPUTE_DERO
    ENDIF
    rsq_grc   = ( f3_mrn - 2.0_DP * f2_mrn * d_ero + f1_mrn * d_ero**2.0_DP ) / f1_mrn

    !-----------------------------------------------------
    ! Grains properties - 2
    !-----------------------------------------------------
    ! - Compute the initial grain mantles mass (per unit 
    !   of gas volume) Mgrm0 based on the * species 
    !   initial abundances.
    ! - Initialize total grain mass - Mgrain = Mgrc + Mgrm
    ! - Initialize d_ads (adsorption zone size) based on
    !   the * species initial abundances (NR algorithm)
    !   Two formalisms are possible depending on whether
    !   we consider that mantles are built on the eroded
    !   grains (more realistic) or not (simpler)
    ! - Initialize rsq_grm
    !   Again two formalisms are possible
    !-----------------------------------------------------
    Mgrm0   = DOT_PRODUCT( speci(b_gra:e_gra)%density, speci(b_gra:e_gra)%mass )! / amu
    Mgrm    = Mgrm0
    Mgrain  = Mgrc + Mgrm
    ab_ads  = SUM( speci(b_gra:e_gra)%density )
    CALL COMPUTE_DADS
    d_site  = ( Mgrm / ab_ads / rho_grm )**(1.0_DP/3.0_DP)
    Nlayers = d_ads / d_site
    ! 1st formalism
    rsq_grm = ( f3_mrn + 2.0_DP * f2_mrn * d_ads + f1_mrn * d_ads**2.0_DP ) / f1_mrn
    r_grm   = sqrt(rsq_grm)
    ! 2nd formalism
    ! rsq_grm  = ( f3_mrn &
    !            + 2.0_DP * f2_mrn * ( d_ads - d_ero ) &
    !            + f1_mrn * ( d_ero**2.0_DP + d_ads**2.0_DP - 2.0_DP * d_ero * d_ads ) ) / f1_mrn
    ! r_grm    = sqrt(rsq_grm)

    WRITE(*,*)
    WRITE(*,'(A35,ES10.3)') "Erosion zone size ............... ", d_ero
    WRITE(*,'(A35,ES10.3)') "Adsorption zone size ............ ", d_ads
    WRITE(*,'(A35,ES10.3)') "Initial core square radius ...... ", rsq_grc
    WRITE(*,'(A35,ES10.3)') "Initial mantles square radius ... ", rsq_grm
    WRITE(*,'(A35,ES10.3)') "Grain initial mass .............. ", Mgrain
    WRITE(*,*)

    ratio_GRAIN_gas = Mgrc / (1.4 * nH * mP)

    !-----------------------------------------------------
    ! Grains properties - 3
    !-----------------------------------------------------
    ! - Check the compatibility between G0, G-, and G+
    !   abundances and that derived from ** species above
    ! - renormalize grain density G0 + G- + G+ = dens_grc
    ! - initialize G0, G+, and G- mass
    !-----------------------------------------------------
    dum = speci(ind_G)%density + speci(ind_Gplus)%density + speci(ind_Gminus)%density

    PRINT *, ''
    PRINT *, 'Grains'
    PRINT *, '------'
    PRINT *, 'G  ', speci(ind_G)%density
    PRINT *, 'G+ ', speci(ind_Gplus)%density
    PRINT *, 'G- ', speci(ind_Gminus)%density
    PRINT *, '  -> n(G) + n(G+) + n(G-) =', dum
    PRINT *, '  -> Total ng from ** species   =', dens_grc

    ! check that input grain number densities match total grain density
    IF (ABS(dum/dens_grc - 1._DP) > 5e-2) then 
      PRINT *, ' Grain number does not match grain mass and size distribution'
      PRINT *, ' Run with shock type "S", Nfluids=1 to update species.in'
      ! STOP
    ENDIF

    ! renormalisation of grain densities
    speci(ind_G)%density      = speci(ind_G)%density      * dens_grc / dum
    speci(ind_Gplus)%density  = speci(ind_Gplus)%density  * dens_grc / dum
    speci(ind_Gminus)%density = speci(ind_Gminus)%density * dens_grc / dum

    ! recompute electron density
    speci(Nspec)%density = SUM( speci(1:Nspec-1)%density * speci(1:Nspec-1)%charge )
    speci(Nspec)%density_init = speci(Nspec)%density
    IF (speci(Nspec)%density < density_limit) THEN
       speci(Nspec)%density = density_limit
    ENDIF

    !!!! !-----------------------------------------------------
    !!!! ! initialize G0, G+, and G- mass
    !!!! ! Try to add the mass of the grains in the computation
    !!!! ! of mass and momentum fluxes between different flows
    !!!! ! This development doesn't work so far - BG2017
    !!!! ! (two other parts are commented in evolution.f90
    !!!! !-----------------------------------------------------
    !!!!
    !!!! elements(ind_elem_G)%mass = Mgrain / dens_grc
    !!!! speci(ind_G)%mass         = DOT_PRODUCT( DBLE(elements(1:Nelements)%mass), &
    !!!!                                          DBLE(speci(ind_G)%formula(1:Nelements)) ) &
    !!!!                           - speci(ind_G)%charge * me
    !!!! speci(ind_Gplus)%mass     = DOT_PRODUCT( DBLE(elements(1:Nelements)%mass), &
    !!!!                                          DBLE(speci(ind_Gplus)%formula(1:Nelements)) ) &
    !!!!                           - speci(ind_Gplus)%charge * me
    !!!! speci(ind_Gminus)%mass    = DOT_PRODUCT( DBLE(elements(1:Nelements)%mass), &
    !!!!                                          DBLE(speci(ind_Gminus)%formula(1:Nelements)) ) &
    !!!!                           - speci(ind_Gminus)%charge * me
    !!!! mass_G                    = speci(ind_G)%mass
    !!!! mass_Gplus                = speci(ind_Gplus)%mass 
    !!!! mass_Gminus               = speci(ind_Gminus)%mass

    !=====================================================================================
    ! 2 - initialize fluids
    !=====================================================================================

    !-----------------------------------------------------
    ! number density (cm-3) of each fluid = sum(density)
    !-----------------------------------------------------
    DensityN   = SUM(DBLE(speci(b_neu:e_neu)%density)) ! neutrals
    DensityI   = SUM(DBLE(speci(b_ion:e_ion)%density)) ! ions > 0
    DensityA   = SUM(DBLE(speci(b_ani:e_ani)%density)) ! ions < 0
    DensityNEG = SUM(DBLE(speci(b_neg:e_neg)%density)) ! species < 0 (including electrons)

    !-----------------------------------------------------
    ! mass density (g cm-3) of each fluid = sum(mass*dens)
    ! include the contributions of the grains
    !-----------------------------------------------------
    RhoN   = DOT_PRODUCT(DBLE(speci(b_neu:e_neu)%density), DBLE(speci(b_neu:e_neu)%mass))
    RhoI   = DOT_PRODUCT(DBLE(speci(b_ion:e_ion)%density), DBLE(speci(b_ion:e_ion)%mass))
    RhoA   = DOT_PRODUCT(DBLE(speci(b_ani:e_ani)%density), DBLE(speci(b_ani:e_ani)%mass))
    RhoNEG = DOT_PRODUCT(DBLE(speci(b_neg:e_neg)%density), DBLE(speci(b_neg:e_neg)%mass))

    !-----------------------------------------------------
    ! mean mass (g) of each fluid mu = rho / density
    !-----------------------------------------------------
    muN   = RhoN   / DensityN
    muI   = RhoI   / DensityI
    muA   = RhoA   / DensityA
    muNEG = RhoNEG / DensityNEG

    !-----------------------------------------------------
    ! temperatures : all the same at the beginning
    !-----------------------------------------------------
    Ti = Tn
    Te = Tn

    !-----------------------------------------------------
    ! velocities (cm/s) of neutrals and ions
    !-----------------------------------------------------
    Vn = Vs_cm
    IF (Nfluids == 1) THEN
       DeltaV = 0.0_DP
    ELSE
       ! = Vi - Vn, to start the shock
       DeltaV = DeltaVmin
    ENDIF
    ABS_DeltaV = ABS(DeltaV)
    Vi = Vn - DeltaV

    !-----------------------------------------------------
    ! velocity gradient (km.s-1.cm-1)
    !-----------------------------------------------------
    ! Vgrad = Vs_km/distance
    Vgrad = Vs_km / Tout_V

    !-----------------------------------------------------
    ! Compute Vmagnet, Vsound, and Valfven
    !    Add grains in charges density and mass for the
    !    computation of Vmagnet and Vi. This is not done 
    !    in RhoI, RhoNEG, DensityI or DensityNEG because 
    !    momentum and energy exchanges are based on 
    !    ions-neutral collision cross sections - BG2017
    !-----------------------------------------------------
    DensityCharges = DensityI &
                   + DensityNEG &
                   + speci(ind_G)%density
    RhoCharges     = RhoI + RhoNEG + Mgrain
    muCharges      = RhoCharges / DensityCharges

    Vmagnet        = Gamma * kB * DensityCharges * (Ti+Te) / RhoCharges &
                   + ( Bfield * Vs_cm / Vi )**2._DP / ( 4._DP * pi * RhoCharges )
    Vmagnet        = SQRT(Vmagnet)
    Vsound         = SQRT ( Gamma * kB * Tn / muN )
    Valfven        = SQRT ( (Bfield * Vs_cm / Vi)**2._DP / (4._DP * pi * RhoN) )
    !-----------------------------------------------------
    ! Debug
    !-----------------------------------------------------
    ! WRITE(*,'(11(A,ES10.2,/))') &
    !                    "Vmagnet = ", Vmagnet,          &
    !                    "Gamma   = ", Gamma,            &
    !                    "Ti+Te   = ", (Ti+Te),          &
    !                    "n+*mu+  = ", DensityI*muI,     &
    !                    "n-*mu-  = ", DensityNEG*muNEG, &
    !                    "n+      = ", DensityI,         &
    !                    "n-      = ", DensityNeg,       &
    !                    "B*Vs/Vi = ", Bfield*Vs_cm/Vi,  &
    !                    "RhoI    = ", RhoI,             &
    !                    "RhoNeg  = ", RhoNEG,           &
    !                    "Mgrain  = ", Mgrain
    !-----------------------------------------------------

    !-----------------------------------------------------
    ! Check entrance velocity compare to sound and 
    ! magnetosonic velocity. Stop the code or set up 
    ! a J-type shock
    !-----------------------------------------------------
    PRINT *, 'Vmagnet=', Vmagnet, ' Vsound=', Vsound, ' Valfven=', Valfven
    IF ( shock_type_init == 'C' ) THEN
       !---------------------------------------
       ! Check we are in the proper range for
       ! C-shocks: Vsound < Vn < Vmagnet
       !---------------------------------------
       ! 1 - if Vn > Vmagnet -> set a J shock
       !---------------------------------------
       IF ( Vmagnet < Vn ) THEN
          PRINT *, 'Entrance speed Vn=', Vn, ' is bigger than Vmagnet...'
          shock_type_lock = "J"
          shock_type_fin  = "J"
          viscosity_lock  = .TRUE.
          Nfluids         = 1
          grad_V          = 1.0e-2_dp * Vs_cm / XLL
       ENDIF
       !---------------------------------------
       ! 2 - if Vn < Vsound  -> stop the code
       !---------------------------------------
       IF ( Vn < Vsound ) THEN
          PRINT*,'Entrance speed Vn=',Vn,' is sub-sonic...'
          STOP
       ENDIF
    ENDIF

    !=====================================================================================
    ! 3 - initialize chemical reactions
    !=====================================================================================

    WRITE(*,'("reactions : ")',ADVANCE='NO')
    ! read in file_chemistry
    CALL READ_REACTIONS
    WRITE(*,'(I4," ... ")',ADVANCE='NO') Nreact

    ! find reaction type and re-order the entire set
    ! calculate DE according to the type of the reaction
    CALL REACTION_TYPE

    ! check if the set of reactions is correct
    ! to do after REACTION_TYPE, because this subroutine needs react%type
    WRITE(*,'("check ... ")', ADVANCE='NO')
    CALL CHECK_REACTIONS

    ! compute energy defect. Update Tabone 09/18. DE was previously computed in READ_TYPE
    ! (without compution DE for REVERSE_REACTIONS). Now: loop of ALL REACTIONS
    ! => now, it's easier to check for which reaction DE is computed
    ! ENERGY_DEFEC ahs to be called here so that DE can be used by ADD_REVERSE_REACTIONS
    CALL ENERGY_DEFECT_REACTION
    
    ! add the missing endothermic reactions
    IF (do_we_add_reactions) CALL ADD_REVERSE_REACTIONS
    WRITE(*,'("+",I3," added.")') Nrever

    ! compute energy defect. Update Tabone 09/18. DE was previously computed in READ_TYPE
    ! (without compution DE for REVERSE_REACTIONS). Now: loop of ALL REACTIONS
    ! => now, it's easier to check for which reaction DE is computed
    CALL ENERGY_DEFECT_REACTION

    ! Read photoreaction whose rates are integrated over their cross sections
    CALL READ_PHOTODES

    ! Initialize tables for computing the photoelectric rates
    CALL INITIALIZE_PHOTOELE

    !=====================================================================================
    ! 4 - initialize specific molecules (H2, SiO, CO)
    !=====================================================================================

    !-----------------------------------------------------
    !---                      H2                      ----
    !-----------------------------------------------------

    WRITE(*,'("H2 molecule : initialization.")')

    ! H2 levels
    CALL READ_H2_LEVELS

    ! ortho / para levels indices (important to reduce 
    ! the computational time - look FGKDATA subroutines)
    CALL LLEV_H2

    ! population of the levels, according to the choice in ortho:para
    CALL INITIALIZE_ROVIB_H2

    ! read collision rates for H-H2, He-H2, H2-H2
    CALL READ_H2_RATES

    ! read quadrupolar lines of H2
    CALL READ_H2_LINES

    ! read H2 electronic transition data (+initialization)
    CALL READ_H2_ELECTRONIC

    !-----------------------------------------------------
    !---                      SiO                     ----
    !-----------------------------------------------------
    WRITE(*,'("SiO molecule : initilization.")')

    ! collision rates for SiO-H2
    CALL READ_SIO_RATES

    ! Einstein coefficients for SiO
    CALL EINSTEIN_COEFF_SiO

    ! read in atomic data for Fe+
    CALL READ_FE_DATA

    !-----------------------------------------------------
    !---                      CO                      ----
    !-----------------------------------------------------
    WRITE(*,'("CO molecule : initilization.")')

    ! CO levels
    CALL READ_CO_LEVELS

    ! Initialize populations of CO levels
    CALL BOLTZMANN_ROVIB_CO

    ! read CO electronic transition data (+initialization)
    CALL READ_CO_ELECTRONIC

    !=====================================================================================
    ! 5 - initialize lines emissivities ans integrated intensities
    !=====================================================================================
    emihat    = 0.0_DP
    emihat_o  = 0.0_DP
    inthat    = 0.0_DP
    emicat    = 0.0_DP
    emicat_o  = 0.0_DP
    intcat    = 0.0_DP
    eminat    = 0.0_DP
    eminat_o  = 0.0_DP
    intnat    = 0.0_DP
    emioat    = 0.0_DP
    emioat_o  = 0.0_DP
    intoat    = 0.0_DP
    emisat    = 0.0_DP
    emisat_o  = 0.0_DP
    intsat    = 0.0_DP
    emisiat   = 0.0_DP
    emisiat_o = 0.0_DP
    intsiat   = 0.0_DP
    emicpl    = 0.0_DP
    emicpl_o  = 0.0_DP
    intcpl    = 0.0_DP
    eminpl    = 0.0_DP
    eminpl_o  = 0.0_DP
    intnpl    = 0.0_DP
    emiopl    = 0.0_DP
    emiopl_o  = 0.0_DP
    intopl    = 0.0_DP
    emispl    = 0.0_DP
    emispl_o  = 0.0_DP
    intspl    = 0.0_DP
    emisipl   = 0.0_DP
    emisipl_o = 0.0_DP
    intsipl   = 0.0_DP
    emifepl   = 0.0_DP
    emifepl_o = 0.0_DP
    intfepl   = 0.0_DP

    !-----------------------------------------------------
    !---          Set the names of transition          ---
    !-----------------------------------------------------
    CALL BUILD_TR_NAME

    !=====================================================================================
    ! 6 - initialize interstellar radiation field and the radiation field grid
    !=====================================================================================
    CALL READ_RAD_GRID

    CALL READ_DRAINE_FILES

    CALL INIT_WAVELENGTH

    CALL INIT_RADIATION

    CALL INTERP_SIGMA_G_RF

    !-----------------------------------------------------
    ! Compute the grains coefficients interpolate over
    ! the wavelength grid. These coef are used for
    ! - compute the transfert (if F_RAD_COUP = 1)
    ! - compute the photoejection rate from grains
    ! - compute the grain temperatures
    !-----------------------------------------------------
    CALL COMPUTE_QCOEF
    CALL INTERP_QABSO(nwlg,qabso_D_rf,wlg)

    !-----------------------------------------------------
    ! Set the radiation field used to compute the
    ! - pumping of H2
    ! - the photodestruction rates
    ! - the photoejction rate from grains
    ! - the grain temperature
    ! - the photoelectric effect
    !
    ! IF F_COUP_RAD = 1
    !    => Compute the initial optical depth
    !    => do not include gas absorption because the 
    !       column densities in the buffer are not known
    ! IF F_COUP_RAD = 2
    !    => The radiation field is set to that
    !       given in the grid file
    !-----------------------------------------------------
    IF      (F_COUP_RAD == 1) THEN
       doptdpth(1:nwlg) = pi * rsq_grm * qabso_D_rf(1:nwlg) * dens_grc * conv_coldens * Av / nH
       ! optdpth(1:nwlg)  = optdpth_old(1:nwlg) + doptdpth(1:nwlg)
       CALL SIMPLE_TRANSFER
    ELSE IF (F_COUP_RAD == 2) THEN
       CALL SET_RAD_GRID(1)
    ENDIF

    !-----------------------------------------------------
    ! Compute H2 and CO UV lines absorption coefficients
    ! used to calculate the pumping of H2 and the photo-
    ! dissociation of H2 and CO
    ! WARNING: 
    ! - pumping coefficients of H2 and dissociation rates
    !   of H2 and CO are calculated later when diffun is
    !   called
    !-----------------------------------------------------
    cdinf = 1.e23_dp

    bt_s_bd = vturb / SQRT(vturb*vturb + 2.0_dp * kB * Tn / speci(ind_H2)%mass)
    dumy = SUM(H2_lev(1:NH2_lev)%density)
    DO i = 1, NH2_lev
       H2_lev(i)%cd_l     = N_H2_0 * bt_s_bd * H2_lev(i)%density / dumy
       H2_lev(i)%cd_r     = cdinf  * bt_s_bd * H2_lev(i)%density / dumy
       H2_lev(i)%cd_l_old = H2_lev(i)%cd_l
    ENDDO

    bt_s_bd = vturb / SQRT(vturb*vturb + 2.0_dp * kB * Tn / speci(ind_CO)%mass)
    dumy = SUM(CO_lev(1:NCO_lev)%density)
    DO i = 1, NCO_lev
       CO_lev(i)%cd_l     = N_CO_0 * bt_s_bd * CO_lev(i)%density / dumy
       CO_lev(i)%cd_r     = cdinf  * bt_s_bd * CO_lev(i)%density / dumy
       CO_lev(i)%cd_l_old = CO_lev(i)%cd_l
    ENDDO

    CALL FGKCOEF

    !-----------------------------------------------------
    ! Compute the rates of
    ! - photodestruction (for species with cross section)
    ! - heating by photodestruction (Tabone 01/18)
    ! WARNING: 
    ! - photoejection rate from grains & grain temperature
    !   are computed in diffun
    !-----------------------------------------------------
    DO i = 1, Ncross
       phdest_rate(i)    = SECT_INT(phdest(i))
       phheating_rate(i) = HEATING_PHOTOCHEM_SECT_INT(phdest(i))
       ! -- debug --
       ! WRITE(*,*) speci(phdest(i)%ispe)%name, phheating_rate(i), phdest_rate(i)
    ENDDO

    !=====================================================================================
    ! 7 - allocate table to save the shock profils
    !=====================================================================================
    CALL ALLOCATE_TABLES
    CALL INITIALIZE_TABLES(Nstep_max, traj_main)
    CALL INITIALIZE_TABLES(Nstep_max, traj_high)
    CALL INITIALIZE_TABLES(Nstep_max, traj_down)
    CALL INITIALIZE_TABLES(Nstep_max, traj_curr)

  END SUBROUTINE INITIALIZE

END MODULE MODULE_INITIALIZE
