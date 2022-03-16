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

MODULE MODULE_DRIVE

  IMPLICIT NONE
  INCLUDE "precision.f90"

  !-------------------------------------------------------
  ! values saved before call to DRIVE
  !-------------------------------------------------------
  REAL (KIND=DP)            :: Vn_old          = 0.0_DP
  REAL (KIND=DP)            :: Vi_old          = 0.0_DP
  REAL (KIND=DP)            :: dVn_old         = 0.0_DP
  REAL (KIND=DP)            :: dVi_old         = 0.0_DP
  REAL (KIND=DP)            :: old_Tn          = 0.0_DP
  REAL (KIND=DP)            :: Tinit           = 0.0_DP ! Detects whether Tn decreases.
  LOGICAL                   :: increase        = .false.
  LOGICAL                   :: maximum_reached = .false.

  !-------------------------------------------------------
  ! initial time scales
  !-------------------------------------------------------
  REAL (KIND=DP)            :: time_scale
  REAL (KIND=DP)            :: tiny = 1.0e-40_DP
  INTEGER                   :: i_time

  !-------------------------------------------------------
  ! Computation timescales
  !-------------------------------------------------------
  REAL (KIND=DP)            :: time1
  REAL (KIND=DP)            :: time2

  !-------------------------------------------------------
  ! Variables for CJ-type shocks
  !-------------------------------------------------------
  ! find jump conditions
  LOGICAL                   :: cj_active   = .false.
  LOGICAL                   :: cs_active   = .false.
  REAL (KIND=DP)            :: t_min_jump
  REAL (KIND=DP)            :: t_max_jump
  REAL (KIND=DP)            :: t_jump
  INTEGER                   :: i_min_jump
  INTEGER                   :: i_max_jump
  INTEGER                   :: i_jump
  ! adjust trajectory
  LOGICAL                   :: adjust_traj  = .false.
  LOGICAL                   :: adjust_start = .false.
  LOGICAL                   :: test_sonic   = .false.
  LOGICAL                   :: sonic_point  = .false.
  LOGICAL                   :: found_high   = .false.
  LOGICAL                   :: found_down   = .false.
  REAL (KIND=DP), PARAMETER :: eps_traj = 1.e-2_dp
  REAL (KIND=DP)            :: t_limit
  INTEGER                   :: i_limit
  INTEGER                   :: i_adjst                ! QUESTION - REMOVE i_adjst ? -> ONLY DEAL WITH i_limit
  INTEGER                   :: i_down
  INTEGER                   :: nadj = 0
  REAL (KIND=DP)            :: alpha
  REAL (KIND=DP)            :: angle_traj = 0.0_dp
  REAL (KIND=DP)            :: dum
  REAL (KIND=DP)            :: dum1
  REAL (KIND=DP)            :: dum2
  REAL (KIND=DP)            :: coef1
  REAL (KIND=DP)            :: coef2
  REAL (KIND=DP)            :: Vn_high
  REAL (KIND=DP)            :: Vn_down
  REAL (KIND=DP)            :: Vi_high
  REAL (KIND=DP)            :: Vi_down
  REAL (KIND=DP)            :: step_cross = 0.0_dp
  REAL (KIND=DP)            :: sonic_jump = 0.0_dp
  ! recoupling of fluids
  LOGICAL                   :: fluid_decoup = .false.
  LOGICAL                   :: fluid_recoup = .false.
  REAL (KIND=DP)            :: decoup_dist  = 0.0_dp
  REAL (KIND=DP)            :: decoup_strg  = 0.0_dp
  REAL (KIND=DP)            :: corr_flux_dum= 0.0_dp

  INTEGER                   :: countersave

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


CONTAINS


  SUBROUTINE RESET_INTEGRATION(is, again)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    reset integration for the next slab
     ! subroutine/function needed :
     ! input variables :
     !    is    - integer of current slabl
     !    again - boolean
     ! ouput variables :
     !    again - boolean : false if the shock code must stop
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS
     USE MODULE_PHYS_VAR
     USE MODULE_GRAINS
     USE MODULE_H2
     USE MODULE_CHEMICAL_SPECIES
     USE MODULE_CHEM_REACT
     USE MODULE_DUST_TREATMENT
     USE MODULE_RADIATION
     USE MODULE_EVOLUTION
     USE MODULE_VAR_VODE

     IMPLICIT none

     INTEGER,        INTENT(in)    :: is
     LOGICAL,        INTENT(inout) :: again
     INTEGER                       :: i

     LOGICAL,        SAVE          :: first=.true.

     ! -------------------------------------------------
     ! set up the origin of distances
     ! -------------------------------------------------
     IF (first) THEN
        first    = .false.
        distance = posgrd(is-1)
     ENDIF

     ! -------------------------------------------------
     ! set the type of shock to run
     ! -------------------------------------------------
     IF ( posgrd(is-1) < 0 ) THEN
        shock_type = "S1"
     ELSE
        shock_type = shock_type_lock
     ENDIF
     IF      ( shock_type(1:1) == "S" .OR.&
               shock_type(1:1) == "P" .OR.&
               shock_type(1:1) == "W" ) THEN
        Nfluids    = 1
        viscosity  = .FALSE.
     ELSE IF ( shock_type(1:1) == "J" ) THEN
        Nfluids    = 1
        viscosity  = viscosity_lock
     ELSE IF ( shock_type(1:1) == "C" ) THEN
        Nfluids    = 3
        viscosity  = .FALSE.
     ENDIF

     ! -------------------------------------------------
     ! IF F_COUP_RAD = 2
     !    Set the radiation field to its value at the
     !    current slab and reset the FGK coef and
     !    photodestruction rates with this radiation
     ! -------------------------------------------------
     IF (F_COUP_RAD == 2) THEN
        CALL SET_RAD_GRID(is)
        CALL FGKCOEF

        DO i = 1, Ncross
           phdest_rate(i)    = SECT_INT(phdest(i))
           phheating_rate(i) = HEATING_PHOTOCHEM_SECT_INT(phdest(i))
          ! -- debug --
          ! WRITE(*,*) speci(phdest(i)%ispe)%name, phheating_rate(i), phdest_rate(i)
        ENDDO
     ENDIF

     ! -------------------------------------------------
     ! Reset dvode arguments and call diffun
     ! -------------------------------------------------
     again         = .TRUE.
     IF      (integ_type == 0) THEN
        T0_V         = 0.0_dp                   ! distance : value of the integration variable
        distance_old = 0.0_dp                   ! distance old
        Tout_V       = posgrd(is)-distance      ! step to reach
     ELSE IF (integ_type == 1) THEN
        T0_V         = 0.0_dp                   ! distance : value of the integration variable
        distance_old = EXP(0.0_dp)              ! distance old
        Tout_V       = LOG(posgrd(is)-distance) ! step to reach
     ENDIF

     Hnext         = 0.0d0
     H0_V          = XLL * 1.d-5             ! step length - not really useful
     MF_V          = 22                      ! chosen numerical procedure
     Itask_V       = 1                       ! Which Task DVODE should do
     Istate_V      = 1                       ! How was it done (2 on normal return, 1 to initialise dvode)
     count_istat4  = 0

     CALL diffun(d_v_var,T0_v,v_lvariab,v_dvariab)
     dVn = v_dvariab(iv_Vn)
     dVi = v_dvariab(iv_Vi)

     ! -------------------------------------------------
     ! Info on chemical rates
     ! -------------------------------------------------
     ! PRINT *,'Species name,   density,    Chemical rate (1/s):'
     ! DO i = 1, Nspec
     !    PRINT *, speci(i)%name, speci(i)%density, YN(i)
     ! ENDDO
     !-------------------------------------------------------

     !-------------------------------------------------------
     ! Print info on initial length scales
     !-------------------------------------------------------
     PRINT *
     time_scale =  1.0_DP / MAXVAL( abs(v_dvariab(1:d_v_var) ) + tiny )
     i_time = MAXLOC( abs(v_dvariab), dim = 1 )
     PRINT *,' Length scale estimates AT START:'
     PRINT "(' lenght_scale=',1pe13.3,' cm')", time_scale
     PRINT "(' Max_var number i=',i4)", i_time
     IF (i_time.ge.bv_H2_lev.and.i_time.le.ev_H2_lev) THEN
        PRINT *,' this corresponds to H2 level number ',i_time-bv_H2_lev+1
     ENDIF
     IF (i_time.ge.bv_speci.and.i_time.le.ev_speci) THEN
        PRINT *,' this corresponds to species ',speci(i_time-bv_speci+1)%name
     ENDIF
     PRINT "(' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
          exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
     PRINT *
     IF      (integ_type == 0) THEN
        rwork_V(5) = MIN(1e-4*ABS(time_scale),1e-4_dp) ! H0_V           ! step size on first step
        rwork_V(6) = Tout_V                            ! 0 !1.0D17      ! HMAX
        rwork_V(7) = MIN(1e-8*ABS(time_scale),1e-8_dp) ! H0_V * 1.0D-10 ! HMIN
     ELSE IF (integ_type == 1) THEN
        rwork_V(5) = 1e-4_dp   ! H0_V           ! step size on first step
        rwork_V(6) = 1.0_dp    ! 0 !1.0D17      ! HMAX
        rwork_V(7) = 1e-16_dp  ! H0_V * 1.0D-10 ! HMIN
     ENDIF

  END SUBROUTINE RESET_INTEGRATION



  SUBROUTINE TEST_CONSERVATION_LAWS
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    test if the elements, the charges, and the level populations
     !    of H2 are conserved. Force conservation laws if necessary
     !      In subroutines MAIN_CARRIER & INDEX_CONSERVATION
     !      the main carrier of all elements are calculated
     !      CONSISTENCY solves a system where the main 
     !      carrier are adjusted as much as necessary to
     !      insure conservation of elements
     !      - IN PROGRESS - NOT SATISFACTORY CRUDE PATCH
     !      - ALSO ADD CHARGE CONSERVATION
     !      DOESNT WORK YET
     ! subroutine/function needed :
     ! input variables :
     ! ouput variables :
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_PHYS_VAR
     USE MODULE_CHEMICAL_SPECIES
     USE MODULE_H2

     IMPLICIT none

     CHARACTER(LEN=8)                :: fail
     REAL(KIND=DP),    PARAMETER     :: epscon = 1e-05_dp
     LOGICAL                         :: cons
     REAL(KIND=DP)                   :: dum
     REAL (KINd=DP)                  :: densHini
     REAL (KINd=DP)                  :: neucompr
     REAL (KINd=DP)                  :: ioncompr
     REAL(KIND=DP)                   :: somme
     INTEGER                         :: i
     INTEGER                         :: ii

     IF      ( F_CONS == 0 ) THEN
        RETURN
     ELSE IF ( F_CONS == 1 ) THEN
        IF ( ABS( v_variab (bv_specy + ind_H2) - SUM( v_variab(bv_H2_lev:ev_H2_lev) ) ) &
             / v_variab (bv_specy + ind_H2) > epscon) THEN
           v_variab (bv_specy + ind_H2) = SUM( v_variab(bv_H2_lev:ev_H2_lev) )
           v_lvariab(bv_specy + ind_H2) = LOG( v_variab(bv_specy + ind_H2)   )
        ENDIF

        IF ( ABS( v_variab (bv_specy + Nspec) - SUM( v_variab(bv_specy+1:bv_specy+Nspec-1)*speci(1:Nspec-1)%charge ) ) &
             / v_variab (bv_specy + Nspec) > epscon) THEN
           v_variab (bv_specy + Nspec) = SUM( v_variab(bv_specy+1:bv_specy+Nspec-1)*speci(1:Nspec-1)%charge )
           v_lvariab(bv_specy + Nspec) = LOG( v_variab(bv_specy + Nspec)   )
        ENDIF

        RETURN
     ENDIF

     ! -------------------------------------------
     ! all conservation equation are relative to
     ! the total H initial abundance
     ! -------------------------------------------
     densHini = elements(ind_elem_H)%Dens_init

     ! -------------------------------------------
     ! neutral and ion compression factors
     ! -------------------------------------------
     neucompr = v_variab(iv_compr_n)
     ioncompr = v_variab(iv_compr_i)

     ! -------------------------------------------
     ! Check if the chemical abundances verify
     ! conservation equations
     ! - H2 populations
     ! - elemental abundances
     !   (careful with differential compression)
     ! - abundances of PAHs
     ! - charge
     ! -------------------------------------------
     cons = .TRUE.
     fail = ""

     somme = SUM( H2_lev(1:NH2_lev)%density )
     dum   = ABS( speci(ind_H2)%density - somme ) / speci(ind_H2)%density
     IF ( dum  > epscon ) THEN
        cons = .FALSE.
        fail = "pop H2  "
     ENDIF
     ! WRITE(*,'(A9,2(A8),3(ES10.3,2X))') "coucou b ", "popH2   ", "H2      ", &
     !                                    speci(ind_H2)%density, somme, dum

     DO i = 1, Nelements
        IF ( elements(i)%icon == 0 ) CYCLE
        somme = 0.0_dp
        DO ii = 1, Nspec
           IF ( speci(ii)%fluid == 1 ) THEN
              somme = somme + DBLE(speci(ii)%formula(i)) * speci(ii)%density / neucompr
           ELSE
              somme = somme + DBLE(speci(ii)%formula(i)) * speci(ii)%density / ioncompr
           ENDIF
        ENDDO
        somme = somme / densHini
        dum   = ABS( elements(i)%ab_init - somme ) / elements(i)%ab_init
        IF ( dum  > epscon ) THEN
           cons = .FALSE.
           fail = elements(i)%name
        ENDIF
        ! WRITE(*,'(A9,2(A8),3(ES10.3,2X))') "coucou c ", elements(i)%name, speci(elements(i)%icon)%name, &
        !                                    elements(i)%ab_init, somme, dum
     ENDDO

     somme = speci(ind_PAH0)%density / neucompr  &
           + speci(ind_PAHp)%density / ioncompr  &
           + speci(ind_PAHm)%density / ioncompr
     somme = somme / densHini
     dum   = ABS( PAH_abinit - somme ) / PAH_abinit
     IF ( dum  > epscon ) THEN
        cons = .FALSE.
        fail = "PAH"
     ENDIF
     ! WRITE(*,'(A9,2(A8),3(ES10.3,2X))') "coucou d ", "PAH     ", speci(PAH_icon)%name, &
     !                                    PAH_abinit, somme, dum

     somme = SUM( speci(1:Nspec-1)%density * speci(1:Nspec-1)%charge )
     dum   = ABS( speci(Nspec)%density - somme ) / speci(Nspec)%density
     IF ( dum  > epscon ) THEN
        cons = .FALSE.
        fail = "charge"
     ENDIF

     IF( fail /= "" ) THEN
        ! WRITE(*,*) "conservation law broken for ... ", fail
        ! READ(*,*)
     ENDIF

     ! -------------------------------------------
     ! If one of the conservation equation is 
     ! broken, find the main carriers and the 
     ! species used for conservations and apply 
     ! the consistency procedure
     ! -------------------------------------------
     IF ( .NOT.cons ) THEN
        ! CALL MAIN_CARRIER
        ! CALL INDEX_CONSERVATION
        CALL CONSISTENCY

        ! WRITE(*,*)
        ! DO i = 1, Nspec
        !    WRITE(*,'(A8,2X,3(ES10.3,2X))') speci(i)%name, speci(i)%density, v_variab (bv_specy+i), &
        !                                    abs(speci(i)%density - v_variab (bv_specy+i)) / v_variab (bv_specy+i)
        ! ENDDO

        ! -------------------------------------------
        ! update other chemistry dependent
        ! physical variables
        ! -------------------------------------------
        DensityN   = SUM(DBLE(speci(b_neu:e_neu)%density))
        DensityI   = SUM(DBLE(speci(b_ion:e_ion)%density))
        DensityA   = SUM(DBLE(speci(b_ani:e_ani)%density))
        DensityNEG = SUM(DBLE(speci(b_neg:e_neg)%density))
        RhoN       = DOT_PRODUCT(DBLE(speci(b_neu:e_neu)%density), DBLE(speci(b_neu:e_neu)%mass))
        RhoI       = DOT_PRODUCT(DBLE(speci(b_ion:e_ion)%density), DBLE(speci(b_ion:e_ion)%mass))
        RhoA       = DOT_PRODUCT(DBLE(speci(b_ani:e_ani)%density), DBLE(speci(b_ani:e_ani)%mass))
        RhoNEG     = DOT_PRODUCT(DBLE(speci(b_neg:e_neg)%density), DBLE(speci(b_neg:e_neg)%mass))

        ! -------------------------------------------
        ! update dvode variables and take logs
        ! -------------------------------------------
        v_variab(bv_speci:ev_speci) = speci(1:Nspec)%density
        v_variab(iv_DensityN  )     = DensityN
        v_variab(iv_DensityI  )     = DensityI
        v_variab(iv_DensityA  )     = DensityA
        v_variab(iv_DensityNeg)     = DensityNeg
        v_variab(iv_RhoN      )     = RhoN
        v_variab(iv_RhoI      )     = RhoI
        v_variab(iv_RhoA      )     = RhoA
        v_variab(iv_RhoNEG    )     = RhoNEG

        v_lvariab(bv_speci:ev_speci)  = LOG( v_variab(bv_speci:ev_speci) )
        v_lvariab(iv_DensityN  )      = LOG( v_variab(iv_DensityN  )     )
        v_lvariab(iv_DensityI  )      = LOG( v_variab(iv_DensityI  )     )
        v_lvariab(iv_DensityA  )      = LOG( v_variab(iv_DensityA  )     )
        v_lvariab(iv_DensityNeg)      = LOG( v_variab(iv_DensityNeg)     )
        v_lvariab(iv_RhoN      )      = LOG( v_variab(iv_RhoN      )     )
        v_lvariab(iv_RhoI      )      = LOG( v_variab(iv_RhoI      )     )
        v_lvariab(iv_RhoA      )      = LOG( v_variab(iv_RhoA      )     )
        v_lvariab(iv_RhoNEG    )      = LOG( v_variab(iv_RhoNEG    )     )
     ENDIF

  END SUBROUTINE TEST_CONSERVATION_LAWS



  SUBROUTINE CONSISTENCY
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    solves a system where the abundances of the main carriers of all 
     !    elements, charge, PAHs, and H2 populations are adjusted in order 
     !    to insure conservation laws
     ! subroutine/function needed :
     !    DGESVX
     ! input variables :
     ! output variables :
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_PHYS_VAR
     USE MODULE_CHEMICAL_SPECIES
     USE MODULE_H2

     IMPLICIT none

     LOGICAL,        DIMENSION(Nelements)        :: skipele
     INTEGER                                     :: c_nvar
     LOGICAL,        ALLOCATABLE, DIMENSION(:)   :: sp_var
     INTEGER,        ALLOCATABLE, DIMENSION(:)   :: c_ivar
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: c_am
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: c_bbmn
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: c_bbmx
     INTEGER,        ALLOCATABLE, DIMENSION(:)   :: c_indx
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: c_AF
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)   :: c_R
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)   :: c_C
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)   :: c_FERR
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)   :: c_BERR
     REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)   :: c_WORK
     INTEGER,        ALLOCATABLE, DIMENSION(:)   :: c_IWORK
     REAL (KIND=dp)                              :: c_RCOND
     CHARACTER (LEN=1)                           :: c_FACT  = 'E'
     CHARACTER (LEN=1)                           :: c_TRANS = 'N'
     CHARACTER (LEN=1)                           :: c_EQUED
     INTEGER                                     :: c_nrhs  = 1
     INTEGER                                     :: c_INFO
     INTEGER                                     :: i, k, kk, ll
     
     REAL (KINd=DP)                              :: densHini
     REAL (KINd=DP)                              :: neucompr
     REAL (KINd=DP)                              :: ioncompr
     REAL (KINd=DP)                              :: somme

     ! -------------------------------------------
     ! all conservation equation are relative to
     ! the total H abundance
     ! -------------------------------------------
     densHini = elements(ind_elem_H)%Dens_init

     ! -------------------------------------------
     ! neutral and ion compression factors
     ! -------------------------------------------
     neucompr = v_variab(iv_compr_n)
     ioncompr = v_variab(iv_compr_i)

     ! -------------------------------------------
     ! Count the number of elements to include
     ! -------------------------------------------
     c_nvar = 0
     DO k = 1, Nelements
        ! IF ( k == ind_elem_Mg .OR. k == ind_elem_Fe ) CYCLE
        ! IF ( elements(k)%icon == 0 ) THEN
        IF ( elements(k)%icon == 0 .OR. k == ind_elem_Mg .OR. k == ind_elem_Fe .OR. k == ind_elem_H ) THEN
           skipele(k) = .true.
        ELSE
           skipele(k) = .false.
           c_nvar = c_nvar + 1
        ENDIF
     ENDDO

     ! -------------------------------------------
     ! add conservation of H2 pop, charge and PAH
     ! -------------------------------------------
     c_nvar = c_nvar + 3
 
     ! -------------------------------------------
     ! Allocate tables to solve the system
     ! -------------------------------------------
     ALLOCATE ( sp_var (1:Nspec)           )
     ALLOCATE ( c_ivar (1:c_nvar)          )
     ALLOCATE ( c_am   (1:c_nvar,1:c_nvar) )
     ALLOCATE ( c_bbmn (1:c_nvar,1:c_nrhs) )
     ALLOCATE ( c_bbmx (1:c_nvar,1:c_nrhs) )
     ALLOCATE ( c_indx (1:c_nvar)          )
     ALLOCATE ( c_AF   (1:c_nvar,1:c_nvar) )
     ALLOCATE ( c_R    (1:c_nvar)          )
     ALLOCATE ( c_C    (1:c_nvar)          )
     ALLOCATE ( c_FERR (1:1)               )
     ALLOCATE ( c_BERR (1:1)               )
     ALLOCATE ( c_WORK (1:4*c_nvar)        )
     ALLOCATE ( c_IWORK(1:c_nvar)          )
 
     ! -------------------------------------------
     ! Find the species used for each conserv law
     ! - H2 population => use H2
     ! - charge        => use electron
     ! -------------------------------------------
     c_ivar(1:c_nvar) = 0
     sp_var(1:Nspec)  = .FALSE.

     kk = 1
     c_ivar(kk) = ind_H2
     sp_var(ind_H2) = .TRUE.
     DO k = 1, Nelements
        IF ( skipele(k) ) THEN
           CYCLE
        ELSE
           kk = kk + 1
           c_ivar(kk) = elements(k)%icon
           sp_var(elements(k)%icon) = .TRUE.
        ENDIF
     ENDDO
     kk = kk + 1
     c_ivar(kk) = PAH_icon
     sp_var(PAH_icon) = .TRUE.
     kk = kk + 1
     c_ivar(kk) = charge_icon
     sp_var(charge_icon) = .TRUE.

     ! -------------------------------------------
     ! Write the system to solve
     ! -------------------------------------------
     c_am   (1:c_nvar,1:c_nvar) = 0.0_dp
     c_bbmn (1:c_nvar,1:c_nrhs) = 0.0_dp

     kk = 1
     c_am(kk,kk)  = 1.0_dp
     c_bbmn(kk,1) = SUM( H2_lev(1:NH2_lev)%density ) / neucompr / densHini
     DO k = 1, Nelements
        IF ( skipele(k) ) THEN
           CYCLE
        ELSE
           kk = kk + 1
           DO ll = 1, c_nvar
              c_am(kk,ll) = speci(c_ivar(ll))%formula(k)
           ENDDO
           somme = 0.0_dp
           DO i = 1, Nspec
              IF ( sp_var(i) ) CYCLE
              IF ( speci(i)%fluid == 1 ) THEN
                 somme = somme + DBLE(speci(i)%formula(k)) * speci(i)%density / neucompr
              ELSE
                 somme = somme + DBLE(speci(i)%formula(k)) * speci(i)%density / ioncompr
              ENDIF
           ENDDO
           somme = somme / densHini
           c_bbmn(kk,1) = elements(k)%ab_init - somme
        ENDIF
     ENDDO
     kk = kk + 1
     c_am(kk,kk)  = 1.0_dp
     IF ( charge_icon == ind_PAHm ) THEN
        c_am(kk,kk+1) = 1.0_dp
     ENDIF
     somme = 0.0_dp
     IF ( .NOT.sp_var(ind_PAH0) ) somme = somme + speci(ind_PAH0)%density / neucompr
     IF ( .NOT.sp_var(ind_PAHp) ) somme = somme + speci(ind_PAHp)%density / ioncompr
     IF ( .NOT.sp_var(ind_PAHm) ) somme = somme + speci(ind_PAHm)%density / ioncompr
     somme = somme / densHini
     c_bbmn(kk,1) = PAH_abinit - somme
     kk = kk + 1
     DO ll = 1, c_nvar
        c_am(kk,ll) = speci(c_ivar(ll))%charge
     ENDDO
     somme = 0.0_dp
     DO i = 1, Nspec
        IF ( sp_var(i) ) CYCLE
        somme = somme - speci(i)%charge * speci(i)%density / ioncompr
     ENDDO
     somme = somme / densHini
     c_bbmn(kk,1) = somme

     ! -------------------------------------------
     ! Write the system on screen - for debug
     ! -------------------------------------------
     ! DO kk = 1, c_nvar
     !    DO ll = 1, c_nvar
     !       WRITE(*,'(ES10.3,2X)',advance='no') c_am(kk,ll)
     !    ENDDO
     !    WRITE(*,'(5X,ES10.3)') c_bbmn(kk,1)
     ! ENDDO
     ! DO i = 1, Nspec
     !    IF (.NOT. sp_var(i) ) CYCLE
     !    IF ( speci(i)%fluid == 1 ) THEN
     !       WRITE(*,*) speci(i)%name, speci(i)%density / neucompr / densHini, sp_var(i)
     !    ELSE
     !       WRITE(*,*) speci(i)%name, speci(i)%density / ioncompr / densHini, sp_var(i)
     !    ENDIF
     ! ENDDO
     ! WRITE(*,*) 
     ! DO k = 1, Nelements
     !    WRITE(*,'(A8,2X,ES10.3)') elements(k)%Name, elements(k)%ab_init
     ! ENDDO
     ! WRITE(*,*)
     ! -------------------------------------------

     ! -------------------------------------------
     ! Solve system c_am * X = c_bbmn
     ! -------------------------------------------
     CALL DGESVX (c_FACT, c_TRANS, c_nvar, c_nrhs, c_am, c_nvar, c_AF, c_nvar, c_indx, &
                  c_EQUED, c_R, c_C, c_bbmn, c_nvar, c_bbmx, c_nvar, c_RCOND, c_FERR, &
                  c_BERR, c_WORK, c_IWORK, c_INFO)
     ! CALL DGESV (c_nvar, c_nrhs, c_am, c_nvar, c_indx, c_bbmn, c_nvar, c_INFO)

     ! WRITE(*,*) c_INFO

     ! DO kk = 1, c_nvar
     !    DO ll = 1, c_nvar
     !       WRITE(*,'(ES10.3,2X)',advance='no') c_am(kk,ll)
     !    ENDDO
     !    WRITE(*,'(5X,ES10.3,2X,ES10.3)') c_bbmn(kk,1)!, c_bbmx(kk,1)
     ! ENDDO
     ! READ(*,*)

     IF ( c_INFO == 0 ) THEN
        DO kk = 1, c_nvar
           i = c_ivar(kk)
           IF(c_bbmx(kk,1) <= 0.0_dp) THEN
              WRITE(*,*) "Error in forcing conservation equations - see dffun"
              WRITE(*,*) speci(i)%name, c_bbmx(kk,1)
              ! STOP
           ELSE IF ( speci(i)%fluid == 1 ) THEN
              speci(i)%density = c_bbmx(kk,1) * densHini * neucompr
           ELSE
              speci(i)%density = c_bbmx(kk,1) * densHini * ioncompr
           ENDIF
        ENDDO
     ENDIF

     ! -------------------------------------------
     ! Deallocate tables to solve the system
     ! -------------------------------------------
     DEALLOCATE ( sp_var  )
     DEALLOCATE ( c_ivar  )
     DEALLOCATE ( c_am    )
     DEALLOCATE ( c_bbmn  )
     DEALLOCATE ( c_bbmx  )
     DEALLOCATE ( c_indx  )
     DEALLOCATE ( c_AF    )
     DEALLOCATE ( c_R     )
     DEALLOCATE ( c_C     )
     DEALLOCATE ( c_FERR  )
     DEALLOCATE ( c_BERR  )
     DEALLOCATE ( c_WORK  )
     DEALLOCATE ( c_IWORK )

  END SUBROUTINE CONSISTENCY


  SUBROUTINE TEST_CONTINUE(again, message)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    test if the shock code must stop
     ! subroutine/function needed :
     ! input variables :
     !    again - boolean
     ! ouput variables :
     !    again - boolean : false if the shock code must stop
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_PHYS_VAR,          ONLY : Nstep_max, counter, stop_code, timeN,&
                                          coldens_h2, N_H2_0, coldens_co, N_CO_0,&
                                          shock_type, distance, Av, inv_Av_fac, Tn
     USE MODULE_VAR_VODE,          ONLY : duration_max, length_max

     IMPLICIT none

     LOGICAL,           INTENT(inout) :: again
     CHARACTER(LEN=80), INTENT(out)   :: message

     !----------------------------------------------------
     ! STOP integration if :
     !----------------------------------------------------

     ! (1) Maximal number of steps has been reached
     IF (counter == Nstep_max) THEN
        again = .FALSE.
        WRITE (message, '("Maximal number of steps has been reached.")' )
        stop_code = 1
     END IF

     ! (2) Maximal evolution time has been reached
     IF (timeN >= duration_max) THEN
        again = .FALSE.
        WRITE (message,'("Maximal evolution time has been reached. NH2=",1pe13.3," +",1pe13.3)') coldens_h2,N_H2_0
        stop_code = 2
     END IF

     ! (3) Maximal shock length has been reached
     IF (distance >= length_max) THEN
        again = .FALSE.
        PRINT *,'length_max',length_max
        PRINT *,'distance',distance
        WRITE (message,'("Maximal shock length has been reached. NH2=",1pe13.3,A1,1pe13.3)') coldens_h2,'+',N_H2_0
        stop_code = 3
     END IF

     ! (4) drift velocity < DeltaVmin
     ! IF (ABS_DeltaV < DeltaVmin * 1.0d-1) THEN
     ! ! IF ((Vn - Vi) < 10.0_DP) THEN
     !    again = .FALSE.
     !    message_end = "Equilibrium has been reached."
     !    stop_code = 4
     ! END IF

     ! (5) we reach a sonic point
     ! IF (Vn < (Vsound + 1.0D+1)) THEN
     !    again = .FALSE.
     !    message_end = "Sonic Point!"
     !    stop_code = 5
     ! END IF

     ! (6) Temperature is below 50.0 K
     ! IF (counter > 1000 .AND. Tn < 50.0_DP) THEN
     !    again = .FALSE.
     !    message_end = "temperature lower than 50 K"
     !    stop_code = 6
     ! END IF

     ! (6) Temperature decreases below initial T
     IF (old_Tn ==0 ) THEN ! detect first step
        Tinit  = Tn
        old_Tn = Tn
     ENDIF
     IF (Tn > old_Tn) THEN
        increase=.true.
     ENDIF
     IF (Tn < old_Tn .and. increase) THEN
        maximum_reached=.true.
        print *,'Reach maximum !'
        increase = .false.
     ENDIF
     IF ( (.false.) .and. (shock_type(1:1) /= 'S') .and. (shock_type(1:1) /= 'P') ) THEN
        again = .FALSE.
        WRITE (message, '("temperature decreases below 0.9*initial T.")')
        PRINT *, "T=", Tinit, 'Tn=', Tn
        stop_code = 6
     ELSE
        old_Tn = Tn
     ENDIF
 
     ! (7) Numerical instabilities prevent conservations of ions
     !!!!!!! IF (abs((DensityI-SUM(v_variab(bv_ion:ev_ion)))/DensityI) > 1.0e-2_DP) THEN
     !!!!!!!    again = .FALSE.
     !!!!!!!    message_end = "Ions are NOT conserved"
     !!!!!!!    stop_code = 7
     !!!!!!! ENDIF

     ! (8) Av is too small and we are integrating "backward"
     IF (Av<0_DP.and.inv_Av_fac<0_DP) THEN
        again = .false.
        WRITE (message, '("Av<0, that s a bit weird...")')
        stop_code = 8
     ENDIF

     ! (9) NH2 is too small as a result of backward integration
     IF ( (inv_Av_fac<0_DP ).and.(coldens_h2+N_H2_0<1e10_DP) ) THEN
        again = .false.
        PRINT *, 'NH2+N_H2_0=', coldens_h2+N_H2_0
        WRITE (message,'("NH2 residual < 1e10, we will soon reach NH2=0. Relative age=",1pe13.3)') timeN/duration_max
        stop_code = 9
     ENDIF

     ! (10) NCO is too small as a result of backward integration
     IF ( (inv_Av_fac<0_DP ).and.(coldens_co+N_CO_0<1e05_DP) ) THEN
        again = .false.
        PRINT *, 'NCO+N_CO_0=', coldens_co+N_CO_0
        WRITE (message,'("NCO residual < 1e05, we will soon reach NCO=0. Relative age=",1pe13.3)') timeN/duration_max
        stop_code = 10
     ENDIF

  END SUBROUTINE TEST_CONTINUE


  SUBROUTINE TEST_COUPLING(restart)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    test if ions and neutrals decouple or recouple or if a non-stationary
     !    CJ type shock needs to be run
     ! subroutine/function needed :
     ! input variables :
     ! ouput variables :
     !    restart - boolean : send code to restart the integration
     !                         with a J-type shock
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS
     USE MODULE_PHYS_VAR
     USE MODULE_RADIATION
     USE MODULE_VAR_VODE
     USE MODULE_H2,            ONLY : H2_lev
     USE MODULE_CO,            ONLY : CO_lev
     USE MODULE_EVOLUTION
     USE MODULE_ENERGETICS,    ONLY : mag_flux_corr
     USE MODULE_PROFIL_TABLES
     USE MODULE_SWITCH_CJ

     IMPLICIT none

     LOGICAL, INTENT(out) :: restart
     INTEGER              :: i

     restart = .false.

     !====================================================================================
     ! Test fluid decoupling strength and fluid recoupling hypothesis
     !====================================================================================
     IF ( shock_type == 'C' ) THEN
        IF     ( .NOT.fluid_decoup ) THEN
           CALL ESTIMATE_DECOUPLING( counter-1, eps_traj, decoup_strg )
           IF ( ABS(decoup_strg) / Vs_cm > eps_traj ) THEN
              fluid_decoup = .true.
              decoup_dist  = distance
           ENDIF
        ELSE IF ( ( ( distance > decoup_dist * (1.0_dp + eps_traj) ).AND.( .NOT.cj_active ).AND.( Vn > Vsound ) ).OR.&
                  ( ( distance > sonic_jump  * (1.0_dp + eps_traj) ).AND.( sonic_point    ) ) ) THEN
           CALL ESTIMATE_DECOUPLING( counter-1, eps_traj, decoup_strg )
           IF ( ABS(decoup_strg) / Vs_cm < eps_traj ) THEN
              fluid_recoup = .true.
           ENDIF
        ENDIF

        IF ( ( timeI >= timeJ .OR. fluid_recoup ).AND.&
             ( .NOT.cj_active .OR. sonic_point  ) ) THEN
           ! -------------------------------------
           ! set up a J-type shock
           ! -------------------------------------
           fluid_recoup = .true.

           ! prevent from switching back to C shocks
           shock_type_lock = "J"
           IF ( timeI >= timeJ ) THEN
              shock_type_fin = "CJ"
           ENDIF

           ! -------------------------------------
           ! reinitialize physical variables
           ! -------------------------------------
           CALL REINIT_PHYS_VARIABLES(traj_main(counter-1))
           dVn_old = dVn
           dVi_old = dVi
           IF (F_AV /= 0) THEN
              DO i = 1, nangle
                 irf(i,1:nwlg) = irf_old(i,1:nwlg)
              ENDDO
              ! optdpth(1:nwlg) = optdpth_old(1:nwlg)
           ENDIF
           DO i = 1, NH2_lev
              H2_lev(i)%cd_l = H2_lev(i)%cd_l_old
           ENDDO
           DO i = 1, NCO_lev
              CO_lev(i)%cd_l = CO_lev(i)%cd_l_old
           ENDDO

           ! -------------------------------------
           ! set neutral velocity to ion velocity
           ! assuming neutrals recouple with ions
           ! and not the opposite
           ! -------------------------------------
           ! CALL FALSE_JUMP(Vi)
           ! Vn = Vi
           ! v_variab(iv_Vn) = Vn
           ! v_lvariab(iv_Vn) = LOG(Vn)
           ! -------------------------------------
           ! set ion velocity to neutral velocity
           ! assuming ions recouple with neutrals
           ! and not the opposite
           ! -------------------------------------
           corr_flux_dum = 1.0_DP / Vi
           Vi = Vn
           v_variab(iv_Vi) = Vi
           v_lvariab(iv_Vi) = LOG(Vi)
           corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
           mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum

           ! -------------------------------------
           ! restart integration
           ! -------------------------------------
           counter = counter - 1
           restart = .true.
        ENDIF
     ENDIF

  END SUBROUTINE TEST_COUPLING


  SUBROUTINE TEST_STAT_CJ(restart)
     !---------------------------------------------------------------------------
     ! called by :
     !    MHD
     ! purpose :
     !    test if a stationary C* or CJ type shock needs to be run
     ! subroutine/function needed :
     ! input variables :
     ! ouput variables :
     !    restart - boolean : send code to restart the integration
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS
     USE MODULE_PHYS_VAR
     USE MODULE_RADIATION
     USE MODULE_VAR_VODE
     USE MODULE_EVOLUTION
     USE MODULE_ENERGETICS,    ONLY : mag_flux_corr
     USE MODULE_PROFIL_TABLES
     USE MODULE_SWITCH_CJ

     IMPLICIT none

     LOGICAL, INTENT(out) :: restart
     INTEGER              :: i

     restart = .false.

     !====================================================================================
     ! Switch to CJ-type or C*-type
     !   The procedure to reconstruct the trajectory of these shocks is used if
     !     - dVn changes sign or
     !     - dVi changes sign and we are already in a CJ or C* shock, or
     !     - Vn smaller than Vi and we are already in a CJ or C* shock, or
     !   Note : the procedure is not built for integrating Av
     !====================================================================================
     IF ( ( shock_type == 'C'                  ).AND.&
          ( .NOT.sonic_point                   ).AND.&
          ( ( Vn - Vsound ) / Vs_cm < eps_traj ).AND.&
          ( ( dVn*dVn_old < 0.0_dp ).OR.( dVi*dVi_old < 0.0_dp ) ) ) THEN

        IF ( F_AV == 1 ) THEN
           WRITE(*,*) "This procedure doesn't work if we integrate Av"
           WRITE(*,*) "-> code stops"
           STOP
        ENDIF

        ! -------------------------------------------
        ! Find minimum and maximum times to jump
        ! -------------------------------------------
        IF ( .NOT.cj_active .AND. .NOT.cs_active ) THEN
           ! ----------------------------------------
           ! save iteration number of original trajec
           ! ----------------------------------------
           counter_main = counter - 1
           ! ----------------------------------------
           ! Find the potential jump conditions
           ! ----------------------------------------
           CALL FIND_JUMP_RANGE(t_min_jump, t_max_jump, i_min_jump, i_max_jump, cs_active)
           cj_active = .true.
           ! ----------------------------------------
           ! Shortcut jump if C* type shock
           ! Useful to debug
           ! ----------------------------------------
           ! IF ( cs_active ) THEN
           !    i_min_jump = i_max_jump - 1
           ! ENDIF
        ELSE
           ! ----------------------------------------
           ! Dichotomie - converge on jumping time
           ! ----------------------------------------
           IF( dVn*dVn_old < 0.0_dp ) THEN
              t_max_jump = t_jump
              IF( i_max_jump - i_min_jump > 1 ) THEN
                 i_max_jump = i_jump
                 cs_active  = .false.
                 shock_type_fin  = "CJ"
              ENDIF
              ! -------------------------------------
              ! Save upper trajectory
              ! -------------------------------------
              traj_high(1:i_limit)           = traj_main(1:i_limit)
              traj_high(i_limit+1:counter-1) = traj_curr(i_limit+1:counter-1)
              found_high = .true.
              alpha = alpha - 1.0_dp / 2.0_dp**(nadj+1)
              counter_high = counter - 1
           ELSE
              t_min_jump = t_jump
              IF( i_max_jump - i_min_jump > 1 ) THEN
                 i_min_jump = i_jump
              ENDIF
              ! -------------------------------------
              ! Save lower trajectory
              ! -------------------------------------
              traj_down(1:i_limit)           = traj_main(1:i_limit)
              traj_down(i_limit+1:counter-1) = traj_curr(i_limit+1:counter-1)
              found_down = .true.
              alpha = alpha + 1.0_dp / 2.0_dp**(nadj+1)
              counter_down = counter - 1
           ENDIF
        ENDIF

        ! -------------------------------------------
        ! Test - recouple fluids as soon as possible
        ! if no possible jump conditions
        ! SHOULD BE REMOVED - NOT PHYSICAL
        ! -------------------------------------------
        IF (i_min_jump >= i_max_jump) THEN
           ! ----------------------------------------
           ! set up a J-type shock
           ! ----------------------------------------
           fluid_recoup = .true.

           ! prevent from switching back to C shock
           shock_type_lock = "J"
           shock_type_fin  = "CJ"

           !--- Find the jump index where ion and neutral recouple
           CALL FIND_RECOUPLING_JUMP(i_jump)
           ! i_jump = i_min_jump
           
           !--- reinitialize physical variables
           CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
           t_jump = distance
           dVn_old = dVn
           dVi_old = dVi
           i_limit = i_jump
           t_limit = t_jump
           !--- Apply Rankine-Hugoniot relations for an adiabatic jump
           !--- hydrodynamical shock - on neutral fluid only
           WRITE(*,*) "pre jump = ", Vn, Vi, Tn, Ti
           CALL RH_JUMP_HYDRO
           WRITE(*,*) "postjump = ", Vn, Vi, Tn, Ti
           ! READ(*,*) 
           ! ----------------------------------------
           ! set neutral velocity to ion velocity
           ! assuming neutrals recouple with ions
           ! and not the opposite
           ! ----------------------------------------
           ! Vn = Vi
           ! v_variab(iv_Vn) = Vn
           ! v_lvariab(iv_Vn) = LOG(Vn)
           ! -------------------------------------
           ! set ion velocity to neutral velocity
           ! assuming ions recouple with neutrals
           ! and not the opposite
           ! -------------------------------------
           corr_flux_dum = 1.0_DP / Vi
           Vi = Vn
           v_variab(iv_Vi) = Vi
           v_lvariab(iv_Vi) = LOG(Vi)
           corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
           mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
        ENDIF

        ! -------------------------------------------
        ! Reinitialize physical variables at i_jump
        ! or interpolate between i_jump and i_jump+1
        ! -------------------------------------------
        IF( .NOT.adjust_traj .AND. .NOT.adjust_start .AND. .NOT.fluid_recoup ) THEN
           IF ( i_max_jump - i_min_jump > 1 ) THEN
              i_jump = ( i_min_jump + i_max_jump ) / 2
              !--- reinitialize physical variables
              CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
              t_jump = distance
              dVn_old = dVn
              dVi_old = dVi
           ELSE
              t_jump = ( t_min_jump + t_max_jump ) / 2.0_dp
              i_jump  = i_min_jump
              CALL INTERP_JUMP_CONDITIONS(t_jump,i_jump)
           ENDIF
           i_limit = i_jump
           t_limit = t_jump
           !--- Apply Rankine-Hugoniot relations for an adiabatic jump
           !--- hydrodynamical shock - on neutral fluid only
           CALL RH_JUMP_HYDRO
           i_adjst = i_limit
        ENDIF

        ! -------------------------------------------
        ! Treatment of C* shocks
        ! -------------------------------------------
        IF ( ( ( ( i_max_jump - i_min_jump == 1 ).AND.( cs_active ) ).OR.( adjust_start ) ).AND.( .NOT.adjust_traj ) ) THEN
           ! ----------------------------------------
           ! scan backward the main trajectory and 
           ! select a point where to apply small 
           ! variations of Vn and Vi
           ! ----------------------------------------
           counter_main_old = counter_main
           IF ( ( .NOT. adjust_start ).OR.&
                ( adjust_start .AND. (.NOT.found_down .OR. .NOT.found_high) .AND. (alpha > 0.9_dp .OR. alpha < 0.1_dp) ) ) THEN
              IF( (counter_main - i_jump) / 5 /= 0 ) THEN
                 counter_main = counter_main - (counter_main - i_jump) / 5
              ELSE
                 counter_main = counter_main - 1
              ENDIF
           ENDIF
           i_limit = counter_main
           i_adjst = i_limit
           
           ! ----------------------------------------
           ! stop the code if we reach i_jump, the
           ! minimal index where the perturbation
           ! method can work
           ! ----------------------------------------
           IF (i_limit == i_jump) THEN
              WRITE(*,*) "Problem in starting C* trajectory"
              STOP
           ENDIF

           ! ----------------------------------------
           ! reinitialize all quantities at i_limit
           ! ----------------------------------------
           CALL REINIT_PHYS_VARIABLES(traj_main(i_limit))
           corr_flux_dum = 1.0_DP / Vi
           t_limit = distance
           dVn_old = dVn
           dVi_old = dVi

           ! ----------------------------------------
           ! Compute the maximum perturbations on 
           ! ion and neutral velocities
           ! ----------------------------------------
           IF ( counter_main_old /= counter_main ) THEN
              ! -------------------------------------
              ! reinitialize identification of up and
              ! down trajectories
              ! -------------------------------------
              found_high = .false.
              found_down = .false.
              ! -------------------------------------
              ! switch to an adjusted trajectory
              ! -------------------------------------
              adjust_start = .true.
              nadj =  0
              angle_traj = 0.0_dp
              CALL COMPUTE_ORTHO_TRAJ(i_limit,eps_traj,coef1,coef2)
              
              dum = 1.0_dp
              DO
                 IF( dum < 1e-7_dp ) THEN
                    CALL ROTATE_ORTHO_TRAJ(coef1,coef2, -angle_traj)
                    IF ( angle_traj < 0.0_dp ) THEN
                       angle_traj = - angle_traj
                    ELSE
                       angle_traj = - angle_traj - pi / 12.0_dp
                    ENDIF
                    IF ( ABS(angle_traj) > pi / 2.0_dp ) THEN
                       WRITE(*,*) "Problem in starting C* trajectory (trajectory rotation)"
                       STOP
                    ENDIF
                    CALL ROTATE_ORTHO_TRAJ(coef1,coef2, angle_traj)
                    dum = 1.0_dp
                 ENDIF
                 Vn_high = traj_main(i_limit)%Vn + dum * 1.5_dp * eps_traj * coef1
                 Vi_high = traj_main(i_limit)%Vi - dum * 1.5_dp * eps_traj * coef2
                 Vn = Vn_high
                 Vi = Vi_high
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 IF      (integ_type == 0) THEN
                    CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 ELSE IF (integ_type == 1) THEN
                    CALL diffun(d_v_var,LOG(t_limit),v_lvariab,v_dvariab)
                 ENDIF
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO
              DO
                 IF( dum < 1e-7_dp ) THEN
                    WRITE(*,*) "Problem in starting C* trajectory (variation amplitude)"
                    STOP
                 ENDIF
                 Vn_down = traj_main(i_limit)%Vn - dum * 1.5_dp * eps_traj * coef1
                 Vi_down = traj_main(i_limit)%Vi + dum * 1.5_dp * eps_traj * coef2
                 Vn = Vn_down
                 Vi = Vi_down
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 IF      (integ_type == 0) THEN
                    CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 ELSE IF (integ_type == 1) THEN
                    CALL diffun(d_v_var,LOG(t_limit),v_lvariab,v_dvariab)
                 ENDIF
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO

              Vn_high = traj_main(i_limit)%Vn + dum * 1.5_dp * eps_traj * coef1
              Vi_high = traj_main(i_limit)%Vi - dum * 1.5_dp * eps_traj * coef2
              Vn_down = traj_main(i_limit)%Vn - dum * 1.5_dp * eps_traj * coef1
              Vi_down = traj_main(i_limit)%Vi + dum * 1.5_dp * eps_traj * coef2
           ENDIF

           ! ----------------------------------------
           ! Apply perturbations on velocities
           ! ----------------------------------------
           IF (nadj == 0) THEN
              alpha = 1.0_dp / 2.0_dp**(nadj+1)
              nadj = nadj + 1
           ELSE
              nadj = nadj + 1
           ENDIF
           Vn = alpha * Vn_high + (1.0_dp - alpha) * Vn_down
           Vi = alpha * Vi_high + (1.0_dp - alpha) * Vi_down
           v_variab(iv_Vn) = Vn
           v_variab(iv_Vi) = Vi
           v_lvariab(iv_Vn) = LOG(Vn)
           v_lvariab(iv_Vi) = LOG(Vi)
           corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
           mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
        ENDIF

        ! -------------------------------------------
        ! Adjust trajectory after the jump to find
        ! the only solution : either the second sonic
        ! point or the downstream point
        ! -------------------------------------------
        IF ( ( ( i_max_jump - i_min_jump == 1 ).OR.( adjust_traj ) ).AND.( .NOT. sonic_point ) ) THEN
           i = i_limit + 1
           IF (found_high .AND. found_down) THEN
              DO
                 CALL FIND_DISTANCE(traj_high(i)%distance, Nstep_max, traj_down, counter_down, i_down)
                 WRITE(*,*) traj_high(i-1)%distance,      traj_high(i)%distance,      traj_high(i+1)%distance,&
                            traj_down(i_down-1)%distance, traj_down(i_down)%distance, traj_down(i_down+1)%distance
                 i_down = MAX(i_down, i_limit+1)
                 ! i_down = MAX(i_down, i_adjst+1)
                 dum = traj_down(i_down+1)%distance - traj_down(i_down)%distance
                 IF ( dum /= 0.0_dp .AND. i_down /= counter_down) THEN
                    dum1 = (traj_high(i)%distance        - traj_down(i_down)%distance) / dum &
                         * (traj_down(i_down+1)%Vn - traj_down(i_down)%Vn) + traj_down(i_down)%Vn
                    dum2 = (traj_high(i)%distance        - traj_down(i_down)%distance) / dum &
                         * (traj_down(i_down+1)%Vi - traj_down(i_down)%Vi) + traj_down(i_down)%Vi
                 ELSE
                    dum1 = traj_down(i_down)%Vn
                    dum2 = traj_down(i_down)%Vi
                 ENDIF
                 IF ( ( ABS( traj_high(i)%Vn - dum1 ) / traj_high(i)%Vn > eps_traj ).OR.&
                      ( ABS( traj_high(i)%Vi - dum2 ) / traj_high(i)%Vi > eps_traj ).OR.&
                      ( i == counter_high ) ) THEN
                    EXIT
                 ENDIF
                 i = i + 1
              ENDDO
           ENDIF
           i = i - 1

           IF ( i > i_limit ) THEN
              ! -------------------------------------
              ! switch to an adjusted trajectory
              ! -------------------------------------
              adjust_traj = .true.
              ! -------------------------------------
              ! average up and down trajectory and
              ! copy in main trajectory
              ! -------------------------------------
              CALL AVERAGE_TRAJEC ( traj_high, traj_down, traj_main, Nstep_max, counter_down, i_limit+1, i, i_adjst)
              counter_main = i
              nadj  = 0
           ENDIF

           IF ( adjust_traj .AND. (nadj == 0 .OR. nadj == 8) ) THEN
              ! -------------------------------------
              ! reinitialize identification of up and
              ! down trajectories
              ! reinitialize or change angle
              ! reinitialize alpha
              ! reinitialize nadj
              ! Careful - order is important
              ! -------------------------------------
              found_high = .false.
              found_down = .false.
              IF ( nadj == 0 ) THEN
                 angle_traj = 0.0_dp
              ENDIF
              IF ( nadj == 8 ) THEN
                 IF ( angle_traj < 0.0_dp ) THEN
                    angle_traj = - angle_traj
                 ELSE
                    angle_traj = - angle_traj - pi / 12.0_dp
                 ENDIF
              ENDIF
              IF ( ABS(angle_traj) > pi / 2.0_dp ) THEN
                 WRITE(*,*) "Problem in rotating trajectory"
                 STOP
              ENDIF
              nadj  = 0
              alpha = 1.0_dp / 2.0_dp**(nadj+1)

              ! -------------------------------------
              ! Compute the maximum perturbations on 
              ! ion and neutral velocities
              ! -------------------------------------
              CALL REINIT_PHYS_VARIABLES(traj_main(i))
              dVn_old = dVn
              dVi_old = dVi
              CALL COMPUTE_ORTHO_TRAJ(i,eps_traj,coef1,coef2)
              CALL ROTATE_ORTHO_TRAJ(coef1,coef2,angle_traj)

              dum = 1.0_dp
              DO
                 Vn_high = traj_main(i)%Vn + dum * 1.5_dp * eps_traj * coef1
                 Vi_high = traj_main(i)%Vi - dum * 1.5_dp * eps_traj * coef2

                 Vn = Vn_high
                 Vi = Vi_high
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 IF      (integ_type == 0) THEN
                    CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 ELSE IF (integ_type == 1) THEN
                    CALL diffun(d_v_var,LOG(t_limit),v_lvariab,v_dvariab)
                 ENDIF
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO
              DO
                 Vn_down = traj_main(i)%Vn - dum * 1.5_dp * eps_traj * coef1
                 Vi_down = traj_main(i)%Vi + dum * 1.5_dp * eps_traj * coef2

                 Vn = Vn_down
                 Vi = Vi_down
                 v_lvariab(iv_Vn) = LOG(Vn)
                 v_lvariab(iv_Vi) = LOG(Vi)
                 IF      (integ_type == 0) THEN
                    CALL diffun(d_v_var,t_limit,v_lvariab,v_dvariab)
                 ELSE IF (integ_type == 1) THEN
                    CALL diffun(d_v_var,LOG(t_limit),v_lvariab,v_dvariab)
                 ENDIF
                 dVn = v_dvariab(iv_Vn)
                 dVi = v_dvariab(iv_Vi)
                 IF ( ( dVn*dVn_old < 0.0_dp ).OR.( Vn > Vsound ).OR. ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
                    dum = dum / 2.0_dp
                 ELSE
                    EXIT
                 ENDIF
              ENDDO

              Vn_high = traj_main(i)%Vn + dum * 1.5_dp * eps_traj * coef1
              Vi_high = traj_main(i)%Vi - dum * 1.5_dp * eps_traj * coef2
              Vn_down = traj_main(i)%Vn - dum * 1.5_dp * eps_traj * coef1
              Vi_down = traj_main(i)%Vi + dum * 1.5_dp * eps_traj * coef2
           ENDIF

           IF ( adjust_traj ) THEN
              nadj = nadj + 1
              CALL REINIT_PHYS_VARIABLES(traj_main(i))
              corr_flux_dum = 1.0_DP / Vi
              dVn_old = dVn
              dVi_old = dVi
              t_limit = traj_main(i)%distance
              ! -------------------------------------
              ! slightly change ion and neutral
              ! velocities between Vn_high & Vn_down
              ! and between Vi_high & Vi_down
              ! -------------------------------------
              Vn = alpha * Vn_high + (1.0_dp - alpha) * Vn_down
              Vi = alpha * Vi_high + (1.0_dp - alpha) * Vi_down
              v_variab(iv_Vn) = Vn
              v_variab(iv_Vi) = Vi
              v_lvariab(iv_Vn) = LOG(Vn)
              v_lvariab(iv_Vi) = LOG(Vi)
              corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
              mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
           ENDIF
           i_limit = i
        ENDIF

        ! -------------------------------------------
        ! Criteria to end the trajectory - either
        ! crossing the sonic point or recoupling the
        ! fluid
        ! -------------------------------------------
        IF ( adjust_traj ) THEN
           IF ( ( ABS(traj_main(i_limit)%Vn - traj_main(i_limit)%Vsound) / Vs_cm < eps_traj ).AND.&
                ( traj_main(i_limit)%Vn < traj_main(i_limit)%Vsound ) ) THEN
              ! -------------------------------------
              ! test if we can cross the sonic point
              ! or not
              ! -------------------------------------
              CALL TEST_CROSS_SONIC(i_limit, eps_traj, test_sonic)
              IF ( test_sonic ) THEN
                 ! ----------------------------------
                 ! reinitialize variables
                 ! ----------------------------------
                 sonic_point = .true.
                 CALL REINIT_PHYS_VARIABLES(traj_main(i_limit))
                 dVn_old = dVn
                 dVi_old = dVi
                 ! ----------------------------------
                 ! extrapolate all variables at
                 ! d_extrapo = d + 2*(d_cross - d)
                 ! ----------------------------------
                 step_cross = 1.0_dp
                 CALL CROSS_SONIC_POINT(i_limit, eps_traj, t_limit, step_cross)
                 distance  = t_limit
                 dist_step = distance - traj_main(i_limit)%distance
                 sonic_jump = distance
                 ! ----------------------------------
                 ! recouple the fluids is necessary
                 ! DOES NOT WORK - recoupling just
                 ! after the jump is numerically
                 ! unstable - dont know why yet
                 ! ----------------------------------
                 ! IF( Vn < Vi .OR. ABS(Vn-Vi) / Vs_cm < eps_traj ) THEN
                 !    ! -------------------------------------
                 !    ! set neutral velocity to ion velocity
                 !    ! assuming neutrals recouple with ions
                 !    ! and not the opposite
                 !    ! -------------------------------------
                 !    Vn = Vi
                 !    v_variab(iv_Vn) = Vn
                 !    v_lvariab(iv_Vn) = LOG(Vn)
                 !    ! -------------------------------------
                 !    ! set up a J-type shock
                 !    ! -------------------------------------
                 !    ! prevent from switching back to C shocks
                 !    shock_type_lock = "J"
                 ! ENDIF
                 DO WHILE ( Vn < Vi )
                    step_cross = step_cross / 2.0_dp
                    CALL CROSS_SONIC_POINT(i_limit, eps_traj, t_limit, step_cross)
                    distance  = t_limit
                    dist_step = distance - traj_main(i_limit)%distance
                    sonic_jump = distance
                    IF( step_cross < eps_traj ) THEN
                       WRITE(*,*) "Problem in jumping sonic point"
                       STOP
                    ENDIF
                 ENDDO
              ENDIF
           ELSE IF ( ABS(traj_main(i_limit)%Vn - traj_main(i_limit)%Vi) / Vs_cm < eps_traj ) THEN
              ! -------------------------------------
              ! set up a J-type shock
              ! -------------------------------------
              fluid_recoup = .true.
              ! prevent from switching back to C shocks
              shock_type_lock = "J"
              ! -------------------------------------
              ! reinitialize variables
              ! -------------------------------------
              CALL REINIT_PHYS_VARIABLES(traj_main(i_limit))
              dVn_old = dVn
              dVi_old = dVi
              t_limit = traj_main(i_limit)%distance
              distance = t_limit
              ! -------------------------------------
              ! set neutral velocity to ion velocity
              ! assuming neutrals recouple with ions
              ! and not the opposite
              ! -------------------------------------
              ! CALL FALSE_JUMP(Vi)
              ! Vn = Vi
              ! v_variab(iv_Vn) = Vn
              ! v_lvariab(iv_Vn) = LOG(Vn)
              ! -------------------------------------
              ! set ion velocity to neutral velocity
              ! assuming ions recouple with neutrals
              ! and not the opposite
              ! -------------------------------------
              corr_flux_dum = 1.0_DP / Vi
              Vi = Vn
              v_variab(iv_Vi) = Vi
              v_lvariab(iv_Vi) = LOG(Vi)
              corr_flux_dum = corr_flux_dum - 1.0_DP / Vi
              mag_flux_corr = mag_flux_corr + Bfield**2.0_DP / (4._DP*pi) * Vs_cm**2._DP * corr_flux_dum
          ENDIF
        ENDIF

        ! ----------------------------------------
        ! restart integration
        ! ----------------------------------------
        counter = i_limit
        restart = .true.
     ENDIF

  END SUBROUTINE TEST_STAT_CJ


END MODULE MODULE_DRIVE
