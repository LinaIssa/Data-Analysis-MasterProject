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

PROGRAM MHD

  !#############################################################################
  !#############################################################################
  !##                                                                         ##
  !##                  Main program for the MHD shock model                   ##
  !##                  ------------------------------------                   ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !##                                                                         ##
  !## update :                                                                ##
  !## --------                                                                ##
  !##    * 2018 : Benjamin Godard                                             ##
  !##   - computation of C* and CJ stationary shocks                          ##
  !##                                                                         ##
  !##    * 2015-2017 : Benjamin Godard                                        ##
  !##   - change grain treatments - adsorption / erosion widths               ##
  !##   - add adsorption and desorption mechanisms                            ##
  !##   - update chemical network                                             ##
  !##   - add photoelectric effect for grain charge                           ##
  !##   - update photoreaction rates with those given by Heays et al. (2017)  ##
  !##   - update the computation of secondary photon processes                ##
  !##   - include radiation field spectrum and radiative transfer             ##
  !##   - include dust absorption properties and heat capacities              ##
  !##   - integrate photoreaction using cross sections                        ##
  !##   - compute the grain temperature through thermal balance               ##
  !##   - add radiative pumping of H2 electronic lines                        ##
  !##   - compute H2 self-shielding using the FGK approximation               ##
  !##   - change the way of computing H2 heating or cooling                   ##
  !##   - modification of dvode for absolute error control on log variables   ##
  !##   - possibility of sorting reactions when computing the derivatives     ##
  !##   - possibility of enforcing elements conserv for unbalanced network    ##
  !##   - possibility of reading H2 and CO dissociation rate in a PDR grid    ##
  !##   - remove Force_I_C parameter                                          ##
  !##   - add variables V_S and V_SH (individual species velocities)          ##
  !##   - introduce an energy criterium to compute the shock size             ##
  !##   - change input / output file - plug ISM services                      ##
  !##                                                                         ##
  !##    * fev 2010-2013 : Pierre Lesaffre                                    ##
  !##   - add variables NH2, NH for self-shielding                            ##
  !##   - add varibale v_CH to account for potential CH individual velocity   ##
  !##   - some debugging.                                                     ##
  !##   - revive photo-reactions.                                             ##
  !##   - Include self-shielding.                                             ##
  !##   - Interface with dumses.                                              ##
  !##   - python fitting routines: in fit/*py                                 ##
  !##                                                                         ##
  !##    * oct 2000 : David Wilgenbus                                         ##
  !##          translation into fortran 90 and revision of the mhdC code      ##
  !##          from G. Pineau des Forets and D. Flower                        ##
  !##                                                                         ##
  !##    * nov 2000 : David Wilgenbus                                         ##
  !##        - correction of grain compression                                ##
  !##        - H2O and CO cooling taken from Neufeld & Kaufman 1993           ##
  !##                                                                         ##
  !##    * mai/juin 2001 : Jacques Le Bourlot                                 ##
  !##        - Remplacement de GEAR par DVODE (d'apres Pierre Hily-Blant)     ##
  !##        - Revision generale (avec allegement)                            ##
  !##        - Correction de quelques bugs mineurs                            ##
  !##        - Suppression des variables 11 et 12                             ##
  !##          (Masse et rayon des grains)                                    ##
  !##                                                                         ##
  !##   * hiver/printemps 2001/2002 : DRF - GPdF - JLB                        ##
  !##        - Correction de nombreux bugs                                    ##
  !##        - Chocs J                                                        ##
  !##                                                                         ##
  !#############################################################################
  !#############################################################################

  USE MODULE_TECHCONFIG
  USE MODULE_TOOLS,             ONLY : GET_FILE_NUMBER, JACO
  USE MODULE_CONSTANTS,         ONLY : YEARsec, pi, Zero, parsec
  USE MODULE_INITIALIZE,        ONLY : READ_METADATA_NBL, READ_METADATA_DAT, INITIALIZE
  USE MODULE_CHEMICAL_SPECIES
  USE MODULE_CHEM_REACT,        ONLY : SORT_REACTIONS
  USE MODULE_H2,                ONLY : H2_lev, Index_VJ_H2, H2_lines
  USE MODULE_CO,                ONLY : CO_lev, Index_VJ_CO
  USE MODULE_PHYS_VAR
  USE MODULE_VAR_TH_BALANCE
  USE MODULE_DUST_TREATMENT,    ONLY : COMPUTE_QCOEF, DEALLOCATE_DUST
  USE MODULE_RADIATION
  USE MODULE_PROFIL_TABLES
  USE MODULE_VAR_VODE,          ONLY : integ_type, Eps_V, T0_V, &
                                       Tout_V, MF_V, Itask_V, Istate_V, Hdone, Hnext, &
                                       liw_V, lrw_V, itol_V, iopt_V, atol_V, rtol_V, &
                                       iwork_V, rwork_V, ipar_V, rpar_V,MAXORD_V,MXSTEP_V
  USE mvode
  USE MODULE_OUTPUTS
  USE MODULE_EVOLUTION
  USE MODULE_ENERGETICS,        ONLY : ENERGETIC_FLUXES
  USE MODULE_MOLECULAR_COOLING, ONLY : err_cool, f_err_cool, n_err_cool
  USE MODULE_DEBUG_JLB
  USE MODULE_LINE_EXCIT
  USE OUTPUT_HDF5,              ONLY : WRITE_ASCII
  USE MODULE_SWITCH_CJ
  USE MODULE_DRIVE
  ! USE NUM_REC,                  ONLY : gamminc0

  IMPLICIT NONE

  !=======================================================================================
  ! variables declaration
  !=======================================================================================

  INTEGER (KIND=LONG)                        :: i
  LOGICAL                                    :: its_OK  = .FALSE.  ! integration step is OK
  LOGICAL                                    :: again   = .TRUE.   ! continue integration while it is .TRUE.
  LOGICAL                                    :: restart = .false.  ! continue integration while it is .TRUE.
  CHARACTER (LEN=80)                         :: message_end        ! message written at the end of integration

  REAL (KIND=DP), DIMENSION (:), ALLOCATABLE :: dlv_var
  INTEGER                                    :: iflag

  !=======================================================================================
  ! 1 - Interface / configuration
  !=======================================================================================

   CALL CPU_TIME(time1)

  !-------------------------------------------------------
  ! message written on screen
  !-------------------------------------------------------
  WRITE(*,'(80("-"))')
  WRITE(*,'(20X,"MHD SHOCK")')
  WRITE(*,'(20X,"---------")')

  ! -------------------------------------------
  ! Choose the input file
  ! -------------------------------------------
  IF (iargc() .ne. 1) THEN
     inpfile = TRIM(data_dir) // TRIM("input_mhd.in")
  ELSE
     CALL getarg(1,inpfile)
  ENDIF

  !-------------------------------------------------------
  ! Test if the correct version of python exist
  !-------------------------------------------------------
  IF (F_W_HDF5_STD == 1 .OR. F_W_HDF5_CHE == 1) THEN
     sh_cmd = "which "//TRIM(python_version)
     CALL SYSTEM (TRIM(ADJUSTL(sh_cmd)), python_exists)
     IF (python_exists /= 0) THEN
        WRITE(*,*) '*** Error: python2.7 is not installed on the computer'
        WRITE(*,*) '           -> "which python2.7" = 0'
        WRITE(*,*) '           -> Shock code stops'
        STOP
     ENDIF
  ENDIF



  !=======================================================================================
  ! 2 - Preliminary treatments
  !=======================================================================================

  !-------------------------------------------------------
  ! initializations :
  !    * shock parameters
  !    * chemical species (+ check)
  !    * chemical reactions (+ check)
  !    * molecule H2 (levels, collision rates, lines)
  !    * molecule SiO (collision rates, Einstein coeff)
  !-------------------------------------------------------
  CALL INITIALIZE

  !-------------------------------------------------------
  ! Read the metadata of all output quantities
  ! 1 - read number of lines in metadata files
  ! 2 - allocate the metadata tables
  ! 3 - read the metadata files and store info
  !-------------------------------------------------------
  IF (F_W_HDF5_STD == 1) THEN
     fichier = TRIM(data_dir)//TRIM(metdat_dir)//TRIM(metstd_fil)
     CALL READ_METADATA_NBL(fichier, nhdf5_std)
     ALLOCATE( metatab_std(nhdf5_std) )
     CALL READ_METADATA_DAT(fichier, metatab_std, nhdf5_std)
  ENDIF
  IF (F_W_HDF5_CHE == 1) THEN
     fichier = TRIM(data_dir)//TRIM(metdat_dir)//TRIM(metche_fil)
     CALL READ_METADATA_NBL(fichier, nhdf5_che)
     ALLOCATE( metatab_che(nhdf5_che) )
     CALL READ_METADATA_DAT(fichier, metatab_che, nhdf5_che)
  ENDIF

  !-------------------------------------------------------
  ! HERE, adjust some quantities according to nH
  !-------------------------------------------------------
  XLL = XLL / nH

  !-------------------------------------------------------
  ! open the file for error messages
  !-------------------------------------------------------
  f_err_cool = GET_FILE_NUMBER()
  n_err_cool = TRIM(out_dir) // TRIM(n_err_cool)
  OPEN(f_err_cool, file=n_err_cool, STATUS='REPLACE', &
       access='SEQUENTIAL', form='FORMATTED', action='WRITE')

  !-------------------------------------------------------
  ! computes elemental abundances in (H, C, N, ...)
  ! finds the main carrier of each element
  !-------------------------------------------------------
  CALL ELEMENTAL_ABUNDANCES
  elements(1:Nelements)%ab_init = elements(1:Nelements)%ab
  ! DO i = 1, Nelements
  !    WRITE(*,'(A8, 2X, 2(ES10.3,2X))') elements(i)%name, elements(i)%ab_init, elements(i)%Dens_init
  ! ENDDO
  ! WRITE(*,'(A8,2X,ES10.3)') "PAH     ", PAH_abinit
  ! STOP

  CALL MAIN_CARRIER
  CALL INDEX_CONSERVATION

  !-------------------------------------------------------
  ! put into vectorial form
  !    * alloc of vectors containing physical variables
  !    * dimension = 1:d_v_var
  !    * deallocation is made at the end of MHD
  ! number of MHD variables
  !    * add NH, NH2, and NCO as independent variables
  !    * add CH, S, and SH velocities
  !-------------------------------------------------------
  Nv_MHD  = 22
  d_v_var = Nv_MHD + Nspec + NH2_lev_var

  ALLOCATE ( v_variab     (1:d_v_var) ) ! physical variables (Y)
  ALLOCATE ( v_lvariab    (1:d_v_var) ) ! logarithm of the variables
  ALLOCATE ( v_dz_lvariab (1:d_v_var) ) ! d/dz of the logarithm of the variables
  ALLOCATE ( v_dvariab    (1:d_v_var) ) ! value of Y*dY at the last call to DIFFUN
  ALLOCATE ( v_l_var      (1:d_v_var) ) ! Y at the last call to DIFFUN
  ALLOCATE ( v_l_der      (1:d_v_var) ) ! logarithm of the variables
  ALLOCATE ( dlv_var      (1:d_v_var) ) ! right hand side
  ALLOCATE ( YN           (0:Nspec)   ) ! change in number density (cm-3.s-1), calculated in CHEMISTRY
  ALLOCATE ( YNQ          (0:Nspec)   ) ! change in number density (cm-3.s-1), calculated in CHEMISTRY



  !=======================================================================================
  ! 3 - Initialization of all variables
  !     Initialization of dvode parameters
  !     Write info file
  !=======================================================================================

  !-------------------------------------------------------
  ! initialization
  !-------------------------------------------------------
  v_variab(1:d_v_var)     = Zero
  v_lvariab(1:d_v_var)    = Zero
  v_dz_lvariab(1:d_v_var) = Zero
  v_dvariab(1:d_v_var)    = Zero
  v_l_var(1:d_v_var)      = Zero
  v_l_der(1:d_v_var)      = Zero

  !-------------------------------------------------------
  ! MHD variables, number = Nv_MHD
  !-------------------------------------------------------
  i = 1  ; v_variab(i) = Vn         ; iv_Vn         = i ! neutrals (cm/s)
  i = 2  ; v_variab(i) = Vi         ; iv_Vi         = i ! ions (cm/s)
  i = 3  ; v_variab(i) = RhoN       ; iv_RhoN       = i ! neutrals (g.cm-3)
  i = 4  ; v_variab(i) = RhoI       ; iv_RhoI       = i ! ions (g.cm-3)
  i = 5  ; v_variab(i) = RhoA       ; iv_RhoA       = i ! anions (g.cm-3)
  i = 6  ; v_variab(i) = RhoNEG     ; iv_RhoNEG     = i ! negative ions and electrons (g.cm-3)
  i = 7  ; v_variab(i) = Tn         ; iv_Tn         = i ! neutrals (K)
  i = 8  ; v_variab(i) = Ti         ; iv_Ti         = i ! ions (K)
  i = 9  ; v_variab(i) = Te         ; iv_Te         = i ! electrons (K)
  i = 10 ; v_variab(i) = DensityN   ; iv_DensityN   = i ! neutrals (cm-3)
  i = 11 ; v_variab(i) = DensityI   ; iv_DensityI   = i ! ions (cm-3)
  i = 12 ; v_variab(i) = DensityA   ; iv_DensityA   = i ! anions (cm-3)
  i = 13 ; v_variab(i) = DensityNeg ; iv_DensityNeg = i ! negative ions and electrons (cm-3)
  i = 14 ; v_variab(i) = grad_V     ; iv_gv         = i ! neutral velocity gradient
  i = 15 ; v_variab(i) = coldens_h  ; iv_nh         = i ! H  column density (cm-2)
  i = 16 ; v_variab(i) = coldens_h2 ; iv_nh2        = i ! H2 column density (cm-2)
  i = 17 ; v_variab(i) = coldens_co ; iv_nco        = i ! CO column density (cm-2)
  i = 18 ; v_variab(i) = Vn         ; iv_vCH        = i ! v_CH velocity !CH!
  i = 19 ; v_variab(i) = Vn         ; iv_vS         = i ! v_S  velocity !S !
  i = 20 ; v_variab(i) = Vn         ; iv_vSH        = i ! v_SH velocity !SH!
  i = 21 ; v_variab(i) = compr_n    ; iv_compr_n    = i ! neutral compression factor
  i = 22 ; v_variab(i) = compr_i    ; iv_compr_i    = i ! ion     compression factor
  !!!!! i = 23 ; v_variab(i) = d_ero      ; iv_d_ero      = i ! erosion size
  !!!!! i = 24 ; v_variab(i) = d_ads      ; iv_d_ads      = i ! adsorption size

  !-------------------------------------------------------
  ! density of species (cm-3)
  ! number = Nspec (does not include added species)
  !-------------------------------------------------------
  bv_speci = Nv_MHD + 1
  bv_specy = bv_speci - 1
  ev_speci = bv_speci + Nspec - 1
  v_variab(bv_speci:ev_speci) = speci(1:Nspec)%density

  !-------------------------------------------------------
  ! index of neutrals, species on grain mantles,
  ! species on cores, positive ions, negative ions
  !-------------------------------------------------------
  bv_neu = b_neu + bv_speci - 1
  ev_neu = e_neu + bv_speci - 1
  bv_gra = b_gra + bv_speci - 1
  ev_gra = e_gra + bv_speci - 1
  bv_cor = b_cor + bv_speci - 1
  ev_cor = e_cor + bv_speci - 1
  bv_ion = b_ion + bv_speci - 1
  ev_ion = e_ion + bv_speci - 1
  bv_ani = b_ani + bv_speci - 1
  ev_ani = e_ani + bv_speci - 1
  bv_neg = b_neg + bv_speci - 1
  ev_neg = e_neg + bv_speci - 1

  !-------------------------------------------------------
  ! write the names of the species associated to grains
  ! debug purposes
  !-------------------------------------------------------
  ! WRITE(*,'(A,2(I3,2X))') 'bv_neu, ev_neu = ', bv_neu, ev_neu
  ! WRITE(*,'(A,2(I3,2X))') 'bv_gra, ev_gra = ', bv_gra, ev_gra
  ! WRITE(*,'(A,2(I3,2X))') 'bv_cor, ev_cor = ', bv_cor, ev_cor
  ! WRITE(*,'(A,2(I3,2X))') 'bv_ion, ev_ion = ', bv_ion, ev_ion
  ! WRITE(*,'(A,2(I3,2X))') 'bv_neg, ev_neg = ', bv_neg, ev_neg
  ! WRITE(*,*)
  ! DO i = bv_neu,ev_neu
  !    WRITE(*,*) i, bv_speci, i-bv_speci+1, speci(i-bv_speci+1)%name
  ! ENDDO
  ! WRITE(*,*)
  ! DO i = bv_gra,ev_gra
  !    WRITE(*,*) i, bv_speci, i-bv_speci+1, speci(i-bv_speci+1)%name
  ! ENDDO
  ! WRITE(*,*)
  ! DO i = bv_cor,ev_cor
  !    WRITE(*,*) i, bv_speci, i-bv_speci+1, speci(i-bv_speci+1)%name
  ! ENDDO
  ! WRITE(*,*)
  ! DO i = bv_ion,ev_ion
  !    WRITE(*,*) i, bv_speci, i-bv_speci+1, speci(i-bv_speci+1)%name
  ! ENDDO
  ! WRITE(*,*)
  ! DO i = bv_neg,ev_neg
  !    WRITE(*,*) i, bv_speci, i-bv_speci+1, speci(i-bv_speci+1)%name
  ! ENDDO
  ! WRITE(*,*)
  ! WRITE(*,*) ind_GRAIN, speci(ind_GRAIN)%name
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! density of H2 levels (cm-3)
  ! number = NH2_lev_var
  !-------------------------------------------------------
  bv_H2_lev = Nv_MHD + Nspec + 1
  ev_H2_lev = bv_H2_lev + NH2_lev_var - 1
  v_variab(bv_H2_lev:ev_H2_lev) = H2_lev(1:NH2_lev_var)%density

  !-------------------------------------------------------
  ! compute the log of v_variab
  !-------------------------------------------------------
  WHERE (v_variab > 0)
     v_lvariab = LOG(v_variab)
  ELSEWHERE
     v_lvariab = minus_infinity
  END WHERE

  !-------------------------------------------------------
  ! initial values of mass, momentum, and energy fluxes
  ! BG 2016 : initialize to 0 the source terms
  !-------------------------------------------------------
  Bn   = 0.0_DP
  Bi   = 0.0_DP
  Bneg = 0.0_DP
  CALL ENERGETIC_FLUXES

  !-------------------------------------------------------
  ! Initialization of dvode parameters / input
  !-------------------------------------------------------
  liw_V = d_v_var+30                     ! taille du tableau iwork_V
  lrw_V = 23 + 9*d_v_var + 2*d_v_var**2  ! taille du tableau rwork_V
  ALLOCATE ( iwork_V(liw_V)  )
  ALLOCATE ( rwork_V(lrw_V)  )
  ALLOCATE ( rtol_V(d_v_var) )
  ALLOCATE ( atol_V(d_v_var) )
  iwork_V = 0
  rwork_V = 0.0_DP
  rtol_V  = 0.0_DP
  atol_V  = 0.0_DP

  iwork_V(5) = MAXORD_V                  ! 2       max order of the integration
  iwork_V(6) = MXSTEP_V                  ! 1000000 max number of steps in one call to dvode

  !-------------------------------------------------------
  ! error control
  ! - 1 : purely relative => error = rtol_V * abs[ log(y(i)) ]
  ! - 2 : purely absolute => error = atol_V
  !-------------------------------------------------------
  ! 1 :
  ! ---
  ! atol_V  = tiny * 10_DP + 0.e-10_DP
  ! rtol_V  = Eps_V * 10.0_DP
  ! 2 :
  ! ---
  atol_V = Eps_V
  rtol_V = 0.0_DP

  !-------------------------------------------------------
  ! 1 - write informations about the current model 
  !     (parameters, species, chemistry, H2, grains)
  !      in the file file_info
  ! 2 - write small table on where are each elements on
  !     screen - not necessary, written in the info file
  !-------------------------------------------------------
  CALL WRITE_INFO
  ! CALL info_elements(6)



  !=======================================================================================
  ! 4 - Integration
  !     stops when all slabs have been run and again = .FALSE.
  !=======================================================================================

  WRITE(*,*) '--- start of the integration ---'

  DO slab = 1, nposgrd

     !-------------------------------------------------------
     ! Restart flag (e.g. run CJ Shocks)
     !-------------------------------------------------------
     999 CONTINUE
     CALL RESET_INTEGRATION(slab, again)

     DO WHILE ( again .and. t0_v < Tout_V )

        IF ( t0_v + Hnext < Tout_V ) THEN
           itask_v = 2
        ELSE
           itask_v = 1
           ! test
           ! itask_v = 4
           ! rwork_V(1) = Tout_V
           PRINT *, 'slab will terminate next !'
        ENDIF
        !----------------------------------------------------
        ! counts the number of integration steps
        !----------------------------------------------------
        counter = counter + 1

        !----------------------------------------------------
        ! Save the important radiative transfer quantities 
        ! from previous iteration
        !----------------------------------------------------
        IF (F_AV /= 0) THEN
           DO i = 1, nangle
              irf_old(i,1:nwlg) = irf(i,1:nwlg)
           ENDDO
           ! optdpth_old(1:nwlg) = optdpth(1:nwlg)
        ENDIF
        DO i = 1, NH2_lev
           H2_lev(i)%cd_l_old = H2_lev(i)%cd_l
        ENDDO
        DO i = 1, NCO_lev
           CO_lev(i)%cd_l_old = CO_lev(i)%cd_l
        ENDDO

        !----------------------------------------------------
        ! Sort the reaction rate from smallest to largest
        ! for each species
        !----------------------------------------------------
        CALL SORT_REACTIONS

        !----------------------------------------------------
        ! save some values before call to DRIVE
        !----------------------------------------------------
        IF      (integ_type == 0) THEN
           distance_old = T0_v
        ELSE IF (integ_type == 1) THEN
           distance_old = EXP(T0_v)
           ! WRITE(*,*) "T0    = ", T0_v
           ! WRITE(*,*) "d old = ", distance_old
           ! WRITE(*,*) "Tout  = ", Tout_V
           ! WRITE(*,*) "rwork_V(5) = ", rwork_V(5)
           ! WRITE(*,*) 
        ENDIF
        Tn_old       = Tn
        Vn_old       = Vn
        Vi_old       = Vi
        dVn_old      = dVn
        dVi_old      = dVi
        speci(1:Nspec)%Dens_old             = speci(1:Nspec)%density
        H2_lev(1:NH2_lev)%Dens_old          = H2_lev(1:NH2_lev)%density
        H2_lines(1:NH2_lines_out)%emiss_old = H2_lines(1:NH2_lines_out)%emiss
        CO_lev(1:NCO_lev)%Dens_old          = CO_lev(1:NCO_lev)%density

        !=================================================================================
        ! 4a - Numerical Integration - call DVODE, as modified by Pierre Hily-Blant
        !         most of variables are calculated in DIFFUN
        !         so they are available after call to DRIVE.
        !=================================================================================

        DO WHILE (.not. its_OK)

           !-------------------------------------------------
           ! Advance the solution from state v_lvariab
           ! at time T0_V to time Tout_V
           ! d_v_var is the number of ODE to solve
           !-------------------------------------------------
           CALL dvode( diffun, d_v_var, v_lvariab, T0_V, Tout_V, itol_V, &
                       rtol_V, atol_V, Itask_V, Istate_V, iopt_V, rwork_V, lrw_V, &
                       iwork_V, liw_V, jaco, MF_V, rpar_V, ipar_V )

           !-------------------------------------------------
           ! Raffraîchit les variables (PL)
           !-------------------------------------------------
           CALL diffun(d_v_var,T0_V,v_lvariab,v_dvariab)
           dVn = v_dvariab(iv_Vn)
           dVi = v_dvariab(iv_Vi)

           !-------------------------------------------------
           ! In MAIN_CARRIER and INDEX_CONSERVATION, the main
           ! carrier of every elements are calculated.
           !-------------------------------------------------
           CALL MAIN_CARRIER
           CALL INDEX_CONSERVATION

           !-------------------------------------------------
           ! test chemical conservation laws and apply
           ! corrections if needed
           !-------------------------------------------------
           CALL TEST_CONSERVATION_LAWS

           !-------------------------------------------------
           ! Info on dvode step
           !-------------------------------------------------
           IF (.false.) THEN
              PRINT *,'>','After DVODE:'
              PRINT *,'>','ISTATE=',Istate_v
              PRINT *,'>','Hdone       Hnext       Tcur       tolsf'
              PRINT "(10(1pe10.3,2x))",rwork_V(11:14)
              PRINT *,'>',' NST        NFE          NJE        NORD       NORD_NEXT     IMXER '
              PRINT "(10(i5,7x))",iwork_v(11:16)
              PRINT *,'>',' NLU        NNI          NCFN       NETF '
              PRINT "(10(i5,7x))",iwork_v(19:22)
           ENDIF

           !-------------------------------------------------
           ! Print info on length scales
           !-------------------------------------------------
           time_scale = 1.0_DP / MAXVAL( abs(v_dvariab) + tiny )
           i_time     = MAXLOC ( abs(v_dvariab),dim = 1 )
           ! PRINT *,' Length scale estimates:'
           ! PRINT "(' lenght_scale=',1pe13.3,' cm Max_var number i=',i4)",time_scale,i_time
           IF (i_time.ge.bv_H2_lev.and.i_time.le.ev_H2_lev) THEN
              PRINT "('Max_var number i=',i4,' H2 lev:',i4,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                    i_time,i_time-bv_H2_lev+1,exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
              ! PRINT *,' this corresponds to H2 level number ',i_time-bv_H2_lev+1
           ENDIF
           IF (i_time.ge.bv_speci.and.i_time.le.ev_speci) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,speci(i_time-bv_speci+1)%name,exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
              ! PRINT *,' this corresponds to species ', speci(i_time-bv_speci+1)%name
           ENDIF
           IF      ( i_time ==  1 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"Vn",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  2 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"Vi",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  3 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"RhoN",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  4 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"RhoI",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  5 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"RhoA",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  6 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"RhoNeg",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  7 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"Tn",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  8 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"Ti",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time ==  9 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"Te",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 10 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"DensityN",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 11 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"DensityI",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 12 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"DensityA",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 13 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"DensityNeg",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 14 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"gradV",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 15 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"NH",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 16 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"NH2",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 17 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"NCO",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 18 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"vCH",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 19 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"vS",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 20 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"vSH",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 21 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"compr_n",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ELSE IF ( i_time == 22 ) THEN
              PRINT "('Max_var number i=',i4,' speci: ',a,' y(i)=',1pe13.3,'  y(i)/dzy(i)=',1pe13.3,' cm')", &
                   i_time,"compr_i",exp(v_lvariab(i_time)),1d0/v_dvariab(i_time)
           ENDIF

           Hdone = rwork_V(11)
           Hnext = rwork_V(12)

           !-------------------------------------------------
           ! Success or failure of dvode
           !    - stop the code when necessary
           !    - output messages / information
           !-------------------------------------------------
           IF (Istate_V == 2) THEN
              its_OK = .TRUE.
              CALL DVINDY (T0_V, 1, rwork_V(21), d_v_var, dlv_var, iflag)
              dVn = dlv_var(1) * Vn
              dVi = dlv_var(2) * Vi
           ELSE IF (Istate_V == -1) THEN
              PRINT *, " "
              PRINT *, " Istate_V = ", Istate_V
              PRINT *, " MXSTEP = ", iwork_V(6)
              PRINT *, " "
              Tout_V = T0_V
              Hnext = Hnext * 0.1_DP
              iwork_V(6) = iwork_V(6) + 100
              Istate_V = 3
           ELSE IF (Istate_V == -2) THEN
              PRINT *, " Istate_V = ", Istate_V
              PRINT *, "  To be done"
              STOP
           ELSE IF (Istate_V == -3) THEN
              PRINT *, " Istate_V = ", Istate_V
              PRINT *, "  To be done"
              STOP
           ELSE IF (Istate_V == -4) THEN
              PRINT *, " "
              PRINT *, " Istate_V = ", Istate_V
              PRINT *, " IMXER = ", iwork_V(16), v_variab(iwork_V(16))
              PRINT *, " T = ", T0_V
              PRINT *, " Tout = ", Tout_V
              PRINT *, " Hnext = ", Hnext
              PRINT *, " "
              IF (counter /= countersave) THEN
                 count_istat4 = 1
                 countersave = counter
              ELSE
                 count_istat4 = count_istat4 + 1
              ENDIF
              IF (count_istat4 > 100) THEN
                 its_OK = .TRUE.
                 Tout_V = T0_V
                 WRITE(*,*) count_istat4
              ENDIF
              Hnext = Hnext * 0.1_DP
              Istate_V = 2
           ELSE IF (Istate_V == -5) THEN
              PRINT *, " "
              PRINT *, " Istate_V = ", Istate_V
              PRINT *, " IMXER = ", iwork_V(16), v_variab(iwork_V(16))
              PRINT *, " T = ", T0_V
              PRINT *, " Tout = ", Tout_V
              PRINT *, " Hnext = ", Hnext
              count_istat4 = count_istat4 + 1
              IF (count_istat4 > 100) THEN
                 its_OK = .TRUE.
                 Tout_V = T0_V
                 WRITE(*,*) count_istat4
              ENDIF
              Hnext = Hnext * 0.1_DP
              Istate_V = 2
           ELSE IF (Istate_V == -6) THEN
              PRINT *, " Istate_V = ", Istate_V
              PRINT *, "  To be done"
              STOP
           ELSE IF (Istate_V ==0) THEN
              PRINT *,'[solver = MEBDFI]'
           ELSE
              PRINT *, " Istate_V = ", Istate_V
              PRINT *, "  Unknown error"
              STOP
           ENDIF

        ENDDO

        its_OK = .FALSE.

        ! write (70,'(30ES20.11E3)') Tout_V, Hdone, Tn, Ti, Te, &
        !                            tr_1, tr_2, tr_3, tr_4, tr_5, &
        !                            tr_6, tr_7, tr_8, tr_9, tr_10, &
        !                            tr_11, tr_12, tr_13, tr_14, tr_15
        ! ------------------------------------------------------------

        !----------------------------------------------------
        ! calculate variables not already computed in DIFFUN
        !----------------------------------------------------

        !--- distance (cm) ---
        ! The following allows different origin of distance,
        ! hence to restart at T0_v = 0 at JC transition.
        IF      (integ_type == 0) THEN
           dist_step  = T0_v - distance_old
        ELSE IF (integ_type == 1) THEN
           dist_step  = EXP(T0_v) - distance_old
        ENDIF
        distance  = distance + dist_step

        !--- flow times for neutral and ions (year) ---
        timeN = timeN + 0.5_DP*dist_step*abs(1._DP/Vn_old+1._DP/Vn)/YEARsec
        timeI = timeI + 0.5_DP*dist_step*abs(1._DP/Vi_old+1._DP/Vi)/YEARsec

        !--- velocity gradient (km.s-1.cm-1) ---
        !--- lower limit : 1 km s-1 / pc
        !   Vgrad = ABS(v_l_der(iv_Vn)) / 1.D5  ! Computed in DIFFUN (evolution.f90)
        Vgrad = SQRT(v_l_der(iv_Vn)*v_l_der(iv_Vn) + 1.0D10 / (parsec*parsec)) * 1.0D-5

        !----------------------------------------------------
        ! calculate the column density (cm-2) of each
        ! chemical species and each H2 energy level
        ! note : trapezian rule
        !----------------------------------------------------
        ! a) chemical species
        speci(1:Nspec)%Col_dens = speci(1:Nspec)%Col_dens + 0.5_DP * dist_step &
                                * (speci(1:Nspec)%Dens_old + speci(1:Nspec)%density)
        ! b) H2 levels
        H2_lev(1:NH2_lev)%Col_dens = H2_lev(1:NH2_lev)%Col_dens + 0.5_DP * dist_step &
                                   * (H2_lev(1:NH2_lev)%Dens_old + H2_lev(1:NH2_lev)%density)

        ! b) CO levels
        CO_lev(1:NCO_lev)%Col_dens = CO_lev(1:NCO_lev)%Col_dens + 0.5_DP * dist_step &
                                   * (CO_lev(1:NCO_lev)%Dens_old + CO_lev(1:NCO_lev)%density)

        ! c) local Av
        Col_Dens_nH = speci(ind_H)%Col_Dens+ 2._DP * speci(ind_H2)%Col_Dens + speci(ind_Hplus)%Col_Dens

        !----------------------------------------------------
        ! Calculate the integrated intensity (erg/s/cm2/sr)
        ! of each H2 line and each fine structure line.
        !----------------------------------------------------
        ! a) H2 lines
        H2_lines(1:NH2_lines_out)%intensity = H2_lines(1:NH2_lines_out)%intensity + 0.5_DP * dist_step / (4._DP*pi)  &
                                            * (H2_lines(1:NH2_lines_out)%emiss_old + H2_lines(1:NH2_lines_out)%emiss)
        ! b) lines integrated intensities
        CALL LINE_INTEG

        !----------------------------------------------------
        ! calculate mass, momentum and energy fluxes
        ! and check the conservation of those quantities
        !----------------------------------------------------
        CALL ENERGETIC_FLUXES


        !=================================================================================
        ! 4c - Perform tests to stop, continue, or restart the shock code
        !=================================================================================

        !----------------------------------------------------
        ! Shall we stop the integration
        !----------------------------------------------------
        CALL TEST_CONTINUE(again,message_end)

        !----------------------------------------------------
        ! Test fluid decoupling strength and
        ! recoupling hypothesis
        ! Start a non-stationary CJ type shock if needed
        !----------------------------------------------------
        CALL TEST_COUPLING(restart)

        IF ( restart ) THEN
           GOTO 999
        ENDIF

        !----------------------------------------------------
        ! Test the existence of stationary CJ or C* shocks
        !----------------------------------------------------
        CALL TEST_STAT_CJ(restart)

        IF ( restart ) THEN
           GOTO 999
        ENDIF


        !=================================================================================
        ! 4d - Numerical Integration - end
        !      save shock profiles
        !      write output ascii files
        !=================================================================================

        !----------------------------------------------------
        ! save shock profil
        !----------------------------------------------------
        IF ( (.NOT.cj_active) .OR. (sonic_point) .OR. (fluid_recoup) ) THEN
           CALL FILL_PHYS_TABLES(traj_main(counter))
        ELSE
           CALL FILL_PHYS_TABLES(traj_curr(counter))
        ENDIF

        !----------------------------------------------------
        ! write output ascii files
        !----------------------------------------------------
        IF (F_W_ASCII == 1) THEN
           IF( counter == 1 ) THEN
              CALL WRITE_OUTPUT("first")
           ENDIF
           IF( cj_active.AND.(counter == i_limit+1) ) THEN
              CALL WRITE_OUTPUT("blank")
           ENDIF
           IF (MOD(counter-1,Nstep_w) == 0) THEN
              !--- write to several output file ---
              CALL WRITE_OUTPUT("usual")
              !--- write H2 excitation diagram
              CALL WRITE_EXCIT
           ENDIF
           ! IF(MOD(counter-1,Nstep_w*20) == 0) THEN
           !    CALL WRITE_RADIATION_FIELD("usual")
           ! ENDIF
           IF(.NOT.again) THEN
              CALL WRITE_OUTPUT("final")
           ENDIF
        ENDIF
        !--- write final abundances in input format
        CALL WRITE_SPECIES
        CALL WRITE_H2ABUND

        !--- write basic information on screen ---
        CALL CHEMISTRY
        CALL WRITE_SCREEN

        !----------------------------------------------------
        ! Switch off artificial viscosity after discontinuity
        ! in J shocks
        ! If (viscosity = .TRUE.), as initially, then the 
        ! derivative of the flow  velocity v is determined by
        ! integrating a second order equation for v, which
        ! arises from the artificial viscosity. This 
        ! integration proceeds through the J-shock 
        ! 'discontinuity', where v is varying rapidly, up to
        ! the point at which the velocity gradient becomes 
        ! equal to that calculated neglecting artificial 
        ! viscosity, to within a certain tolerance (see 
        ! below). itest is then set equal to 1, and the 
        ! integration is pursued neglecting the artificial
        ! viscosity terms (which have become negligible).
        ! Furthermore, the logarithmic derivatives are
        ! integrated beyond this point.
        !----------------------------------------------------
        IF ( ( shock_type == "J" )   .AND. &
             ( viscosity )           .AND. &
             ( ABS((save_dv + grad_V) / grad_V) < 2.0e-5_dp ) ) THEN
           viscosity      = .FALSE.
           viscosity_lock = .FALSE.
           PRINT *, 'Switch off viscosity'
        ENDIF

     ENDDO

  !-------------------------------------------------------
  ! end of integration
  !-------------------------------------------------------
  ENDDO

  !=======================================================================================
  ! 5 - Code ends / final treatments
  !     write messages on screen
  !     compute shock size, column densities, and line emissions
  !     write HDF5 outputs files
  !     deallocate tables
  !=======================================================================================

  !-------------------------------------------------------
  ! deal with shock profile
  ! - deallocate unused profils
  ! - reconstruct chemical rates along the main trajectory
  !-------------------------------------------------------
  CALL DEALLOCATE_TABLE(Nstep_max, traj_curr)
  DEALLOCATE(traj_curr)
  CALL DEALLOCATE_TABLE(Nstep_max, traj_down)
  DEALLOCATE(traj_down)
  CALL DEALLOCATE_TABLE(Nstep_max, traj_high)
  DEALLOCATE(traj_high)
  CALL FILL_CHEM_TABLE(Nstep_max, traj_main, counter)

  !-------------------------------------------------------
  ! compute shock integrated properties
  !    - shock size
  !    - column densities
  !    - line emissions
  !-------------------------------------------------------
  IF ( shock_type(1:1) /= 'S' .AND. &
       shock_type(1:1) /= 'P' .AND. &
       shock_type(1:1) /= 'W' ) THEN
     CALL COMPUTE_CHOC_SIZE(counter)
  ELSE IF ( shock_type(1:1) == 'P') THEN
     isize_shock = counter
     size_shock  = traj_main(isize_shock)%distance
     time_shock  = traj_main(isize_shock)%timeN
   ELSE
     isize_shock = 0
     size_shock  = 0.0_DP
  ENDIF

  IF (isize_shock == 0) THEN
     speci(1:Nspec)%Col_dens             = 0.0_DP
     H2_lev(1:NH2_lev)%Col_dens          = 0.0_DP
     H2_lines(1:NH2_lines_out)%intensity = 0.0_DP
     CO_lev(1:NCO_lev)%Col_dens          = 0.0_DP
     inthat (1:ntrhat)                   = 0.0_DP
     intcat (1:ntrcat)                   = 0.0_DP
     intoat (1:ntrnat)                   = 0.0_DP
     intnat (1:ntroat)                   = 0.0_DP
     intsat (1:ntrsat)                   = 0.0_DP
     intsiat(1:ntrsiat)                  = 0.0_DP
     intcpl (1:ntrcpl)                   = 0.0_DP
     intnpl (1:ntrnpl)                   = 0.0_DP
     intopl (1:ntropl)                   = 0.0_DP
     intspl (1:ntrspl)                   = 0.0_DP
     intsipl(1:ntrsipl)                  = 0.0_DP
  ELSE
     speci(1:Nspec)%Col_dens             = traj_main(isize_shock)%coldens(1:Nspec)
     H2_lev(1:NH2_lev)%Col_dens          = traj_main(isize_shock)%cold_h2lev(1:NH2_lev)
     H2_lines(1:NH2_lines_out)%intensity = traj_main(isize_shock)%inta_h2lin(1:NH2_lines_out)
     CO_lev(1:NCO_lev)%Col_dens          = traj_main(isize_shock)%cold_colev(1:NCO_lev)
     inthat (1:ntrhat)                   = traj_main(isize_shock)%inta_h  (1:ntrhat) 
     intcat (1:ntrcat)                   = traj_main(isize_shock)%inta_c  (1:ntrcat) 
     intoat (1:ntrnat)                   = traj_main(isize_shock)%inta_n  (1:ntrnat) 
     intnat (1:ntroat)                   = traj_main(isize_shock)%inta_o  (1:ntroat) 
     intsat (1:ntrsat)                   = traj_main(isize_shock)%inta_s  (1:ntrsat) 
     intsiat(1:ntrsiat)                  = traj_main(isize_shock)%inta_si (1:ntrsiat)
     intcpl (1:ntrcpl)                   = traj_main(isize_shock)%inta_cp (1:ntrcpl) 
     intnpl (1:ntrnpl)                   = traj_main(isize_shock)%inta_np (1:ntrnpl) 
     intopl (1:ntropl)                   = traj_main(isize_shock)%inta_op (1:ntropl) 
     intspl (1:ntrspl)                   = traj_main(isize_shock)%inta_sp (1:ntrspl) 
     intsipl(1:ntrsipl)                  = traj_main(isize_shock)%inta_sip(1:ntrsipl)
  ENDIF


  !-------------------------------------------------------
  ! write messages to screen and in info_mhd.out
  !-------------------------------------------------------
  WRITE(*,*) message_end
  WRITE(*,'(80("-"))')
  IF (err_cool) THEN
     WRITE(*,'("*** WARNING, check the file:", A)') n_err_cool
  ENDIF
  CALL WRITE_INFO_END

  !-------------------------------------------------------
  ! file closure
  !-------------------------------------------------------
  CLOSE(f_err_cool)


  !----------------------------------------------------
  ! Write standard hdf5 file
  !----------------------------------------------------
  IF (F_W_HDF5_STD == 1) THEN
     ! ------------------------------------------------
     ! Write standard ascii files in a dummy directory
     ! ------------------------------------------------
     ascii_dir_std = TRIM(out_dir)//TRIM(ascii_std)
     CALL WRITE_ASCII(ascii_dir_std, metatab_std, nhdf5_std,'s')
     ! ------------------------------------------------
     ! Call python scrpt to build HDF5 from ascii files
     ! ------------------------------------------------
     sh_cmd = TRIM(python_version)//" "//&
              TRIM("src/")//&
              TRIM(write_hdf5_file)//" "//&
              TRIM(ascii_dir_std)//" "//&
              TRIM(out_dir)//" "//&
              TRIM(ADJUSTL(modele))//"_s"//&
              " "//"1"
     PRINT *,"Launch command : ", TRIM(ADJUSTL(sh_cmd))
     CALL SYSTEM(sh_cmd)
     ! ------------------------------------------------
     ! Remove both the ascii files and the directory
     ! ------------------------------------------------
     ! sh_cmd = "rm -r "//TRIM(ascii_dir_std)
     ! CALL SYSTEM(sh_cmd)
  ENDIF
  ! ---------------------------------------------------
  ! Write the chemical hdf5 output files
  ! ---------------------------------------------------
  IF (F_W_HDF5_CHE == 1) THEN
     ! ------------------------------------------------
     ! Write chemical ascii files in a dummy directory
     ! ------------------------------------------------
     ascii_dir_che = TRIM(out_dir)//TRIM(ascii_che)
     CALL WRITE_ASCII(ascii_dir_che, metatab_che, nhdf5_che,'c')
     ! ------------------------------------------------
     ! Call python scrpt to build HDF5 from ascii files
     ! ------------------------------------------------
     sh_cmd = TRIM(python_version)//" "//&
              TRIM("src/")//&
              TRIM(write_hdf5_file)//" "//&
              TRIM(ascii_dir_che)//" "//&
              TRIM(out_dir)//" "//&
              TRIM(ADJUSTL(modele))//"_c"//&
              " "//"1"
     PRINT *,"Launch command : "
     PRINT *, sh_cmd
     CALL SYSTEM(sh_cmd)
     ! ------------------------------------------------
     ! Remove both the ascii files and the directory
     ! ------------------------------------------------
     ! sh_cmd = "rm -r "//TRIM(ascii_dir_che)
     ! CALL SYSTEM(sh_cmd)
  ENDIF


  ! ------------------------------------------------------
  ! Arrays deallocations
  ! ------------------------------------------------------

  !--- arrays allocated in MHD
  DEALLOCATE(v_variab)
  DEALLOCATE(v_lvariab)
  DEALLOCATE(v_dvariab)
  DEALLOCATE(v_l_var)
  DEALLOCATE(v_l_der)
  DEALLOCATE(YN)

  !--- arrays linked to metatable files
  IF (F_W_HDF5_STD == 1) DEALLOCATE(metatab_std)
  IF (F_W_HDF5_CHE == 1) DEALLOCATE(metatab_che)

  !--- arrays linked to dust properties
  CALL DEALLOCATE_DUST

  !--- arrays containing shock profiles
  CALL DEALLOCATE_TABLE(Nstep_max, traj_main)
  DEALLOCATE(traj_main)

  !--- arrays allocated in anywhere else
  DEALLOCATE(speci)
  DEALLOCATE(index_VJ_H2)
  DEALLOCATE(H2_lines)
  DEALLOCATE(index_VJ_CO)

  !--- deallocation of radiation variables
  DEALLOCATE (wlg       )
  DEALLOCATE (irf_free  )
  DEALLOCATE (irf_secpho)
  DEALLOCATE (irf       )
  DEALLOCATE (irf_old   )
  DEALLOCATE (urf       )
  DEALLOCATE (nrf       )

  !#############################################################################
  !#############################################################################
  !##                                                                         ##
  !##                 END OF THE MHD SHOCK CODE                               ##
  !##                                                                         ##
  !#############################################################################
  !#############################################################################

   CALL CPU_TIME(time2)
   WRITE(*,'(A10,ES10.3,A1)') 'step   1 : ', time2 - time1, 's'

END PROGRAM MHD
