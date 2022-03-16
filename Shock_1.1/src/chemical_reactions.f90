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

MODULE MODULE_CHEM_REACT
  !***************************************************************************************
  !**     The module 'MODULE_CHEM_REACT' contains variables, types, and subroutines     **
  !**     related to chemical reactions network: read datafiles, identify type of       **
  !**     reactions, check for identical reactions, compute endothermicities and        **
  !**     include reverse reactions                                                     **
  !***************************************************************************************

  USE MODULE_TECHCONFIG
  USE MODULE_CHEMICAL_SPECIES

  IMPLICIT NONE

  INCLUDE "precision.f90"

  !---------------------------------------------------------------------------------!
  ! input chemical files                                                            !
  !---------------------------------------------------------------------------------!
  CHARACTER(LEN=lenfilename)                       :: name_file_chemistry

  !---------------------------------------------------------------------------------!
  ! data type of one reaction                                                       !
  !---------------------------------------------------------------------------------!
  TYPE TYPE_REACTION
     INTEGER(KIND=LONG), DIMENSION(3)              :: R                             ! index of the 3 reactants ! ***Modif Tabone 12/2017 3BODY
     INTEGER(KIND=LONG), DIMENSION(4)              :: P                             ! index of the 4 products
     REAL(KIND=DP)                                 :: gamma                         ! Arrhenius coefficient
     REAL(KIND=DP)                                 :: alpha                         ! Arrhenius coefficient
     REAL(KIND=DP)                                 :: beta                          ! Arrhenius coefficient
     REAL(KIND=DP)                                 :: DE                            ! exo(if >0)/endo(if <0)-thermicity
     INTEGER(KIND=LONG),DIMENSION(3)               :: Nreac_3F                      ! number of reactant in each fluid n i e used to compute heat_exchange
     INTEGER(KIND=LONG),DIMENSION(3)               :: Nprod_3F                      ! number of products in each fluid n i e used to compute heat_exchange
     REAL(KIND=DP)                                 :: mass_prod                     ! sum of the product's masse (g)
     INTEGER(KIND=LONG)                            :: Nprod_m1                      ! number of products minus one
     CHARACTER(LEN=5)                              :: ref                           ! reference of the reaction
     CHARACTER(LEN=5)                              :: TYPE                          ! type of the reaction (string)
     INTEGER                                       :: itype                         ! type of the reaction (integer)
     REAL(KIND=DP)                                 :: rate                          ! reaction rate (in cm-3 s-1)
     INTEGER(KIND=LONG)                            :: idx_raw_integ                 ! index of raw integration (for photodestruction only so far)
     CHARACTER(LEN=14)                             :: dfile                         ! datafile if needed (for photodestruction so far)
  END TYPE TYPE_REACTION

  !---------------------------------------------------------------------------------!
  ! Type which contains the parameters of photodestruction reactions on which we    !
  ! base calculations of the photodestruction rates. The rates are computed at each ! 
  ! point by integrating the product of the radiation field with the photodestruct  !
  ! cross sections found listed in 'data/section/*.dat'                             !
  ! --------------------------------------------------------------------------------!
  TYPE, PUBLIC :: CROSS_SEC
     INTEGER                                       :: ispe                          ! index of the species
     INTEGER                                       :: npts                          ! number of points
     REAL(KIND=DP),      DIMENSION(:), ALLOCATABLE :: wl                            ! wavelength
     REAL(KIND=DP),      DIMENSION(:), ALLOCATABLE :: sigma                         ! Cross section in cm^2
  END TYPE CROSS_SEC

  !---------------------------------------------------------------------------------!
  ! Type for photoelectric reaction                                                 !
  ! --------------------------------------------------------------------------------!
  TYPE, PUBLIC :: PHELE_REACTION
     INTEGER                                       :: charge                        ! charge of the reactant
     REAL(KIND=DP)                                 :: rate                          ! reaction rate
  END TYPE PHELE_REACTION

  !---------------------------------------------------------------------------------!
  ! Number of reactions and type of reaction calculated in READ_REACTIONS           !
  ! Nreact is modified in ADD_INVERSE_REACTIONS as missing endothermic reactions    !
  !        are added (if the boolean do_we_add_reactions is .TRUE.)                 !
  !---------------------------------------------------------------------------------!
  LOGICAL,               PARAMETER                 :: do_we_add_reactions = .TRUE.  ! add the reverse reactions ?
  INTEGER(KIND=LONG),    PARAMETER                 :: Nreact_max = 2000             ! max. number of reactions
  INTEGER(KIND=LONG)                               :: Nreact                        ! number of reactions
  INTEGER                                          :: Ncross = 0                    ! number of photoreactions integrated over cross sections

  ! number of reactions in each type of reaction (PHOTO, SPUTTER, ...)
  ! calculated in INITIALIZE, exept Nrever_new et Nrever : in ADD_INVERSE_REACTIONS
  INTEGER(KIND=LONG)                               :: Nphoto = 0
  INTEGER(KIND=LONG)                               :: Nphele = 0
  INTEGER(KIND=LONG)                               :: Ngrrec = 0
  INTEGER(KIND=LONG)                               :: Ncr_io = 0
  INTEGER(KIND=LONG)                               :: Ncr_de = 0
  INTEGER(KIND=LONG)                               :: Nh2_fo = 0
  INTEGER(KIND=LONG)                               :: Nthree = 0
  INTEGER(KIND=LONG)                               :: N3body = 0                    ! ***Modif Tabone 12/2017 3BODY
  INTEGER(KIND=LONG)                               :: Nsputt = 0
  INTEGER(KIND=LONG)                               :: Nerosi = 0
  INTEGER(KIND=LONG)                               :: Nadsor = 0
  INTEGER(KIND=LONG)                               :: Nother = 0
  INTEGER(KIND=LONG)                               :: Nrever = 0
  INTEGER(KIND=LONG)                               :: Ndisso = 0
  INTEGER(KIND=LONG)                               :: Nphdes = 0
  INTEGER(KIND=LONG)                               :: Nsecde = 0
  INTEGER(KIND=LONG)                               :: Nthdes = 0

  ! index of beginning and end of each type of reaction
  ! calculated in INITIALIZE, exept b_rever et e_rever : in ADD_INVERSE_REACTIONS
  INTEGER(KIND=LONG)                               :: b_photo = 0, e_photo = 0
  INTEGER(KIND=LONG)                               :: b_phele = 0, e_phele = 0
  INTEGER(KIND=LONG)                               :: b_grrec = 0, e_grrec = 0
  INTEGER(KIND=LONG)                               :: b_cr_io = 0, e_cr_io = 0
  INTEGER(KIND=LONG)                               :: b_cr_de = 0, e_cr_de = 0
  INTEGER(KIND=LONG)                               :: b_h2_fo = 0, e_h2_fo = 0
  INTEGER(KIND=LONG)                               :: b_three = 0, e_three = 0
  INTEGER(KIND=LONG)                               :: b_3body = 0, e_3body = 0      ! ***Modif Tabone 12/2017 3BODY
  INTEGER(KIND=LONG)                               :: b_sputt = 0, e_sputt = 0
  INTEGER(KIND=LONG)                               :: b_erosi = 0, e_erosi = 0
  INTEGER(KIND=LONG)                               :: b_adsor = 0, e_adsor = 0
  INTEGER(KIND=LONG)                               :: b_phdes = 0, e_phdes = 0
  INTEGER(KIND=LONG)                               :: b_secde = 0, e_secde = 0
  INTEGER(KIND=LONG)                               :: b_thdes = 0, e_thdes = 0
  INTEGER(KIND=LONG)                               :: b_disso = 0, e_disso = 0
  INTEGER(KIND=LONG)                               :: b_other = 0, e_other = 0
  INTEGER(KIND=LONG)                               :: b_rever = 0, e_rever = 0

  !---------------------------------------------------------------------------------!
  ! vector of reactions                                                             !
  !---------------------------------------------------------------------------------!
  TYPE(TYPE_REACTION),  DIMENSION(Nreact_max)      :: react                         ! chemical         reactions (rates from input chemical files)
  TYPE(CROSS_SEC),      DIMENSION(:), ALLOCATABLE  :: phdest                        ! photodestruction reactions (with cross section integration)
  REAL(KIND=DP),        DIMENSION(:), ALLOCATABLE  :: phdest_rate                   ! photodestruction rates     (with cross section integration)
  REAL(KIND=DP),        DIMENSION(:), ALLOCATABLE  :: phheating_rate                ! Tabone photoheating
  TYPE(PHELE_REACTION), DIMENSION(:), ALLOCATABLE  :: phele_reac                    ! photoelectric    reactions (photon / grain interactions)

  !---------------------------------------------------------------------------------!
  ! dummy variables for determining wich reactions are correct -> READ_REACTIONS    !
  !---------------------------------------------------------------------------------!
  CHARACTER(LEN=name_length),PARAMETER             :: charact_blank = ''
  INTEGER(KIND=LONG),        PARAMETER             :: RP_undefined  = -1 ! if reactant or product undefined -> reaction is incorrect
  INTEGER(KIND=LONG),        PARAMETER             :: P_none        = 0  ! si product is not present in chemical species

  !---------------------------------------------------------------------------------!
  ! read/write format for chemical reactions                                        !
  !---------------------------------------------------------------------------------!
  CHARACTER(len=*), PRIVATE, PARAMETER  :: format_reaction = &
       '(A5,1X,5(A7,1X),A7,ES8.2,1X,F5.2,1X,F8.1,1X,I5,3X,A14,1X,F5.2)'

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE WRITE_REACTION(num,one_reaction)
    !--------------------------------------------------------------------------
    ! called by :
    !     * WRITE_INFO
    !     * CHECK_REACTIONS
    ! purpose :
    !     write informations about one chemical reaction
    ! subroutine/function needed :
    ! input variables :
    !     * num -> file number where to write the informations
    !     * one_reaction -> (type TYPE_REACTION) chemical reaction
    ! ouput variables :
    ! results :
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG),INTENT(in)::num
    TYPE(TYPE_REACTION),INTENT(in)::one_reaction

    IF (one_reaction%R(3)==0) THEN
       WRITE(num,format_reaction) &
            one_reaction%ref, &
            speci(one_reaction%R(1:2))%name, &
            speci(one_reaction%P(1:4))%name, &
            one_reaction%gamma, &
            one_reaction%alpha, &
            one_reaction%beta, &
            one_reaction%itype, &
            one_reaction%dfile, &
            one_reaction%DE
    ELSE ! Tabone 12/2017: 3BODY 3 reactants + 3 products
       WRITE(num,format_reaction)&
            one_reaction%ref, &
            speci(one_reaction%R(1:3))%name, &
            speci(one_reaction%P(1:3))%name, &
            one_reaction%gamma, &
            one_reaction%alpha, &
            one_reaction%beta, &
            one_reaction%itype, &
            one_reaction%dfile, &
            one_reaction%DE
    ENDIF

  END SUBROUTINE WRITE_REACTION



  FUNCTION REACTIONS_IDENTICAL(reaction1, reaction2) RESULT(res)
    !---------------------------------------------------------------------------
    ! called by :
    !     * CHECK_REACTIONS
    !     * ADD_INVERSE_REACTIONS
    ! purpose :
    !     test if 2 reactions are identical (same reactants and products)
    !     whatever the order in wich reactants and products are written
    ! subroutine/function needed :
    ! input variables :
    !     reaction1, reaction2 -> (type TYPE_REACTION) the 2 reactions
    ! output variables :
    ! results :
    !     (boolean) .TRUE. if the reactions are identical
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(TYPE_REACTION), INTENT(in) :: reaction1 , reaction2
    LOGICAL :: reactants_identical, produits_identical, res
    INTEGER(KIND=LONG) :: k,l
    LOGICAL,DIMENSION(4) :: prod_id1, prod_id2

    ! test if reactants are identical
    reactants_identical = &
         ( ( reaction1%R(1) == reaction2%R(1) .AND. &
             reaction1%R(2) == reaction2%R(2) ) .OR. &
           ( reaction1%R(1) == reaction2%R(2) .AND. &
             reaction1%R(2) == reaction2%R(1) ) )

    ! if reactants identical, test if products are identical
    produits_identical = .FALSE.
    IF (reactants_identical) THEN
       ! test if each product of reaction 1 is in reaction 2
       prod_id1(:) = .FALSE.
       DO k=1,4
          DO l=1,4
             IF (reaction1%P(k) == reaction2%P(l)) prod_id1(k) = .TRUE.
          END DO
       ENDDO
       ! test if each product of reaction 2 is in reaction 1
       prod_id2(:) = .FALSE.
       DO k=1,4
          DO l=1,4
             IF (reaction2%P(k)==reaction1%P(l)) prod_id2(k)=.TRUE.
          END DO
       ENDDO
       ! produits_identical is .TRUE. if all products are common to the 2 reactions
       produits_identical = ALL(prod_id1) .AND. ALL(prod_id2)
    END IF
    ! result is : (same reactants) and (sames products)
    res = reactants_identical .AND. produits_identical

  END FUNCTION REACTIONS_IDENTICAL



  SUBROUTINE READ_REACTIONS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read chemical reactions in the file file_chemistry
    !    and find index of each reactant and each product
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    react, Nreact
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS,     ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    USE MODULE_PHYS_VAR,  ONLY : shock_type_init
    IMPLICIT NONE

    CHARACTER(len=5), PARAMETER :: e_comment = '!end'
    CHARACTER(len=5), PARAMETER :: e_file    = 'END'
    CHARACTER(len=5)            :: ref = ''
    CHARACTER(len=8)            :: R1             ! name of 1st reactant
    CHARACTER(len=8)            :: R2             ! name of 2nd reactant
    CHARACTER(len=8)            :: R3             ! name of 3rd reactant
    CHARACTER(len=8)            :: P1             ! name of 1st product
    CHARACTER(len=8)            :: P2             ! name of 2nd product
    CHARACTER(len=8)            :: P3             ! name of 3rd product
    CHARACTER(len=8)            :: P4             ! name of 4th product
    REAL(KIND=DP)               :: alpha          ! reaction coefficient
    REAL(KIND=DP)               :: beta           ! reaction coefficient
    REAL(KIND=DP)               :: gamma          ! reaction coefficient
    INTEGER                     :: itype          ! type of reaction (integer code)
    CHARACTER(len=14)           :: dfile          ! name of datafile (if needed)
    INTEGER(KIND=LONG)          :: i, error
    LOGICAL                     :: skip1
    LOGICAL                     :: skip2
    LOGICAL                     :: skip3
    INTEGER                     :: file_chemistry ! unit to read chemistry file

    !---------------------------------------------
    ! 1 - Initialization of "react" vector
    !---------------------------------------------
    react(:)%type = ''
    react(:)%ref  = ''
    DO i=1,3 ! Tabone init heat_exch : loop on 3 fluids
       react(:)%Nreac_3F(i) = 0
       react(:)%Nprod_3F(i) = 0
    END DO  
    DO i=1,3
       react(:)%R(i) = RP_undefined
    END DO
    DO i=1,4
       react(:)%P(i) = RP_undefined
    END DO
    react(:)%gamma         = Zero
    react(:)%alpha         = Zero
    react(:)%beta          = Zero
    react(:)%DE            = Zero
    react(:)%Nprod_m1      = -1   ! don't forget the 'minus one' !
    react(:)%mass_prod     = Zero
    react(:)%idx_raw_integ = 0

    !---------------------------------------------
    ! 2 - Open chemical file and 
    !     skip all comments at the beginning
    !---------------------------------------------
    file_chemistry = GET_FILE_NUMBER()
    name_file_chemistry = TRIM(data_dir) // TRIM(chemfile)
    OPEN(file_chemistry,file=name_file_chemistry,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ref   = ''
    error = 0
    DO WHILE (ref /= e_comment .AND. error==0)
       READ(file_chemistry,'(A5)',iostat=error) ref
    ENDDO
    IF ( error /= 0 ) THEN
       STOP "*** WARNING, error in READ_REACTIONS"
    ENDIF

    !---------------------------------------------
    ! 3 - Read all reactions
    !     - stop at the end of file or if error
    !     - skip few reactions in special cases
    !     - read the reac/prod, rates, and types
    !---------------------------------------------
    Nreact = 0
    ref    = ''
    error  = 0
    DO WHILE ( (ref /= e_file) .AND. (Nreact < Nreact_max) .AND. (error==0) )
       Nreact = Nreact + 1

       !------------------------------------------
       ! Read one line of the chemical file
       !------------------------------------------
       R1 = '' ; R2 = '' ; R3 = '' ; P1 = '' ; P2 = '' ; P3 = '' ; P4 = '' ! Tabone 12/2017 3Body
       READ(file_chemistry,format_reaction,iostat=error) &
            ref, R1, R2, P1, P2, P3, P4, gamma, alpha, beta, itype, dfile

       IF ( itype == 10 ) THEN ! Tabone 12/2017 3Body: P1 means R3; P2 means P3.. in the 3BODY format     
          R3 = P1
          P1 = P2
          P2 = P3
          P3 = P4
          P4 = ''
       ENDIF

       !------------------------------------------
       ! Skip reactions involving grains when 
       ! computing a steady-state or a PDR
       ! Only skip erosion reactions - BG 2017
       !------------------------------------------
       skip1=.false.
       if (shock_type_init(1:1) == 'S') THEN
          !IF (R1(1:5)=='GRAIN' .OR. & 
          !    R2(1:5)=='GRAIN' .OR. & 
          !    P1(1:5)=='GRAIN' .OR. & 
          !    P2(1:5)=='GRAIN' .OR. & 
          !    P3(1:5)=='GRAIN' .OR. & 
          !    P4(1:5)=='GRAIN') THEN
          IF (ref(1:5)=='EROSI') THEN
             skip1=.true.
             print*,'Discard reaction: ',R1,R2,P1,P2,P3,P4
             Nreact=Nreact-1
          endif
       endif

       skip2=.false.
       IF (shock_type_init(1:1) == 'P') THEN
          !IF (R1(1:5)=='GRAIN' .OR. & 
          !    R2(1:5)=='GRAIN' .OR. & 
          !    P1(1:5)=='GRAIN' .OR. & 
          !    P2(1:5)=='GRAIN' .OR. & 
          !    P3(1:5)=='GRAIN' .OR. & 
          !    P4(1:5)=='GRAIN') THEN
          IF (ref(1:5)=='EROSI') THEN
             skip2=.true.
             PRINT *, 'Discard reaction: ',R1,R2,P1,P2,P3,P4
             Nreact=Nreact-1
          ENDIF
       ENDIF

       skip3=.false.
       IF (shock_type_init(1:1) == 'W') THEN
          !IF (R1(1:5)=='GRAIN' .OR. & 
          !    R2(1:5)=='GRAIN' .OR. & 
          !    P1(1:5)=='GRAIN' .OR. & 
          !    P2(1:5)=='GRAIN' .OR. & 
          !    P3(1:5)=='GRAIN' .OR. & 
          !    P4(1:5)=='GRAIN') THEN
          IF (ref(1:5)=='EROSI') THEN
             skip3=.true.
             PRINT *, 'Discard reaction: ',R1,R2,P1,P2,P3,P4
             Nreact=Nreact-1
          ENDIF
       ENDIF

       !------------------------------------------
       ! Fill all fields of the "react" vector
       !------------------------------------------
       IF (.not.skip1 .and. .not.skip2 .and. ref /= e_file .AND. (error==0)) THEN
          ! search for indexes of reactants and products
          DO i=0, Nspec_plus ! and not 1, Nspec, as we search also none, photon, crp, grain, SECPHO, VOISIN
             IF (R1 == speci(i)%name) react(Nreact)%R(1) = i
             IF (R2 == speci(i)%name) react(Nreact)%R(2) = i
             IF (R3 == speci(i)%name) react(Nreact)%R(3) = i ! Tabone 12/2017 3Body
             IF (P1 == speci(i)%name) react(Nreact)%P(1) = i
             IF (P2 == speci(i)%name) react(Nreact)%P(2) = i
             IF (P3 == speci(i)%name) react(Nreact)%P(3) = i
             IF (P4 == speci(i)%name) react(Nreact)%P(4) = i
          END DO

          !---------------------------------------
          ! compute the number of reactants & 
          ! products in each fluid
          ! Tabone 09-18 heat_exch
          !---------------------------------------
          DO i = 1, 3
             IF     (speci(react(Nreact)%R(i))%fluid == 1) THEN            ! i is neutral  reactant 
                react(Nreact)%Nreac_3F(1) = react(Nreact)%Nreac_3F(1) + 1
             ELSEIF (speci(react(Nreact)%R(i))%fluid == 2) THEN            ! i is positive reactant
                react(Nreact)%Nreac_3F(2) = react(Nreact)%Nreac_3F(2) + 1
             ELSEIF (speci(react(Nreact)%R(i))%fluid == 3) THEN            ! i is negative reactant
                react(Nreact)%Nreac_3F(3) = react(Nreact)%Nreac_3F(3) + 1
             ENDIF
          ENDDO
          DO i = 1, 4
             IF     (speci(react(Nreact)%P(i))%fluid == 1) THEN            ! i is neutral  reactant 
                react(Nreact)%Nprod_3F(1) = react(Nreact)%Nprod_3F(1) + 1
             ELSEIF (speci(react(Nreact)%P(i))%fluid == 2) THEN            ! i is positive reactant
                react(Nreact)%Nprod_3F(2) = react(Nreact)%Nprod_3F(2) + 1
             ELSEIF (speci(react(Nreact)%P(i))%fluid == 3) THEN            ! i is negative reactant
                react(Nreact)%Nprod_3F(3) = react(Nreact)%Nprod_3F(3) + 1
             ENDIF
          ENDDO

          !---------------------------------------
          ! fill the fields of the reaction
          !---------------------------------------
          react(Nreact)%ref   = ref
          react(Nreact)%gamma = gamma
          react(Nreact)%alpha = alpha
          react(Nreact)%beta  = beta
          react(Nreact)%itype = itype
          react(Nreact)%dfile = dfile

          !---------------------------------------
          ! calculates the number of products and 
          ! the sum of their mass (used in CHEMISTRY)
          !---------------------------------------
          IF (P1 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(1))%mass
          ENDIF
          IF (P2 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(2))%mass
          ENDIF
          IF (P3 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(3))%mass
          ENDIF
          IF (P4 /= charact_blank) THEN
             react(Nreact)%Nprod_m1 = react(Nreact)%Nprod_m1 + 1
             react(Nreact)%mass_prod = react(Nreact)%mass_prod + &
                  speci(react(Nreact)%P(4))%mass
          ENDIF

       END IF
    END DO
    ! if end of file has been reached, we counted 1 reaction too many
    IF (ref == e_file) Nreact = Nreact-1

    !---------------------------------------------
    ! 4 - file closure
    !---------------------------------------------
    CLOSE(file_chemistry)

  END SUBROUTINE READ_REACTIONS



  SUBROUTINE READ_PHOTODES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read ionization and dissociation cross section of 
    !    all the species listed in the photodest.flag file
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    react(j)%idx_raw_integ -> index k of the photodestruction cross section
    !                              to use for the reaction j
    !    phdest(k)%npts         -> number of points of cross sections for the
    !                              photodestruction k
    !    phdest(k)%wl           -> wavelength list of photodestruction k
    !    phdest(k)%sigma        -> photodest cross section of photodestruction k
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_TECHCONFIG, ONLY : fichier, data_dir, phddat_dir
    USE MODULE_TOOLS,      ONLY : GET_FILE_NUMBER

    IMPLICIT NONE
  
    CHARACTER(LEN=lenfilename), DIMENSION(:), ALLOCATABLE :: csect_files
    REAL(KIND=DP)                                         :: dum
    INTEGER                                               :: i, j, k
    INTEGER                                               :: nb_pt

    ! --------------------------------------------
    ! 1 - Get the names of cross section files. 
    !     These filenames are stored in the array 
    !     csect_file
    !     Identify the species involved and the
    !     index of the associate cross section
    ! --------------------------------------------
    ALLOCATE (csect_files(Ncross))
    ALLOCATE (phdest(Ncross))
    j = 0
    DO i=1, Nreact
       IF ( react(i)%type == 'CROSS' ) THEN
          j = j + 1
          csect_files(j) = TRIM(ADJUSTL(react(i)%dfile))
          phdest(j)%ispe = react(i)%R(1)
          react(i)%idx_raw_integ = j
       ENDIF
    ENDDO

    ! --------------------------------------------
    ! 2 - read the files containing the cross 
    !     sections and store the data in "phdest"
    !     Initialize of photodestruction rates
    ! --------------------------------------------
    DO j = 1, Ncross
       fichier = TRIM(data_dir)//TRIM(phddat_dir)//TRIM(csect_files(j))
       OPEN (iwrtmp, file = fichier, status = 'old')
       nb_pt = 0
       DO 
          READ (iwrtmp,*,iostat=k)
          IF ( k /= 0 ) EXIT
          nb_pt = nb_pt + 1 
       END DO
       nb_pt = nb_pt-10
       REWIND(iwrtmp)
       READ (iwrtmp,'(/////////)')

       phdest(j)%npts = nb_pt
       ALLOCATE(phdest(j)%wl(1:nb_pt))
       ALLOCATE(phdest(j)%sigma(1:nb_pt))
       phdest(j)%wl(1:nb_pt)    = 0.0_dp
       phdest(j)%sigma(1:nb_pt) = 0.0_dp
       DO i = 1, nb_pt
          READ(iwrtmp,*) phdest(j)%wl(i), dum, phdest(j)%sigma(i)
          phdest(j)%wl(i) = phdest(j)%wl(i) * 10.0_dp
       ENDDO
       CLOSE(iwrtmp)
    ENDDO

    DEALLOCATE(csect_files)

    ALLOCATE (phdest_rate(Ncross))
    ALLOCATE (phheating_rate(Ncross))
    phdest_rate(1:Ncross)    = 0.0_dp
    phheating_rate(1:Ncross) = 0.0_dp

    ! --------------------------------------------
    ! 3 - Check the index of raw data
    !     integration for each reaction
    !     Only for debug purposes
    ! --------------------------------------------
    ! DO j = 1, Nreact
    !    IF(react(j)%idx_raw_integ /= 0) THEN
    !       WRITE(*,*) speci(react(j)%R(1))%name,&
    !                  speci(react(j)%R(2))%name,&
    !                  speci(react(j)%P(1))%name,&
    !                  speci(react(j)%P(2))%name,&
    !                  speci(react(j)%P(3))%name,&
    !                  speci(react(j)%P(4))%name
    !       WRITE(*,*) j, react(j)%idx_raw_integ
    !    ENDIF
    ! ENDDO
    ! STOP
    ! ------------------------------------

  END SUBROUTINE READ_PHOTODES



  SUBROUTINE INITIALIZE_PHOTOELE
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    initialize the photoelectric structure and save the important
    !    feature (idx, charge of reactant) of the corresponding reaction
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT none

    INTEGER :: i, k

    ALLOCATE(phele_reac(Nphele))

    k = 0
    DO i = 1, Nreact
       IF(react(i)%type == "PHELE") THEN
          k = k + 1
          react(i)%idx_raw_integ = k
          phele_reac(k)%charge   = speci(react(i)%R(1))%charge
          phele_reac(k)%rate     = 0.0_dp
       ENDIF
    ENDDO
  END SUBROUTINE INITIALIZE_PHOTOELE



  SUBROUTINE ENERGY_DEFECT_REACTION
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    Calculate the energy defect (DE) for the reaction i.
    !    If this reaction has one specy with an unkown enthalpy, set DE = Zero
    !    exept if this specy appears (with the same occurence) in BOTH the
    !    reactants and products. In this case, the uncertainty cancels and one
    !    can calculate DE = enthalpy(reactants)-enthalpy(products).
    !    REMARK :
    !       energy defect is set to zero for endothermic reactions
    !       (DE < Zero) => DE = Zero (disabled by BG 2017)
    !    REMARK :
    !        energy defect for radiative recombinaison (PHOTON in the products)
    !        is set to Zero (Tabone 09/2018)
    ! subroutine/function needed :
    ! input variables :
    !    i -> index of the chemical reaction
    ! ouput variables :
    ! results :
    !     react
    ! Modif Tabone 12/2017 3body
    ! Modif Tabone 09/2018 loop on all reactions => allows to target clearl
    ! criteria for not to take DE in a reaction.
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : kCaleV, Zero

    IMPLICIT NONE

    REAL(KIND=DP),      PARAMETER               :: enthalpy_threshold=-99.9_DP
    LOGICAL                                     :: enthalpy_unknown
    INTEGER(KIND=LONG), DIMENSION(0:Nspec_plus) :: occ_reactants
    INTEGER(KIND=LONG), DIMENSION(0:Nspec_plus) :: occ_products
    INTEGER                                     :: j, ireac

    !-------------------------------------
    ! loop on reactions including REVERSE
    !-------------------------------------
    DO ireac = 1, Nreact
       IF ( react(ireac)%type == '3BODY' .OR. &
            react(ireac)%type == 'H2_FO' .OR. &
            react(ireac)%type == 'OTHER' .OR. &
            react(ireac)%type == 'REVER' .OR. &
            react(ireac)%type == 'EXCIT' ) THEN
          !-------------------------------
          ! check if enthalpy of one specy
          ! is unknown and if this specy 
          ! appears with the same occurence
          ! in reactants and in products
          !-------------------------------
          occ_reactants = 0 ! counts how many time appears one specy in the reactants
          occ_products  = 0 ! counts how many time appears one specy in the products
          !-------------------------------
          ! count occurence of reactants
          ! *** modif Tabone 3Body 12/2018
          !-------------------------------
          DO j=1,3
             occ_reactants(react(ireac)%R(j)) = occ_reactants(react(ireac)%R(j))+1
          END DO
          !-------------------------------
          ! count occurence of products
          ! *** modif Tabone 3Body 12/2018
          !-------------------------------
          DO j=1,4
             occ_products(react(ireac)%P(j)) = occ_products(react(ireac)%P(j))+1
          END DO
          enthalpy_unknown = .FALSE.
          !-------------------------------
          ! check enthalpy of reactants
          ! *** modif Tabone 3Body 12/2018
          !-------------------------------
          DO j=1,3
             IF (speci(react(ireac)%R(j))%enthalpy < enthalpy_threshold .AND. &
                  occ_reactants(react(ireac)%R(j)) /= occ_products(react(ireac)%R(j))) &
                  enthalpy_unknown = .TRUE.
          END DO
          !-------------------------------
          ! check enthalpy of products
          ! *** modif Tabone 3Body 12/2018
          !-------------------------------
          DO j=1,4
             IF (speci(react(ireac)%P(j))%enthalpy < enthalpy_threshold .AND. &
                  occ_reactants(react(ireac)%P(j)) /= occ_products(react(ireac)%P(j))) &
                  enthalpy_unknown = .TRUE.
          END DO

          !-------------------------------
          ! Tabone 09/18 DE = 0 for 
          ! radiative recombinaison
          !-------------------------------
          IF ( react(ireac)%P(1) == ind_PHOTON .OR.&
               react(ireac)%P(2) == ind_PHOTON .OR.&
               react(ireac)%P(3) == ind_PHOTON .OR.&
               react(ireac)%P(4) == ind_PHOTON ) THEN
             enthalpy_unknown = .TRUE.
          ENDIF

          !-------------------------------
          ! calculation of DE (eV, whereas
          ! speci%enthaply is in kCal/mol)
          !-------------------------------
          IF (enthalpy_unknown) THEN
             ! DE = 0.0 in this case
             react(ireac)%DE = Zero
          ELSE
             ! DE = enthalpy (reactants) - enthalpy (products)
             react(ireac)%DE = kCaleV &
                             * ( SUM(DBLE(speci(react(ireac)%R(:))%enthalpy)) &
                               - SUM(DBLE(speci(react(ireac)%P(:))%enthalpy)) )
          END IF
       ENDIF

       ! WRITE(*,'(9(A8),(F10.3),2(A8))') &
       !    speci(react(ireac)%R(:))%name, '->', speci(react(ireac)%P(:))%name, &
       !    'DE (eV)', react(ireac)%DE, 'TYPE:', react(ireac)%type
    ENDDO

  END SUBROUTINE ENERGY_DEFECT_REACTION



  SUBROUTINE REACTION_TYPE
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    Find the type of each reaction and re-order the entire set
    !    according to the different types. Compute the
    !    indexes (beginning, end, number) of each reaction type.
    !    Also calculate DE for each reaction, according to it's type.
    !    Order is :
    !        (0) standard chemical reactions     (type='OTHER')
    !        (X) reverse endothermic reactions   (type='REVER')
    !            => these reactions are defined and added in ADD_REVERSE_REACTIONS
    !        (1) H2 and HD formation             (type='H2_FO')
    !        (2) collisional dissociation of H2  (type='DISSO')
    !        (3) reactions with excited H2       (type='EXCIT')
    !        (4) radiative recombination         (type='OTHER')
    !       (10) gaseous 3-body reactions        (type='3BODY')
    !       (20) photo-reactions                 (type='PHOTO')
    !       (21) self-shielded photoreactions    (type='SHILD')
    !       (22) cross section photoreaction     (type='CROSS')
    !       (30) CR ionization or dissociation   (type='CR_IO')
    !       (40) photoelectric effect            (type='PHELE')
    !       (41) recombination on grains         (type='GRREC')
    !       (50) adsorption on grain surfaces    (type='ADSOR')
    !       (51) photodes from grain surfaces    (type='PHDES')
    !       (52) CR induced grain desorption     (type='CR_DE')
    !       (53) sec photon induced desorption   (type='SECDE')
    !       (54) Thermal desorpt. from surface   (type='THDES')
    !       (55) erosion of grain cores          (type='EROSI')
    !       (56) sputtering of grain mantles     (type='SPUTT')
    !       (57) three body reactions on surface (type='THREE')
    !    NOTE PL,5/8/2010 - some keywords ('excit', 'D&B', for example) 
    !         are now being transferred from the ref label to the type label.
    !         This allows to easily tag some reactions for which you wish a
    !         special treatment.
    !
    ! subroutine/function needed :
    !    ENERGY_DEFECT_REACTION
    ! input variables :
    ! ouput variables :
    ! results :
    !    reactions
    !
    !  22 juin 2001 - JLB - Add collisional dissociation of H2
    !       WARNING ! Required order in chemistry file : H2 + X -> X + H + H
    !                 identification is done on R1=H2, P2=H, P3=H, R2 = P1
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : Zero
    USE MODULE_PHYS_VAR,  ONLY : pdr_file_fnd

    IMPLICIT NONE

    INTEGER(KIND=LONG) :: i
    TYPE (TYPE_REACTION), DIMENSION(:), ALLOCATABLE :: react_aux

    !-------------------------------------
    ! Save the unsorted reaction set
    ! react_aux = old set (unsorted)
    !-------------------------------------
    ALLOCATE(react_aux(1:Nreact))
    react_aux(1:Nreact) = react(1:Nreact)

    !--- initialize the vector 'react' ---
    react(:)%type = ''
    react(:)%ref  = ''
    DO i=1,2
       react(:)%R(i) = RP_undefined
    END DO
    DO i=1,4
       react(:)%P(i) = RP_undefined
    END DO
    react(:)%gamma     = Zero
    react(:)%alpha     = Zero
    react(:)%beta      = Zero
    react(:)%DE        = Zero
    react(:)%Nprod_m1  = -1   ! useless here ..
    react(:)%mass_prod = Zero

    !===========================================================
    ! Find the type of each reaction.
    ! Re-order the entire set and calculate the energy defect
    ! according to this type
    ! NOW : 'react' = sorted reactions set
    !===========================================================

    !---------------------------------------------------
    ! (0) gaseous 3-body reactions  (type='3BODY')
    !     *** Modif Tabone 12/2017 3Body
    !---------------------------------------------------
    b_3body = 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 10 ) THEN
          N3body = N3body + 1
          react(N3body+ b_3body-1)      = react_aux(i)
          react_aux(i)%type             = '3BODY'
          react(N3body+ b_3body-1)%type = '3BODY'
       ENDIF
    END DO
    e_3body = b_3body + N3body - 1 ! index of end

    !---------------------------------------------------
    ! (20) photoreactions               (type='PHOTO')
    ! (21) self-shielded photoreactions (type='SHILD')
    ! (22) cross section photoreaction  (type='CROSS')
    !---------------------------------------------------
    b_photo = e_3body + 1
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 20 .OR.&
            react_aux(i)%itype == 21 .OR.&
            react_aux(i)%itype == 22 ) THEN
          Nphoto = Nphoto + 1
          !---------------------------------------------
          ! Treatment of all photoreactions
          !---------------------------------------------
          react(Nphoto + b_photo-1)     = react_aux(i)
          !---------------------------------------------
          ! Treatment of reactions with self-shielding
          !---------------------------------------------
          IF     ( react_aux(i)%itype == 21 ) THEN
             react_aux(i)%type             = 'SHILD'
             react(Nphoto+ b_photo-1)%type = 'SHILD'
          ELSE IF( react_aux(i)%itype == 22 ) THEN
             react_aux(i)%type             = 'CROSS'
             react(Nphoto+ b_photo-1)%type = 'CROSS'
             Ncross = Ncross + 1
          ELSE
             react_aux(i)%type             = 'PHOTO'
             react(Nphoto+ b_photo-1)%type = 'PHOTO'
          ENDIF
          react(Nphoto + b_photo-1)%ref = react_aux(i)%ref
          IF( (react_aux(i)%ref == 'PDR16').AND.(.NOT.pdr_file_fnd) ) THEN
             WRITE(*,*) "You can't use the H2 and CO photodissociation rates computed"
             WRITE(*,*) "with the Meudon PDR Code -> input files are missing         "
             WRITE(*,*) "   ---> code stops"
             STOP
          ENDIF
       ENDIF
    END DO
    e_photo = b_photo + Nphoto - 1 ! index of end

    !---------------------------------------------------
    ! (40) photoelectric effect (type='PHELE')
    !---------------------------------------------------
    b_phele = e_photo + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 40 ) THEN
          Nphele = Nphele + 1
          react(Nphele + b_phele-1)      = react_aux(i)
          react_aux(i)%type              = 'PHELE'
          react(Nphele + b_phele-1)%type = 'PHELE'
       ENDIF
    END DO
    e_phele = b_phele + Nphele - 1 ! index of end

    !---------------------------------------------------
    ! (41) recombination on grains (type='GRREC')
    !---------------------------------------------------
    b_grrec = e_phele + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 41 ) THEN
          Ngrrec = Ngrrec + 1
          react(Ngrrec + b_grrec-1)      = react_aux(i)
          react_aux(i)%type              = 'GRREC'
          react(Ngrrec + b_grrec-1)%type = 'GRREC'
       ENDIF
    END DO
    e_grrec = b_grrec + Ngrrec - 1 ! index of end

    !---------------------------------------------------
    ! (30) CR ionization or dissociation (type='CR_IO')
    !---------------------------------------------------
    b_cr_io = e_grrec + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 30 ) THEN
          Ncr_io = Ncr_io + 1
          react(Ncr_io+ b_cr_io-1)      = react_aux(i)
          react_aux(i)%type             = 'CR_IO'
          react(Ncr_io+ b_cr_io-1)%type = 'CR_IO'
          ! set to beta a high value when beta is small
          IF (react(Ncr_io+ b_cr_io-1)%beta < 1.D-9) THEN
             react(Ncr_io+ b_cr_io-1)%beta = 1.D8
          ENDIF
       ENDIF
    END DO
    e_cr_io = b_cr_io + Ncr_io - 1 ! index of end

    !---------------------------------------------------
    ! (52) CR induced grain desorption (type='CR_DE')
    !---------------------------------------------------
    b_cr_de = e_cr_io + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 52 ) THEN
          Ncr_de = Ncr_de + 1
          react(Ncr_de+ b_cr_de-1)      = react_aux(i)
          react_aux(i)%type             = 'CR_DE'
          react(Ncr_de+ b_cr_de-1)%type = 'CR_DE'
       ENDIF
    END DO
    e_cr_de = b_cr_de + Ncr_de - 1 ! index of end

    !---------------------------------------------------
    ! (1) H2 and HD formation (type='H2_FO')
    !---------------------------------------------------
    b_h2_fo = e_cr_de + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 1 ) THEN
          Nh2_fo = Nh2_fo + 1
          react(Nh2_fo+ b_h2_fo-1)      = react_aux(i)
          react_aux(i)%type             = 'H2_FO'
          react(Nh2_fo+ b_h2_fo-1)%type = 'H2_FO'
       ENDIF
    END DO
    e_h2_fo = b_h2_fo + Nh2_fo - 1 ! index of end

    !---------------------------------------------------
    ! (57) three body reactions on surface (type='THREE')
    !---------------------------------------------------
    b_three = e_h2_fo + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 57 ) THEN
          Nthree = Nthree + 1
          react(Nthree+ b_three-1)      = react_aux(i)
          react_aux(i)%type             = 'THREE'
          react(Nthree+ b_three-1)%type = 'THREE'
       ENDIF
    END DO
    e_three = b_three + Nthree - 1 ! index of end

    !---------------------------------------------------
    ! (56) sputtering of grain mantles (type='SPUTT')
    !---------------------------------------------------
    b_sputt = e_three + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 56 ) THEN
          Nsputt = Nsputt + 1
          react(Nsputt+ b_sputt-1)      = react_aux(i)
          react_aux(i)%type             = 'SPUTT'
          react(Nsputt+ b_sputt-1)%type = 'SPUTT'
       ENDIF
    END DO
    e_sputt = b_sputt + Nsputt - 1 ! index of end

   !---------------------------------------------------
    ! (55) erosion of grain cores (type='EROSI')
    !---------------------------------------------------
    b_erosi =  e_sputt + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 55 ) THEN
          Nerosi = Nerosi + 1
          react(Nerosi+ b_erosi-1)      = react_aux(i)
          react_aux(i)%type             = 'EROSI'
          react(Nerosi+ b_erosi-1)%type = 'EROSI'
       ENDIF
    END DO
    e_erosi = b_erosi + Nerosi - 1 ! index of end

    !---------------------------------------------------
    ! (50) adsorption on grain surfaces (type='ADSOR')
    !---------------------------------------------------
    b_adsor = e_erosi + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 50 ) THEN
          Nadsor = Nadsor + 1
          react(Nadsor+ b_adsor-1)      = react_aux(i)
          react_aux(i)%type             = 'ADSOR'
          react(Nadsor+ b_adsor-1)%type = 'ADSOR'
       ENDIF
    END DO
    e_adsor = b_adsor + Nadsor - 1 ! index of end

    !---------------------------------------------------
    ! (51) photodes from grain surfaces (type='PHDES')
    !---------------------------------------------------
    b_phdes = e_adsor + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 51 ) THEN
          Nphdes = Nphdes + 1
          react(Nphdes+ b_phdes-1)      = react_aux(i)
          react_aux(i)%type             = 'PHDES'
          react(Nphdes+ b_phdes-1)%type = 'PHDES'
       ENDIF
    END DO
    e_phdes = b_phdes + Nphdes - 1 ! index of end

    !--------------------------------------------------
    ! (53) sec photon induced desorption (type='SECDE')
    !--------------------------------------------------
    b_secde = e_phdes + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 53 ) THEN
          Nsecde = Nsecde + 1
          react(Nsecde+ b_secde-1)      = react_aux(i)
          react_aux(i)%type             = 'SECDE'
          react(Nsecde+ b_secde-1)%type = 'SECDE'
       ENDIF
    END DO
    e_secde = b_secde + Nsecde - 1 ! index of end

    !---------------------------------------------------
    ! (54) Thermal desorpt. from surface (type='THDES')
    !---------------------------------------------------
    b_thdes = e_secde + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 54 ) THEN
          Nthdes = Nthdes + 1
          react(Nthdes+ b_thdes-1)      = react_aux(i)
          react_aux(i)%type             = 'THDES'
          react(Nthdes+ b_thdes-1)%type = 'THDES'
       ENDIF
    END DO
    e_thdes = b_thdes + Nthdes - 1 ! index of end

    !---------------------------------------------------
    ! (2) collisional dissociation of H2 (type='DISSO')
    !---------------------------------------------------
    b_disso = e_thdes + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 2 ) THEN
          Ndisso = Ndisso + 1
          react(Ndisso+ b_disso-1)      = react_aux(i)
          react_aux(i)%type             = 'DISSO'
          react(Ndisso+ b_disso-1)%type = 'DISSO'
       ENDIF
    END DO
    e_disso = b_disso + Ndisso - 1 ! index of end

    !---------------------------------------------------
    ! (0) standard chemical reactions (type='OTHER')
    ! (3) reactions with excited H2   (type='EXCIT')
    ! (4) radiative recombination     (type='OTHER')
    !---------------------------------------------------
    b_other = e_disso + 1 ! index of beginning
    DO i=1, Nreact
       IF ( react_aux(i)%itype == 0 .OR. &
            react_aux(i)%itype == 3 .OR. &
            react_aux(i)%itype == 4 ) THEN
          Nother = Nother + 1
          !---------------------------------------------
          ! Differentiate reactions with excited H2 and
          ! all standard chemical reactions
          !---------------------------------------------
          IF ( react_aux(i)%itype == 3 ) THEN
             react(Nother+ b_other-1)      = react_aux(i)
             react_aux(i)%type             = 'EXCIT'
             react(Nother+ b_other-1)%type = 'EXCIT'
          ELSE
             react(Nother+ b_other-1)      = react_aux(i)
             react_aux(i)%type             = 'OTHER'
             react(Nother+ b_other-1)%type = 'OTHER'
          ENDIF
       ENDIF
    END DO
    e_other = b_other + Nother - 1 ! index of end

    !--- Nrever, b_rever, e_rever -> see ADD_REVERSE_REACTIONS ---

    !--- deallocation of the temporary variable ---
    DEALLOCATE(react_aux)

  END SUBROUTINE REACTION_TYPE



  SUBROUTINE CHECK_REACTIONS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    tests if the set of reactions is correct
    !        * each reaction must have at least 2 reactants and 1 product
    !        * reactants and products must be in the set of chemical species
    !        * species have to be conserved (exept for reaction of the type
    !           THREE, EROSI or ADSOR)
    !        * charge has to be conserved
    ! subroutine/function needed :
    !    WRITE_REACTION
    !    REACTIONS_IDENTICAL
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : screen
    
    IMPLICIT NONE
    INTEGER(KIND=LONG)                         :: i,j,ii
    INTEGER(KIND=LONG), DIMENSION(3,Nelements) :: formula_reactants ! Tabone 12/2017 3Body
    INTEGER(KIND=LONG), DIMENSION(4,Nelements) :: formula_products
    INTEGER(KIND=LONG)                         :: charge_reactants
    INTEGER(KIND=LONG)                         :: charge_products
    CHARACTER(LEN=8)                           :: name
    LOGICAL                                    :: incorrect

    DO i=1,Nreact
       ! --- checks if reaction is correct ---
       incorrect = &
            ! at least 2 reactants
       react(i)%R(1) == P_none .OR. &
       react(i)%R(2) == P_none .OR. &
            ! at least 1 product
       react(i)%Nprod_m1 < 0 .OR. &
            ! use only species in the current set
       react(i)%R(1) == RP_undefined .OR. &
       react(i)%R(2) == RP_undefined .OR. &
       react(i)%R(3) == RP_undefined .OR. &  ! modif Tabone 12/2017 3Body
       react(i)%P(1) == RP_undefined .OR. &
       react(i)%P(2) == RP_undefined .OR. &
       react(i)%P(3) == RP_undefined .OR. &
       react(i)%P(4) == RP_undefined
       IF (incorrect) THEN
          WRITE(*,*) "*** WARNING, reaction ",i," is incorrect"
          CALL WRITE_REACTION(screen,react(i))
          STOP
       END IF

       ! --- species conservation ---
       formula_reactants(:,:) = 0
       formula_products(:,:)  = 0
       SELECT CASE (react(i)%type)
       CASE ('THREE', 'EROSI', 'ADSOR')
          ! for these reaction, one can't chek the conservation of species
       CASE DEFAULT
          DO j=1,3 ! chemical formula of reactants
             name=speci(react(i)%R(j))%name
             SELECT CASE (name)
             CASE ('PHOTON', 'CRP', 'ELECTR', 'GRAIN','SECPHO', 'VOISIN')
                ! do nothing
             CASE DEFAULT
                formula_reactants(j,:) = speci(react(i)%R(j))%formula
             END SELECT
          END DO
          DO j=1,4 ! chemical formula of products
             name=speci(react(i)%P(j))%name
             SELECT CASE (name)
             CASE ('PHOTON', 'CRP', 'ELECTR', 'GRAIN')
                ! do nothing
             CASE DEFAULT
                formula_products(j,:) = speci(react(i)%P(j))%formula
             END SELECT
          END DO
          ! checks the conservation of species between reactants and products
          DO j=1,Nelements
             IF (SUM(formula_reactants(:,j)) /= SUM(formula_products(:,j))) THEN
                WRITE(*,*)"*** WARNING, element ",elements(j)%name," is not conserved in the reaction ",i
                CALL WRITE_REACTION(screen,react(i))
                print*,'reactants',formula_reactants(:,j)
                print*,'products',formula_products(:,j)
                do ii=1,nspec
                   print*,speci(ii)%name,speci(ii)%formula(j)
                enddo
                STOP
             END IF
          ENDDO
       END SELECT

       ! --- conservation of the charge ---
       charge_reactants = SUM(speci(react(i)%R(1:2))%charge)
       charge_products  = SUM(speci(react(i)%P(1:4))%charge)
       IF (charge_reactants /= charge_products) THEN
          WRITE(*,*)"*** WARNING, charge isn't conserved in the reaction ",i
          CALL WRITE_REACTION(screen,react(i))
          STOP
       ENDIF

       ! --- look for the same reaction in the rest of the chemical set ---
       DO j=i+1,Nreact
          IF (REACTIONS_IDENTICAL(react(i),react(j))) THEN
             WRITE(*,*)"*** WARNING, reactions ",i," and ",j," are identical"
             CALL WRITE_REACTION(screen,react(i))
             CALL WRITE_REACTION(screen,react(j))
!!!PIERRE!!!             STOP
          END IF
       END DO
    END DO
  END SUBROUTINE CHECK_REACTIONS



  SUBROUTINE ADD_REVERSE_REACTIONS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    add reverse reactions when :
    !        * ion-neutral reaction (no barrier)
    !        * 2 reactants <-> 2 products
    !        * 0< DE <= DE_threshold
    !        * no radiative asspciation
    !        * no recombination
    !        * reverse reaction isn't already included
    !    the reverse reaction has the following :
    !        * same gamma, alpha as direct reaction
    !        * beta=DE(direct)*11600.
    !        * DE=0.0
    ! subroutine/function needed :
    !    REACTIONS_IDENTICAL
    ! input variables :
    ! ouput variables :
    ! results :
    !    react, Nreact, Nrever ...
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE

    REAL(KIND=DP), PARAMETER :: DE_threshold=2._DP ! eV
    INTEGER(KIND=LONG) :: i,j
    LOGICAL :: eliminate
    TYPE(TYPE_REACTION) :: reaction_aux

    ! the reactions are added at the end of the set
    b_rever=Nreact+1

    ! looks for ion-neutral reaction only -> type 'OTHER'
    DO i=b_other,e_other
       ! eliminates the following reactions
       eliminate= &
            (react(i)%Nprod_m1 /= 1)     .OR. &    ! more or less than 2 products
            (react(i)%DE <= Zero)        .OR. &    ! DE <= 0
            (react(i)%DE > DE_threshold) .OR. &    ! DE > DE_threshold
            ((react(i)%R(1)>=b_neu) .AND. (react(i)%R(1)<=e_neu) .AND.&
             (react(i)%R(2)>=b_neu) .AND. (react(i)%R(2)<=e_neu)).OR. & ! neutral-neutral (possibility of a barrier)
            (react(i)%R(1)==ind_e      .OR. react(i)%R(2)==ind_e).OR. & ! recombination
            (react(i)%P(1)==ind_PHOTON .OR. react(i)%P(2)==ind_PHOTON)  ! radiative association

       ! tests if the reverse reaction isn't already in the set
       j=1
       DO WHILE ((j<=e_other) .AND. (.NOT. eliminate))
          ! build the reverse reaction
          reaction_aux%R(1)=react(i)%P(1)
          reaction_aux%R(2)=react(i)%P(2)
          reaction_aux%R(3)= 0              ! Tabone 12/2017 no ion/neutral 3body reactions!!! => to be updated if ion/neutral 3body reac included
          reaction_aux%P(1)=react(i)%R(1)
          reaction_aux%P(2)=react(i)%R(2)
          reaction_aux%P(3:4)=P_none
          reaction_aux%Nreac_3F(1:3) = react(i)%Nprod_3F(1:3) ! Tabone heat_exchange 09-18 react become prod
          reaction_aux%Nprod_3F(1:3) = react(i)%Nreac_3F(1:3) ! Tabone heat_exchange 09-18 prod  become react
          ! test if it doesn't exist already
          eliminate=REACTIONS_IDENTICAL(react(j),reaction_aux)
          j=j+1
       END DO

       ! if the reaction is not eliminated, we add the corresponding reverse reaction
       IF (.NOT. eliminate) THEN
          ! update of Nreact
          Nreact=Nreact+1
          ! Nreact_max must be large enough
          IF (Nreact > Nreact_max) STOP "*** WARNING, Nreact_max is too small"
          ! update of Nrever (e_rever is calculated at the end of the subroutine)
          Nrever=Nrever+1

          ! definition of the reverse reaction
          react(Nreact)=reaction_aux ! reactants and products
          react(Nreact)%gamma=react(i)%gamma
          react(Nreact)%alpha=react(i)%alpha
          react(Nreact)%beta=react(i)%DE*11600._DP
          react(Nreact)%DE=Zero
          react(Nreact)%Nprod_m1=1 ! only 2 products
          react(Nreact)%mass_prod= &
               speci(react(Nreact)%P(1))%mass + &
               speci(react(Nreact)%P(2))%mass
          react(Nreact)%itype=0
          react(Nreact)%type='REVER'
          react(Nreact)%ref='REVER'
          react(Nreact)%idx_raw_integ = 0
       END IF
    END DO
    ! calculate e_rever
    e_rever=b_rever+Nrever-1

  END SUBROUTINE ADD_REVERSE_REACTIONS


  SUBROUTINE SORT_REACTIONS
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD_VODE
    ! purpose :
    !    Sort the reaction rates by increasing order of absolute values for 
    !    each species
    !    If first call, allocate the matrix containing the indices of the
    !    reactions for each species. Else, fill this matrix
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    speci(:)%nbreac
    !    speci(:)%sort_reac
    !---------------------------------------------------------------------------
    IMPLICIT none

    REAL(KIND=dp), DIMENSION(Nreact)  :: ratef
    REAL(KIND=dp), DIMENSION(Nreact)  :: rated
    INTEGER,       DIMENSION(Nreact)  :: indxf
    INTEGER,       DIMENSION(Nreact)  :: indxd
    INTEGER                           :: IR1, IR2, IR3, IP1, IP2, IP3, IP4  ! modif Tabone 12/2017 3Body
    INTEGER                           :: i, j, nbd, nbf
    LOGICAL, SAVE                     :: first = .true.
    

    !----------------------------------------------------
    ! Initial call - find the number of reactions
    !----------------------------------------------------
    IF (first) THEN
       speci(1:nspec)%nbreac_f = 0 
       speci(1:nspec)%nbreac_d = 0 
       DO i = 1, Nreact
          IR1 = react(i)%R(1) ; IR2 = react(i)%R(2)
          IR3 = react(i)%R(3)                        !modif Tabone 12/2017 3Body
          IP1 = react(i)%P(1) ; IP2 = react(i)%P(2)
          IP3 = react(i)%P(3) ; IP4 = react(i)%P(4)
          ! count the number of reaction for each species
          speci(IR1)%nbreac_d = speci(IR1)%nbreac_d + 1
          speci(IR2)%nbreac_d = speci(IR2)%nbreac_d + 1
          speci(IR3)%nbreac_d = speci(IR3)%nbreac_d + 1   !modif Tabone 12/2017 3Body
          speci(IP1)%nbreac_f = speci(IP1)%nbreac_f + 1
          speci(IP2)%nbreac_f = speci(IP2)%nbreac_f + 1
          speci(IP3)%nbreac_f = speci(IP3)%nbreac_f + 1
          speci(IP4)%nbreac_f = speci(IP4)%nbreac_f + 1
       ENDDO
       DO j = 1, nspec
          ALLOCATE(speci(j)%sort_reac_f(speci(j)%nbreac_f))
          ALLOCATE(speci(j)%sort_reac_d(speci(j)%nbreac_d))
          speci(j)%sort_reac_f(1:speci(j)%nbreac_f) = 0
          speci(j)%sort_reac_d(1:speci(j)%nbreac_d) = 0
       ENDDO
       first = .false.
    ENDIF

    !----------------------------------------------------
    ! Sort the reaction rates by increasing order of 
    ! absolute values for each species
    !----------------------------------------------------
    DO j = 1, nspec
       nbf = 0
       nbd = 0
       DO i = 1, Nreact
          IR1 = react(i)%R(1) ; IR2 = react(i)%R(2)
          IR3 = react(i)%R(3)                          !modif Tabone 12/2017 3Body
          IP1 = react(i)%P(1) ; IP2 = react(i)%P(2)
          IP3 = react(i)%P(3) ; IP4 = react(i)%P(4)
          IF(IR1 == j) THEN
             nbd = nbd + 1
             rated(nbd) = -react(i)%rate
             indxd(nbd) = i
          ENDIF
          IF(IR2 == j)  THEN
             nbd = nbd + 1
             rated(nbd) = -react(i)%rate
             indxd(nbd) = i
          ENDIF
          IF(IR3 == j)  THEN                            !modif Tabone 12/2017 3Body
             nbd = nbd + 1
             rated(nbd) = -react(i)%rate
             indxd(nbd) = i
          ENDIF
          IF(IP1 == j)  THEN
             nbf = nbf + 1
             ratef(nbf) = react(i)%rate
             indxf(nbf) = i
          ENDIF
          IF(IP2 == j)  THEN
             nbf = nbf + 1
             ratef(nbf) = react(i)%rate
             indxf(nbf) = i
          ENDIF
          IF(IP3 == j)  THEN
             nbf = nbf + 1
             ratef(nbf) = react(i)%rate
             indxf(nbf) = i
          ENDIF
          IF(IP4 == j)  THEN
             nbf = nbf + 1
             ratef(nbf) = react(i)%rate
             indxf(nbf) = i
          ENDIF
       ENDDO
       CALL SORT_TABLE(nbf,ratef(1:nbf),indxf(1:nbf))
       CALL SORT_TABLE(nbd,rated(1:nbd),indxd(1:nbd))
       ! -----------------------------------------
       ! check the sorting - for debug purposes
       ! -----------------------------------------
       ! IF(speci(j)%name == 'G' .OR. speci(j)%name == 'G-' .OR. speci(j)%name == 'G+') THEN
       !    WRITE(*,*) speci(j)%name, nb
       ! ENDIF
       ! -----------------------------------------
       DO i = 1, nbf
          ! --------------------------------------
          ! check the sorting - for debug purposes
          ! --------------------------------------
          ! IF(speci(j)%name == 'G' .OR. speci(j)%name == 'G-' .OR. speci(j)%name == 'G+') THEN
          !    WRITE(*,*) rate(i), indx(i), speci(React(indx(i))%R(1))%name,  &
          !    speci(React(indx(i))%R(2))%name,  speci(React(indx(i))%P(1))%name,&
          !    speci(React(indx(i))%P(2))%name,  speci(React(indx(i))%P(3))%name,&
          !    speci(React(indx(i))%P(4))%name
          ! ENDIF
          ! --------------------------------------
          speci(j)%sort_reac_f(i) = indxf(i)
          ! WRITE(*,*) speci(j)%sort_reac(i,1)
       ENDDO
       DO i = 1, nbd
          ! --------------------------------------
          ! check the sorting - for debug purposes
          ! --------------------------------------
          ! IF(speci(j)%name == 'G' .OR. speci(j)%name == 'G-' .OR. speci(j)%name == 'G+') THEN
          !    WRITE(*,*) rate(i), indx(i), speci(React(indx(i))%R(1))%name,  &
          !    speci(React(indx(i))%R(2))%name,  speci(React(indx(i))%P(1))%name,&
          !    speci(React(indx(i))%P(2))%name,  speci(React(indx(i))%P(3))%name,&
          !    speci(React(indx(i))%P(4))%name
          ! ENDIF
          ! --------------------------------------
          speci(j)%sort_reac_d(i) = indxd(i)
          ! WRITE(*,*) speci(j)%sort_reac(i,1)
       ENDDO
       ! STOP
    ENDDO

  END SUBROUTINE SORT_REACTIONS


  SUBROUTINE SORT_TABLE(n,table,idx)
    !---------------------------------------------------------------------------
    ! called by :
    !    SORT_REACTIONS
    ! purpose :
    !    Sort a table by increasing order of absolute value
    ! subroutine/function needed :
    ! input variables :
    !    n     -> table size
    !    table -> values
    !    idx   -> indices
    ! ouput variables :
    !    table -> values
    !    idx   -> indices
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT none

    INTEGER,                          INTENT(in)    :: n
    REAL(KIND=dp), DIMENSION(Nreact), INTENT(inout) :: table
    INTEGER,       DIMENSION(Nreact), INTENT(inout) :: idx
    INTEGER                                         :: i, nn, idum
    REAL(KIND=dp)                                   :: tdum
    LOGICAL                                         :: change

    change = .true.
    nn = n
    DO WHILE (change .AND. nn > 1)
       change = .false.
       DO i = 1, nn-1
          IF( ABS(table(i)) > ABS(table(i+1)) ) THEN
             tdum       = table(i)
             table(i)   = table(i+1)
             table(i+1) = tdum
             idum       = idx(i)
             idx(i)     = idx(i+1)
             idx(i+1)   = idum
             change     = .true.
          ENDIF
       ENDDO
       nn = nn - 1
    ENDDO
  END SUBROUTINE SORT_TABLE

END MODULE MODULE_CHEM_REACT
