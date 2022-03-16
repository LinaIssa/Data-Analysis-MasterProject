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

MODULE MODULE_CHEMICAL_SPECIES
  !*****************************************************************************
  !** The module 'MODULE_CHEMICAL_SPECIES' contains variables and subroutines **
  !** related to the set of chemical species :                                **
  !**     * the data type TYPE_SPECY                                          **
  !**     * the vector containing the species                                 **
  !**     * the numbers of species in each fluid                              **
  !**     * the index useful to find one specy in the vector                  **
  !**     * subroutine to write information on one specy                      **
  !**     * function to calculate the chemical formula of one specy           **
  !**     * subroutine to calculate the elemental abundances                  **
  !*****************************************************************************
  ! SC May 06: 
  ! - C60 replaced by G everywhere. 
  ! - Element 'G' with mass = 1amu added in INITIALIZE_ELEMENTS
  ! - Nelements = 12 instead of 11
  USE MODULE_TECHCONFIG

  IMPLICIT NONE
  INCLUDE "precision.f90"

  CHARACTER(len=lenfilename)    :: name_file_speci

  INTEGER(KIND=LONG)            :: Nlines_comment
  INTEGER(KIND=LONG), PARAMETER :: name_length   = 7   ! (7: old, 8: new) length of specy's or element's name
  INTEGER,            PARAMETER :: lenSpeciHname = 27  ! Size of the name of species in human readable format
  INTEGER,            PARAMETER :: lenSpeciID    = 27  ! Size of the ID or of an INCHI key for species
  INTEGER(KIND=LONG), PARAMETER :: Nelements     = 12  ! number of 'basic' elements: 9 + 1(Mg) + 1(D) + 1(G)
  REAL(KIND=DP),      PARAMETER :: fudge         = 1.0

  REAL(KIND = DP),    PARAMETER :: density_limit = 1.0e-20_DP ! smallest density

  ! --------------------------------------------------------------------------------------
  ! data type for one element
  ! - elements are 'basic' components : H, C, N ...
  ! - initialized in INITIALIZE_ELEMENTS
  ! --------------------------------------------------------------------------------------
  TYPE TYPE_ELEMENT
     CHARACTER(len=name_length)  :: name      ! name of the element
     REAL(KIND=DP)               :: mass      ! mass of the element (g)
     REAL(KIND=DP)               :: Dens_init ! initial density (cm-3)
     REAL(KIND=DP)               :: ab        ! abundance
     REAL(KIND=DP)               :: ab_init   ! initial abundance
     REAL(KIND=DP)               :: ab_ref    ! abundance from Anders & Grevesse 1989
     INTEGER,       DIMENSION(3) :: idsp_mn   ! indices   of the two main gas phase carriers that may be used for conservation
     REAL(KIND=DP), DIMENSION(3) :: absp_mn   ! abundance of the two main gas phase carriers  that may be used for conservation
     REAL(KIND=DP)               :: ratio     ! absp_mn(2) / absp_mn(1)
     INTEGER                     :: icon      ! index of the species used for conservation
  END TYPE TYPE_ELEMENT

  TYPE (TYPE_ELEMENT),DIMENSION(Nelements) :: elements

  REAL (KIND=DP)                           :: PAH_abinit
  INTEGER                                  :: PAH_icon
  INTEGER                                  :: charge_icon

  INTEGER(KIND=LONG) :: ind_elem_H  ! index of 'H'  in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_O  ! index of 'O'  in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_C  ! index of 'C'  in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_N  ! index of 'N'  in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_He ! index of 'He' in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_Na ! index of 'Na' in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_Mg ! index of 'Mg' in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_S  ! index of 'S'  in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_Si ! index of 'Si' in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_Fe ! index of 'Fe' in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_D  ! index of 'D'  in the vector 'elements'
  INTEGER(KIND=LONG) :: ind_elem_G  ! index of 'G'  in the vector 'elements'

  ! --------------------------------------------------------------------------------------
  ! data type for one specy
  ! --------------------------------------------------------------------------------------
  TYPE TYPE_SPECY
     CHARACTER(LEN=name_length)   :: name           ! name (ex : 'SiOH+')
     REAL(KIND=DP)                :: density        ! density (cm-3)
     REAL(KIND=DP)                :: density_init   ! density (cm-3)
     REAL(KIND=DP)                :: Dens_old       ! density at last call to DRIVE
     REAL(KIND=DP)                :: Col_dens       ! column density (cm-2)
     INTEGER                      :: charge         ! in elementary charge (+e)
     REAL(KIND=DP)                :: mass           ! mass (g)
     REAL(KIND=DP)                :: enthalpy       ! enthalpy of formation (kCal/mol)
     REAL(KIND=DP)                :: velocity       ! velocity (cm/s)
     REAL(KIND=DP)                :: temperature    ! temperature (K)
     INTEGER(KIND=LONG), DIMENSION(Nelements) :: formula ! chemical formula
     INTEGER(KIND=LONG), DIMENSION(12)        :: useless ! (old system for chemical formula)
     INTEGER(KIND=LONG)           :: index          ! index of the species
     INTEGER(KIND=LONG)           :: varindex       ! index of the species in the vector of variable species.
     INTEGER                      :: variable       ! whether the species is an independent variable.
     INTEGER                      :: icon           ! whether the species is used for the conservation of an element
     INTEGER                      :: fluid          ! to which fluid the species belong (-1,0=none, 1=neutral, 2=positive ions, 3=negative ions)
     CHARACTER(len=lenSpeciHname) :: HName          ! Human readable name
     CHARACTER(len=lenSpeciID)    :: SpecyID        ! ID (CHp for CH+)
     CHARACTER(len=lenSpeciID)    :: INCHI          ! ID (if no INCHI key, then = SpecyID)
     INTEGER                      :: nbreac_f       ! number of chemical reaction (formation)
     INTEGER                      :: nbreac_d       ! number of chemical reaction (destruction)
     INTEGER, DIMENSION(:), ALLOCATABLE :: sort_reac_f ! indices of formation reactions (sorted by increasing absolute rate)
     INTEGER, DIMENSION(:), ALLOCATABLE :: sort_reac_d ! indices of formation reactions (sorted by increasing absolute rate)
  END TYPE TYPE_SPECY

  ! variables read in READ_SPECIES
  INTEGER(KIND=LONG) :: Nspec      ! number of chemical species considered
  INTEGER(KIND=LONG) :: Nspec_var  ! number of chemical species as independent variables
  INTEGER(KIND=LONG) :: Nspec_plus ! Nspec + 5 (photon, CRP, grain, SECPHO, VOISIN)

  ! calculated in INITIALIZE
  INTEGER(KIND=LONG) :: Nneutrals = 0 ! number of neutrals
  INTEGER(KIND=LONG) :: Nions     = 0 ! number of positive ions
  INTEGER(KIND=LONG) :: Nani      = 0 ! number of negative ions
  INTEGER(KIND=LONG) :: Nneg      = 0 ! number of negative species (including electrons)
  INTEGER(KIND=LONG) :: Nongrains = 0 ! number of species in grain mantle (name with a '*')
  INTEGER(KIND=LONG) :: Noncores  = 0 ! number of species in grain cores (name with a '**')

  !-------------------------------------------------------
  ! Vector containing the species, read in READ_SPECIES
  ! - tables start at index zero, e.g. speci(0)%name is ''
  ! - species are between the index 1 and Nspec
  ! - between Nspec+1 and Nspec_plus, we find added 
  !   species (photons, grains,...)
  !-------------------------------------------------------
  TYPE(type_specy), DIMENSION(:), ALLOCATABLE :: speci
  TYPE(type_specy), DIMENSION(:), ALLOCATABLE :: speci_tmp
  INTEGER,          DIMENSION(:), ALLOCATABLE :: speci_index
  INTEGER,          DIMENSION(:), ALLOCATABLE :: speci_index_tmp

  !-------------------------------------------------------
  ! index of some species
  ! initialized in READ_SPECIES
  !-------------------------------------------------------
  INTEGER(KIND=LONG) :: ind_H      = 0
  INTEGER(KIND=LONG) :: ind_H2     = 0
  INTEGER(KIND=LONG) :: ind_He     = 0
  INTEGER(KIND=LONG) :: ind_O      = 0
  INTEGER(KIND=LONG) :: ind_O2     = 0
  INTEGER(KIND=LONG) :: ind_Oplus  = 0
  INTEGER(KIND=LONG) :: ind_N      = 0
  INTEGER(KIND=LONG) :: ind_C      = 0
  INTEGER(KIND=LONG) :: ind_S      = 0
  INTEGER(KIND=LONG) :: ind_Si     = 0
  INTEGER(KIND=LONG) :: ind_H2O    = 0
  INTEGER(KIND=LONG) :: ind_OH     = 0
  INTEGER(KIND=LONG) :: ind_CO     = 0
  INTEGER(KIND=LONG) :: ind_SH     = 0
  INTEGER(KIND=LONG) :: ind_CH     = 0
  INTEGER(KIND=LONG) :: ind_CHp    = 0
  INTEGER(KIND=LONG) :: ind_NH3    = 0
  INTEGER(KIND=LONG) :: ind_SiO    = 0
  INTEGER(KIND=LONG) :: ind_Cplus  = 0
  INTEGER(KIND=LONG) :: ind_Hplus  = 0
  INTEGER(KIND=LONG) :: ind_Siplus = 0
  INTEGER(KIND=LONG) :: ind_Splus  = 0
  INTEGER(KIND=LONG) :: ind_Nplus  = 0
  INTEGER(KIND=LONG) :: ind_Feplus = 0
  INTEGER(KIND=LONG) :: ind_G      = 0
  INTEGER(KIND=LONG) :: ind_Gplus  = 0
  INTEGER(KIND=LONG) :: ind_Gminus = 0
  INTEGER(KIND=LONG) :: ind_e      = 0
  INTEGER(KIND=LONG) :: ind_PHOTON = 0
  INTEGER(KIND=LONG) :: ind_CRP    = 0
  INTEGER(KIND=LONG) :: ind_GRAIN  = 0
  INTEGER(KIND=LONG) :: ind_SECPHO = 0
  INTEGER(KIND=LONG) :: ind_VOISIN = 0
  INTEGER(KIND=LONG) :: ind_D      = 0
  INTEGER(KIND=LONG) :: ind_PAH0   = 0
  INTEGER(KIND=LONG) :: ind_PAHm   = 0
  INTEGER(KIND=LONG) :: ind_PAHp   = 0
  INTEGER(KIND=LONG) :: ind_Ccore  = 0
  INTEGER(KIND=LONG) :: ind_Ocore  = 0

  !-------------------------------------------------------
  ! mass of some species
  ! initialized in READ_SPECIES
  !-------------------------------------------------------
  REAL(KIND=DP) :: mass_H      = 0.0_DP
  REAL(KIND=DP) :: mass_H2     = 0.0_DP
  REAL(KIND=DP) :: mass_He     = 0.0_DP
  REAL(KIND=DP) :: mass_O      = 0.0_DP
  REAL(KIND=DP) :: mass_O2     = 0.0_DP
  REAL(KIND=DP) :: mass_Oplus  = 0.0_DP
  REAL(KIND=DP) :: mass_N      = 0.0_DP
  REAL(KIND=DP) :: mass_C      = 0.0_DP
  REAL(KIND=DP) :: mass_S      = 0.0_DP
  REAL(KIND=DP) :: mass_Si     = 0.0_DP
  REAL(KIND=DP) :: mass_H2O    = 0.0_DP
  REAL(KIND=DP) :: mass_OH     = 0.0_DP
  REAL(KIND=DP) :: mass_CO     = 0.0_DP
  REAL(KIND=DP) :: mass_SH     = 0.0_DP
  REAL(KIND=DP) :: mass_CH     = 0.0_DP
  REAL(KIND=DP) :: mass_CHp    = 0.0_DP
  REAL(KIND=DP) :: mass_NH3    = 0.0_DP
  REAL(KIND=DP) :: mass_SiO    = 0.0_DP
  REAL(KIND=DP) :: mass_Cplus  = 0.0_DP
  REAL(KIND=DP) :: mass_Hplus  = 0.0_DP
  REAL(KIND=DP) :: mass_Siplus = 0.0_DP
  REAL(KIND=DP) :: mass_G      = 0.0_DP
  REAL(KIND=DP) :: mass_Gplus  = 0.0_DP
  REAL(KIND=DP) :: mass_Gminus = 0.0_DP
  REAL(KIND=DP) :: mass_D      = 0.0_DP
  REAL(KIND=DP) :: mass_Splus  = 0.0_DP
  REAL(KIND=DP) :: mass_Nplus  = 0.0_DP
  REAL(KIND=DP) :: mass_Feplus = 0.0_DP
  REAL(KIND=DP) :: mass_e      = 0.0_DP

  !-------------------------------------------------------
  ! density of some species
  ! calculated in DIFFUN
  !-------------------------------------------------------
  REAL(KIND=DP) :: Dens_H
  REAL(KIND=DP) :: Dens_H2
  REAL(KIND=DP) :: Dens_He
  REAL(KIND=DP) :: Dens_O
  REAL(KIND=DP) :: Dens_Oplus
  REAL(KIND=DP) :: Dens_N
  REAL(KIND=DP) :: Dens_C
  REAL(KIND=DP) :: Dens_S
  REAL(KIND=DP) :: Dens_Si
  REAL(KIND=DP) :: Dens_H2O
  REAL(KIND=DP) :: Dens_OH
  REAL(KIND=DP) :: Dens_CO
  REAL(KIND=DP) :: Dens_13CO ! added by Tabone 02/2018 to compute 13CO colling from K&N
  REAL(KIND=DP) :: Dens_CH
  REAL(KIND=DP) :: Dens_SH
  REAL(KIND=DP) :: Dens_NH3
  REAL(KIND=DP) :: Dens_Cplus
  REAL(KIND=DP) :: Dens_Siplus
  REAL(KIND=DP) :: Dens_Hplus
  REAL(KIND=DP) :: Dens_Splus
  REAL(KIND=DP) :: Dens_Nplus
  REAL(KIND=DP) :: Dens_Feplus
  REAL(KIND=DP) :: Dens_G
  REAL(KIND=DP) :: Dens_Gplus
  REAL(KIND=DP) :: Dens_Gminus
  REAL(KIND=DP) :: Dens_cor
  REAL(KIND=DP) :: Dens_e

  !-------------------------------------------------------
  ! index of beginning and end of each type of species
  ! calculated in INITIALIZE
  !-------------------------------------------------------
  INTEGER(KIND=LONG) :: b_neu=0, e_neu=0 ! neutrals
  INTEGER(KIND=LONG) :: b_ion=0, e_ion=0 ! positive ions
  INTEGER(KIND=LONG) :: b_ani=0, e_ani=0 ! negative ions
  INTEGER(KIND=LONG) :: b_neg=0, e_neg=0 ! negative species (including electrons)
  INTEGER(KIND=LONG) :: b_gra=0, e_gra=0 ! species on grain mantles
  INTEGER(KIND=LONG) :: b_cor=0, e_cor=0 ! species on grain mantles

  !-------------------------------------------------------
  ! read/write format for one specy
  !-------------------------------------------------------
  CHARACTER(len=*), PRIVATE, PARAMETER :: format_specy = '(I3,2X,A7,2X,2I2,10I1,ES10.3,F10.3,1x,a1)'

  ! --------------------------------------------------------------------------------------
  ! data type for UV transitions
  ! --------------------------------------------------------------------------------------
  TYPE TYPE_UVDATA
     INTEGER                                  :: nuv  ! nb of UV transitions
     INTEGER,       DIMENSION(:), ALLOCATABLE :: vl   ! v lower level
     INTEGER,       DIMENSION(:), ALLOCATABLE :: jl   ! j lower level
     INTEGER,       DIMENSION(:), ALLOCATABLE :: vu   ! v upper level
     INTEGER,       DIMENSION(:), ALLOCATABLE :: dj   ! Delta J (J_up = J_low + Delta J)
     INTEGER,       DIMENSION(:), ALLOCATABLE :: iw   ! index of line in wavelength grid
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: osc  ! oscillator strength
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: wlg  ! wavelength (Angstrom)
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: ilft ! inverse life time (s-1)
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: pbdi ! Upper level dissociation probability
     REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: fgk  ! fgk radiation field for the transition
  END TYPE TYPE_UVDATA

  ! --------------------------------------------------------------------------------------
  ! variables from "precision.f90" are private
  ! --------------------------------------------------------------------------------------
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE INITIALIZE_ELEMENTS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    initialised the elements and their abundances from
    !    Anders & Grevesse 1989.
    !    elements are 'basic' components of molecules : H, C, N ...
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    elements
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : AMU, Zero

    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i

    ! initialization
    elements(:)%name    = ''
    elements(:)%mass    = Zero
    elements(:)%ab      = Zero
    elements(:)%ab_ref  = Zero
    elements(:)%ab_init = Zero

    ! first with mass in amu
    i = 1 ; elements(i)%name = 'H' ; elements(i)%mass = 1.00797_DP ; elements(i)%ab_ref = 1.00e+00_DP
    ind_elem_H  = i
    i=i+1 ; elements(i)%name = 'D' ; elements(i)%mass = 2.00000_DP ; elements(i)%ab_ref = 2.00e-05_DP
    ind_elem_D  = i
    i=i+1 ; elements(i)%name = 'He'; elements(i)%mass = 4.00260_DP ; elements(i)%ab_ref = 1.00e-01_DP
    ind_elem_He = i
    i=i+1 ; elements(i)%name = 'C' ; elements(i)%mass = 12.0111_DP ; elements(i)%ab_ref = 3.62e-04_DP
    ind_elem_C  = i
    i=i+1 ; elements(i)%name = 'N' ; elements(i)%mass = 14.0067_DP ; elements(i)%ab_ref = 1.12e-04_DP
    ind_elem_N  = i
    i=i+1 ; elements(i)%name = 'O' ; elements(i)%mass = 15.9994_DP ; elements(i)%ab_ref = 8.53e-04_DP
    ind_elem_O  = i
    i=i+1 ; elements(i)%name = 'Na'; elements(i)%mass = 22.9898_DP ; elements(i)%ab_ref = 2.06e-06_DP
    ind_elem_Na = i
    i=i+1 ; elements(i)%name = 'Mg'; elements(i)%mass = 24.3000_DP ; elements(i)%ab_ref = 3.85e-05_DP
    ind_elem_Mg = i
    i=i+1 ; elements(i)%name = 'Si'; elements(i)%mass = 28.0860_DP ; elements(i)%ab_ref = 3.58e-05_DP
    ind_elem_Si = i
    i=i+1 ; elements(i)%name = 'S' ; elements(i)%mass = 32.0640_DP ; elements(i)%ab_ref = 1.85e-05_DP
    ind_elem_S  = i
    i=i+1 ; elements(i)%name = 'Fe'; elements(i)%mass = 55.8470_DP ; elements(i)%ab_ref = 3.23e-05_DP
    ind_elem_Fe = i
    ! SC May 06: replace C60 by element G (generic grain). small mass to keep muI unaffected
    i=i+1 ; elements(i)%name = 'G' ; elements(i)%mass = 1._DP ;elements(i)%ab_ref = 1.00D-10
    ind_elem_G  = i

    ! conversion of mass : AMU -> grammes
    elements(:)%mass = AMU * elements(:)%mass

  END SUBROUTINE INITIALIZE_ELEMENTS

  SUBROUTINE READ_SPECIES
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     read chemical species from the input file and computes indexes of
    !     current species, add e-, photon, CRP, grain, SECPHO
    !     method : the file is read twice (first to count the number of species)
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     speci, mass_X, ind_X (X is a specy)
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS,     ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero, me
    USE MODULE_PHYS_VAR,  ONLY : nH_init

    IMPLICIT NONE
    CHARACTER(len = 80)               :: dummy_line
    CHARACTER(len =  1)               :: variable_char
    CHARACTER(len = name_length)      :: charact
    INTEGER(KIND = LONG)              :: i, error, ii,j
    INTEGER(KIND = LONG)              :: file_speci
    CHARACTER, DIMENSION(name_length) :: letters

    CHARACTER(len=name_length)        :: espece     ! fortran species name
    CHARACTER(len=lenSpeciHname)      :: especeHN   ! human readable name of species
    CHARACTER(len=lenSpeciID)         :: especeID   ! ID of species
    CHARACTER(len=lenSpeciID)         :: inchi = "" ! inchi key of species

    !-----------------------------------------------------
    ! opening file
    !-----------------------------------------------------
    file_speci = GET_FILE_NUMBER()
    name_file_speci = TRIM(data_dir) // TRIM(specfile)
    OPEN(file_speci, file = name_file_speci, status='OLD',&
         access = 'SEQUENTIAL', form = 'FORMATTED', action = 'READ')

    !-----------------------------------------------------
    ! Count number of comment lines 
    !-----------------------------------------------------
    dummy_line = '!'
    nlines_comment = -1
    DO WHILE (dummy_line(1:1) == '!')
       nlines_comment = nlines_comment+1
       READ(file_speci,'(A80)') dummy_line
    ENDDO
    REWIND(file_speci)

    !-----------------------------------------------------
    ! Skip and print comments
    !-----------------------------------------------------
    DO i=1,nlines_comment
       READ(file_speci,'(A80)') dummy_line
       PRINT *, dummy_line
    END DO

    !-----------------------------------------------------
    ! counts the number of chemical species
    !-----------------------------------------------------
    Nspec = 0 ; error = 0
    ! stop at error or end of file
    DO WHILE (error == 0)
       charact = ''
       ! READ(file_speci,'(A1)',iostat=error) charact        ! new
       READ(file_speci,format_specy,iostat=error) ii, charact ! old
       SELECT CASE (charact(1:1))
       ! if the character is a capital letter, this is a specy
       CASE ('A':'Z')
          Nspec = Nspec + 1
       ! if it's not a specy, nor end of file -> error
       CASE DEFAULT
          IF (error >= 0 ) STOP "*** WARNING1, error in READ_SPECIES"
       END SELECT
    ENDDO
    ! add e- as a species
    Nspec = Nspec + 1

    !-----------------------------------------------------
    ! allocation and initialization of the vector species
    !-----------------------------------------------------
    Nspec_plus = Nspec + 5        ! add photon, grains, CRP, SECPHO, VOISIN
    ALLOCATE(speci(0:Nspec_plus)) ! start a zero : undetermined specy
    ALLOCATE(speci_tmp(0:Nspec_plus)) ! start a zero : undetermined specy
    ALLOCATE(speci_index_tmp(1:Nspec))
    speci(0:Nspec_plus)%name          = ''
    speci(0:Nspec_plus)%charge        = 0
    speci(0:Nspec_plus)%mass          = Zero
    speci(0:Nspec_plus)%enthalpy      = Zero
    speci(0:Nspec_plus)%density       = Zero
    speci(0:Nspec_plus)%density_init  = Zero
    speci(0:Nspec_plus)%Dens_old      = Zero
    speci(0:Nspec_plus)%Col_dens      = Zero
    speci(0:Nspec_plus)%velocity      = Zero
    speci(0:Nspec_plus)%temperature   = Zero
    speci(0:Nspec_plus)%index         = 0
    speci(0:Nspec_plus)%varindex      = 0
    speci(0:Nspec_plus)%variable      = 0
    speci(0:Nspec_plus)%HName         = ''
    speci(0:Nspec_plus)%SpecyID       = ''
    speci(0:Nspec_plus)%INCHI         = ''
    speci(0:Nspec_plus)%nbreac_f      = 0
    speci(0:Nspec_plus)%nbreac_d      = 0
    speci(0:Nspec_plus)%fluid         = -99
    DO i = 1, Nelements
       speci(0:Nspec_plus)%formula(i) = 0
    ENDDO
    speci_tmp(0:Nspec_plus) = speci(0:Nspec_plus)

    REWIND(file_speci)

    !-----------------------------------------------------
    ! go directly to the chemical species
    !-----------------------------------------------------
    DO i=1, Nlines_comment
       READ(file_speci,'(A1)') charact
    END DO

    i         = 0
    error     = 0
    Nneutrals = 0 
    Nions     = 0
    Nani      = 0 
    Nneg      = 0 
    Nongrains = 0
    Noncores  = 0
    !-----------------------------------------------------
    ! stop reading after Nspec - 1 or error
    !-----------------------------------------------------
    DO WHILE (i < Nspec - 1 .AND. error == 0)
       i = i + 1
       READ(file_speci,format_specy, iostat=error) &
          & ii, speci_tmp(i)%name, speci_tmp(i)%useless, &
          & speci_tmp(i)%density, speci_tmp(i)%enthalpy, variable_char
       letters = ''
       DO j = 1, name_length
          letters(j) = speci_tmp(i)%name(j:j)
       ENDDO

       !--------------------------------------------------
       ! Save initial abundance
       ! Convert inital densities to cm-3
       !--------------------------------------------------
       speci_tmp(i)%density_init=speci_tmp(i)%density
       speci_tmp(i)%density = speci_tmp(i)%density * nH_init

       IF (error == 0) THEN
          speci_tmp(i)%index = i
          !-----------------------------------------------
          ! - number of each element in the molecule
          ! - mass is determined from the formula and
          !   the mass of each element (g)
          ! - mass is modified depending on the charge of
          !   the species
          !-----------------------------------------------
          speci_tmp(i)%formula = CHEMICAL_FORMULA(speci_tmp(i)%name)
          speci_tmp(i)%charge  = CHARGE_FORMULA(speci_tmp(i)%name)
          speci_tmp(i)%mass    = DOT_PRODUCT ( DBLE(elements(1:Nelements)%mass), &
                                               DBLE(speci_tmp(i)%formula(1:Nelements)) )
          speci_tmp(i)%mass    = speci_tmp(i)%mass - speci_tmp(i)%charge * me
          IF ( any(letters == '+') .and. (fudge.ne.1.0) ) then
             PRINT *, speci_tmp(i)%name, ' IS FUDGED. fudge=', fudge
             speci_tmp(i)%mass = speci_tmp(i)%mass*fudge
          ENDIF

          !-----------------------------------------------
          ! Names for ism services - ID, HRN, INCHI, ...
          !-----------------------------------------------
          espece = speci_tmp(i)%name
          CALL SPEID(espece, especeID, especeHN)      ! Transform in upper case and create an ID
          speci_tmp(i)%HName   = especeHN             ! Name of species in upper case, ie : CH+
          speci_tmp(i)%SpecyID = especeID             ! ID of species                  ie : CHp
          IF ( (TRIM(inchi) == "NOINCHIFOUND").OR.(TRIM(inchi) == "") ) THEN
             speci_tmp(i)%INCHI  = especeID           ! ID of species (home made ID if no INCHI)
          ELSE
             speci_tmp(i)%INCHI  = inchi              ! INCHI key of species
          ENDIF

          !-----------------------------------------------
          ! Possibility to fix the abundance of a species
          ! Be careful -> if the species is in the chemical
          ! network, such assumption breaks conservation
          !-----------------------------------------------
          SELECT CASE(variable_char)
          CASE('0','C','c','d','D','e','E')
             speci_tmp(i)%variable =  0
             PRINT *, 'WARNING: specy ',speci_tmp(i)%name,' is assumed constant.'
          CASE DEFAULT
             speci_tmp(i)%variable = 1
          END SELECT

          !-----------------------------------------------
          ! Identification of the fluid to which the
          ! species is associated 
          ! Number of species in each fluid
          !-----------------------------------------------
          IF      ( speci_tmp(i)%charge == 0 ) THEN
             IF      ( INDEX(speci_tmp(i)%name,'**') > 0 ) THEN
                Noncores  = Noncores + 1
                speci_tmp(i)%fluid = -1
             ELSE IF ( INDEX(speci_tmp(i)%name,'*')  > 0 ) THEN
                NonGrains = NonGrains + 1
                speci_tmp(i)%fluid =  0
             ! CAREFUL - BG2017 -> choosing G as an ionized 
             ! species is not without consequences
             ! In fact, bad choice, it's better to consider G 
             ! as a neutral species and to add the contributions
             ! when needed (on dynamics, momentum exchange, ...
             ! ELSE IF ( TRIM(speci_tmp(i)%name) == 'G' ) THEN
             !    Nions = Nions + 1
             !    speci_tmp(i)%fluid =  2
             ELSE
                Nneutrals = Nneutrals + 1
                speci_tmp(i)%fluid =  1
             ENDIF
          ELSE IF ( speci_tmp(i)%charge > 0 ) THEN
             Nions = Nions + 1
             speci_tmp(i)%fluid = 2
          ELSE IF ( speci_tmp(i)%charge < 0 ) THEN
             Nani  = Nani  + 1
             Nneg  = Nneg  + 1
             speci_tmp(i)%fluid = 3
          ENDIF

       ELSE
          IF (error > 0) THEN
             PRINT *, 'error=',error
             PRINT *, variable_char
             STOP "*** WARNING2, error in READ_SPECIES"
          ENDIF
       ENDIF
    ENDDO

    !-----------------------------------------------------
    ! file closure
    !-----------------------------------------------------
    CLOSE(file_speci)

    !-----------------------------------------------------
    ! Re-ordering of species by fluid type
    !-----------------------------------------------------
    j = 0
    ! Neutrals
    DO i = 1, Nspec - 1
       IF ( speci_tmp(i)%fluid ==  1 ) THEN
          j = j + 1
          speci(j) = speci_tmp(i)
          speci(j)%index = j
       ENDIF
    ENDDO
    e_neu = j
    ! Grain mantles
    DO i = 1, Nspec - 1
       IF ( speci_tmp(i)%fluid ==  0 ) THEN
          j = j + 1
          speci(j) = speci_tmp(i)
          speci(j)%index = j
       ENDIF
    ENDDO
    e_gra = j
    ! Grain cores
    DO i = 1, Nspec - 1
       IF ( speci_tmp(i)%fluid == -1 ) THEN
          j = j + 1
          speci(j) = speci_tmp(i)
          speci(j)%index = j
       ENDIF
    ENDDO
    e_cor = j
    ! Positive ions
    DO i = 1, Nspec - 1
       IF ( speci_tmp(i)%fluid ==  2 ) THEN
          j = j + 1
          speci(j) = speci_tmp(i)
          speci(j)%index = j
       ENDIF
    ENDDO
    e_ion = j
    ! Negative ions
    DO i = 1, Nspec - 1
       IF ( speci_tmp(i)%fluid ==  3 ) THEN
          j = j + 1
          speci(j) = speci_tmp(i)
          speci(j)%index = j
       ENDIF
    ENDDO
    e_ani = j
    e_neg = j

    ! Compute the starting indices
    b_neu = e_neu - Nneutrals + 1
    b_cor = e_cor - Noncores  + 1
    b_gra = e_gra - NonGrains + 1
    b_ion = e_ion - Nions     + 1
    b_ani = e_ani - Nani      + 1
    b_neg = e_neg - Nneg      + 1

    ! No need for speci_tmp anymore
    DEALLOCATE (speci_tmp)

    !-----------------------------------------------------
    ! Save indexes & masses of several chemical species
    !-----------------------------------------------------
    DO i = 1, Nspec - 1
       SELECT CASE (TRIM(speci(i)%name))
       CASE ('H')
          ind_H       = i
          mass_H      = speci(i)%mass
       CASE ('D')
          ind_D       = i
          mass_D      = speci(i)%mass
       CASE ('H2')
          ind_H2      = i
          mass_H2     = speci(i)%mass
       CASE ('He')
          ind_He      = i
          mass_He     = speci(i)%mass
       CASE ('O')
          ind_O       = i
          mass_O      = speci(i)%mass
       CASE ('O2')
          ind_O2      = i
          mass_O2     = speci(i)%mass
       CASE ('O+')
          ind_Oplus   = i
          mass_Oplus  = speci(i)%mass
       CASE ('N')
          ind_N       = i
          mass_N      = speci(i)%mass
       CASE ('C')
          ind_C       = i
          mass_C      = speci(i)%mass
       CASE ('C**')
          ind_Ccore   = i
       CASE ('O**')
          ind_Ocore   = i
       CASE ('Si')
          ind_Si      = i
          mass_Si     = speci(i)%mass
       CASE ('H2O')
          ind_H2O     = i
          mass_H2O    = speci(i)%mass
       CASE ('OH')
          ind_OH      = i
          mass_OH     = speci(i)%mass
       CASE ('CO')
          ind_CO      = i
          mass_CO     = speci(i)%mass
       CASE ('S')
          ind_S       = i
          mass_S      = speci(i)%mass
       CASE ('SH')
          ind_SH      = i
          mass_SH     = speci(i)%mass
       CASE ('CH')
          ind_CH      = i
          mass_CH     = speci(i)%mass
       CASE ('CH+')
          ind_CHp     = i
          mass_CHp    = speci(i)%mass
       CASE ('NH3')
          ind_NH3     = i
          mass_NH3    = speci(i)%mass
       CASE ('SiO')
          ind_SiO     = i
          mass_SiO    = speci(i)%mass
       CASE ('G')
          ind_G       = i
          mass_G      = speci(i)%mass
       CASE ('C+')
          ind_Cplus   = i
          mass_Cplus  = speci(i)%mass
       CASE ('H+')
          ind_Hplus   = i
          mass_Hplus  = speci(i)%mass
       CASE ('Si+')
          ind_Siplus  = i
          mass_Siplus = speci(i)%mass
       CASE ('S+')
          ind_Splus   = i
          mass_Splus  = speci(i)%mass
       CASE ('N+')
          ind_Nplus   = i
          mass_Nplus  = speci(i)%mass
       CASE ('C54H18')
          ind_PAH0    = i
       CASE ('C54H18+')
          ind_PAHp    = i
       CASE ('C54H18-')
          ind_PAHm    = i
       CASE ('Fe+')
          ind_Feplus  = i
          mass_Feplus = speci(i)%mass
       CASE ('G+')
          ind_Gplus   = i
          mass_Gplus  = speci(i)%mass
       CASE ('G-')
          ind_Gminus  = i
          mass_Gminus = speci(i)%mass
       END SELECT
    ENDDO
    
    !-----------------------------------------------------
    ! count number of variable species and maps into array
    !-----------------------------------------------------
    Nspec_var = 0
    DO i = 1, Nspec - 1
       IF ( speci(i)%variable == 0 ) THEN
          speci(i)%varindex = -1 ! Should produce out of bounds in case it is used
       ELSE
          Nspec_var = Nspec_var+1
          speci_index_tmp(Nspec_var) = i  ! maps variable species into speci array.
          speci(i)%varindex = Nspec_var
       ENDIF
    ENDDO

    !-----------------------------------------------------
    ! addition of e- to the list of species
    ! - charge   = -1
    ! - mass     = me
    ! - enthalpy = 0.0
    !-----------------------------------------------------
    i = Nspec
    ind_e = i
    Nspec_var = Nspec_var + 1
    ! e-, initial density is set in INITIALIZE
    speci(ind_e)%name     = 'ELECTR'
    speci(ind_e)%Hname    = "e-"
    speci(ind_e)%SpecyID  = "electr"
    speci(ind_e)%INCHI    = "ELECTR"
    speci(ind_e)%index    = i
    speci(ind_e)%charge   = -1
    speci(ind_e)%mass     = me
    speci(ind_e)%variable = 1
    speci(ind_e)%varindex = Nspec_var
    speci_index_tmp(Nspec_var) = i
    mass_e                = me
    Nneg                  = Nneg  + 1
    e_neg                 = i
    speci(ind_e)%fluid    = 3

    ! -----------------------------------------------------
    ! Check on screen chemical species
    ! BG2017
    ! -----------------------------------------------------
    ! DO i = 1, Nspec
    !    WRITE(*,'(A8,4(I3,2X),2(ES12.5,2X))') &
    !    speci(i)%name, speci(i)%index, speci(i)%variable, speci(i)%varindex, speci(i)%charge, speci(i)%mass, speci(i)%density
    ! ENDDO

    !-----------------------------------------------------
    ! Summary of input chemical species
    !-----------------------------------------------------
    PRINT *, 'Total number of species:', Nspec
    PRINT *, 'Found Nspec_var = ', Nspec_var,' variable species.'
    ! speci_index with the right dimension
    ALLOCATE(speci_index(Nspec_var))
    speci_index = speci_index_tmp(1:Nspec_var)
    DEALLOCATE(speci_index_tmp)

    !-----------------------------------------------------
    ! avoid small numbers : lower limit to the densities
    !-----------------------------------------------------
    WHERE (speci(1:Nspec-1)%density < density_limit)
       speci(1:Nspec-1)%density = density_limit
    END WHERE
    
    !-----------------------------------------------------
    ! Compute initial electron density
    !-----------------------------------------------------
    speci(Nspec)%density = SUM( speci(1:Nspec-1)%density * speci(1:Nspec-1)%charge )
    speci(Nspec)%density_init = speci(Nspec)%density
    IF (speci(Nspec)%density < density_limit) THEN
       speci(Nspec)%density = density_limit
    ENDIF

    !-----------------------------------------------------
    ! addition of photon, CRP, GRAIN, SECPHO, VOISIN
    ! charge, mass, enthalpy = 0.0
    !-----------------------------------------------------
    i = Nspec
    ! GRAIN, initial density is set in INITIALIZE
    i = i + 1
    speci(i)%name     = 'GRAIN'
    speci(i)%Hname    = "grain"
    speci(i)%SpecyID  = "grain"
    speci(i)%INCHI    = "GRAIN"
    speci(i)%index    = i
    speci(i)%density  = 1._DP
    ind_GRAIN         = i
    ! photon, density = constant = 1._DP (cf CHEMISTRY)
    i = i + 1
    speci(i)%name     = 'PHOTON'
    speci(i)%Hname    = "photon"
    speci(i)%SpecyID  = "photon"
    speci(i)%INCHI    = "PHOTON"
    speci(i)%index    = i
    speci(i)%density  = 1._DP
    ind_PHOTON        = i
    ! CRP (cosmique ray proton), density = constant = 1._DP (cf CHEMISTRY)
    i = i + 1
    speci(i)%name     = 'CRP'
    speci(i)%Hname    = "crp"
    speci(i)%SpecyID  = "crp"
    speci(i)%INCHI    = "CRP"
    speci(i)%index    = i
    speci(i)%density  = 1._DP
    ind_CRP           = i
    ! SECPHO (secondary photon), initial = constant = 1._DP (cf CHEMISTRY)
    i = i + 1
    speci(i)%name     = 'SECPHO'
    speci(i)%Hname    = "secpho"
    speci(i)%SpecyID  = "secpho"
    speci(i)%INCHI    = "secpho"
    speci(i)%index    = i
    speci(i)%density  = 1._DP
    ind_SECPHO        = i
    ! VOISIN (dummy species for thermal desorption), initial = constant = 1._DP (cf CHEMISTRY)
    i = i + 1
    speci(i)%name     = 'VOISIN'
    speci(i)%Hname    = "voisin"
    speci(i)%SpecyID  = "voisin"
    speci(i)%INCHI    = "VOISIN"
    speci(i)%index    = i
    speci(i)%density  = 1._DP
    ind_VOISIN        = i

  END SUBROUTINE READ_SPECIES

  SUBROUTINE CHECK_SPECIES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    tests if the set of species doesn't have identical species
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i, j

    DO i=1, Nspec
       DO j=i+1,Nspec
          IF (speci(i)%name==speci(j)%name) THEN
             WRITE(*,*) "*** WARNING, species ",i," and ",j," are identical"
             STOP
          END IF
       END DO
    END DO
  END SUBROUTINE CHECK_SPECIES

  SUBROUTINE WRITE_SPECY(num,one_specy)
    !--------------------------------------------------------------------------
    ! called by :
    !     WRITE_INFO
    ! purpose :
    !     write informations about one chemical specy
    ! subroutine/function needed :
    ! input variables :
    !     * num -> file number where to write the informations
    !     * one_specy -> (type TYPE_SPECY) chemical specy
    ! output variables :
    ! results :
    !--------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : nH
    IMPLICIT NONE
    TYPE(type_specy), INTENT(in) :: one_specy
    INTEGER(KIND=LONG),INTENT(in) :: num
    character(1)::charact
    write(charact,'(i1)')one_specy%variable

    WRITE(num,format_specy) &
         one_specy%index, &
         one_specy%name, &
         one_specy%useless, &
         one_specy%density/nH, &
         one_specy%enthalpy, &
         charact

  END SUBROUTINE WRITE_SPECY


  SUBROUTINE ELEMENTAL_ABUNDANCES
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    computes elemental abundances in (H,C,N,...)
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    elements%ab
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i

    ! computation of abundances, using speci(i)%formula and speci(i)%density
    DO i=1, Nelements
       elements(i)%ab = SUM(DBLE(speci(1:Nspec)%formula(i)) * speci(1:Nspec)%density)
    END DO
    elements(1:Nelements)%Dens_init = elements(1:Nelements)%ab
    ! normalization to H abundance
    elements(1:Nelements)%ab = elements(1:Nelements)%ab / elements(ind_elem_H)%ab

    PAH_abinit = speci(ind_PAH0)%density + speci(ind_PAHp)%density + speci(ind_PAHm)%density
    PAH_abinit = PAH_abinit / elements(ind_elem_H)%Dens_init

    ! check the calculated and the read masses
    ! DO i=1,Nspec
    !    speci(i)%mass2=DOT_PRODUCT(DBLE(elements(:)%mass),DBLE(speci(i)%formula))
    !    IF (abs((speci(i)%mass-speci(i)%mass2)/speci(i)%mass) > 1.D-5) &
    !         write(*,*)i,speci(i)%name,speci(i)%mass,speci(i)%mass2
    ! END DO

  END SUBROUTINE ELEMENTAL_ABUNDANCES


  SUBROUTINE MAIN_CARRIER
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    find the two most abundant gas phase carrier of every element.
    !    PAHs are not considered as gas phase species and are therefore
    !    excluded from the search
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    elements%idsp_mn
    !    elements%absp_mn
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR
    IMPLICIT NONE
    INTEGER(KIND=LONG)       :: i,j
    INTEGER(KIND=LONG)       :: dum
    REAL(KIND=DP), PARAMETER :: tiny_local = 1e-100_dp

    ! computation of abundances and indices using speci(i)%formula and speci(i)%density
    DO i=1, Nelements
       elements(i)%absp_mn(1) = 0.0_dp
       elements(i)%absp_mn(2) = 0.0_dp
       elements(i)%absp_mn(3) = 0.0_dp
       elements(i)%idsp_mn(1) = 0
       elements(i)%idsp_mn(2) = 0
       elements(i)%idsp_mn(3) = 0

       DO j = 1, Nspec
          dum = speci(j)%formula(i)
          IF (dum == 0) CYCLE
          ! skip species on grains (either mantles or cores)
          IF ( F_CONS == 3 ) THEN
             IF (INDEX(speci(j)%name, "*") /= 0) CYCLE
          ENDIF
          ! skip PAHs
          IF ((j == ind_PAH0).OR.(j == ind_PAHp).OR.(j == ind_PAHm)) CYCLE
          ! skip species which are assumed constant
          IF (speci(j)%variable == 0) CYCLE
          IF (speci(j)%density * dum > elements(i)%absp_mn(1)) THEN
             elements(i)%absp_mn(3) = elements(i)%absp_mn(2)
             elements(i)%idsp_mn(3) = elements(i)%idsp_mn(2)
             elements(i)%absp_mn(2) = elements(i)%absp_mn(1)
             elements(i)%idsp_mn(2) = elements(i)%idsp_mn(1)
             elements(i)%absp_mn(1) = speci(j)%density * dum
             elements(i)%idsp_mn(1) = j
          ELSE IF (speci(j)%density * dum > elements(i)%absp_mn(2)) THEN
             elements(i)%absp_mn(3) = elements(i)%absp_mn(2)
             elements(i)%idsp_mn(3) = elements(i)%idsp_mn(2)
             elements(i)%absp_mn(2) = speci(j)%density * dum
             elements(i)%idsp_mn(2) = j
          ELSE IF (speci(j)%density * dum > elements(i)%absp_mn(3)) THEN
             elements(i)%absp_mn(3) = speci(j)%density * dum
             elements(i)%idsp_mn(3) = j
          ENDIF
       ENDDO
       IF (elements(i)%absp_mn(1) > tiny_local) THEN
          elements(i)%ratio = elements(i)%absp_mn(2) / elements(i)%absp_mn(1)
       ELSE
          elements(i)%ratio = 0.0_dp
       ENDIF
    ENDDO

  END SUBROUTINE MAIN_CARRIER


  SUBROUTINE INDEX_CONSERVATION
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    find species whose abundance are computed with conservation equations
    !    H2 is used for the conservation of its level population and therefore
    !    cannot be used for the conservation of elements.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    speci%icon
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: i, j, l, ii, jj
    INTEGER(KIND=LONG) :: isp1, isp2, isp3
    REAL(KIND=DP)      :: dum
    INTEGER(KIND=LONG) :: id1,id2

    ! reset all conservation indices
    DO j = 1, Nspec
       speci(j)%icon = 0
    ENDDO
    DO i=1, Nelements
       elements(i)%icon = 0
    ENDDO
    PAH_icon    = 0
    charge_icon = 0

    ! --------------------------------------------
    ! find the index for conservation of charge
    ! --------------------------------------------
    dum = 0.0_dp
    DO i = 1, Nspec
       IF ( speci(i)%charge >= 0 ) CYCLE
       IF ( speci(i)%density > dum) THEN
          dum = speci(i)%density
          charge_icon = i
       ENDIF
    ENDDO

    DO i=1, Nelements
       isp1 = elements(i)%idsp_mn(1)
       isp2 = elements(i)%idsp_mn(2)
       isp3 = elements(i)%idsp_mn(3)
       IF (isp1 == 0) CYCLE

       ii = 0
       DO l = 1, i-1
          IF (speci(isp1)%icon == l) ii = l
       ENDDO
       jj = 0
       DO l = ii+1, i-1
          IF (speci(isp1)%icon == l) jj = l
       ENDDO
  
       IF (ii == 0) THEN
          IF      ( isp1 /= ind_H2 .AND. isp1 /= charge_icon ) THEN
             speci(isp1)%icon = i
             elements(i)%icon = isp1
          ELSE IF ( isp2 /= ind_H2 .AND. isp2 /= charge_icon ) THEN
             speci(isp2)%icon = i
             elements(i)%icon = isp2
          ELSE
             speci(isp3)%icon = i
             elements(i)%icon = isp3
          ENDIF
       ELSE
          IF (jj == 0) THEN
             ! IF (elements(i)%ratio > elements(ii)%ratio) THEN
             DO l = 1, i-1
                IF (speci(isp2)%icon == l) jj = l
             ENDDO
             IF ( jj == 0 ) THEN
                IF ( isp2 /= charge_icon ) THEN
                   speci(isp2)%icon = i
                   elements(i)%icon = isp2
                ELSE
                   speci(isp3)%icon = i
                   elements(i)%icon = isp3
                ENDIF
             ELSE
                IF ( isp3 /= charge_icon ) THEN
                   speci(isp3)%icon = i
                   elements(i)%icon = isp3
                ELSE
                   WRITE(*,*) "Error - Problem in the determination"
                   WRITE(*,*) "        of conservation equations"
                   WRITE(*,*) "        => code stop"
                   STOP
                ENDIF
             ENDIF
             ! ELSE
             !    speci(isp1)%icon = i
             !    elements(i)%icon = isp1
             !    speci(elements(ii)%idsp_mn(2))%icon = ii
             !    elements(ii)%icon = isp2
             ! ENDIF
          ELSE
             IF ( isp3 /= charge_icon ) THEN
                speci(isp3)%icon = i
                elements(i)%icon = isp3
             ELSE
                WRITE(*,*) "Error - Problem in the determination"
                WRITE(*,*) "        of conservation equations"
                WRITE(*,*) "        => code stop"
                STOP
             ENDIF
          ENDIF
       ENDIF
    ENDDO

    ! --------------------------------------------
    ! find the index for conservation of PAHs
    ! --------------------------------------------
    IF ( speci(ind_PAH0)%density > speci(ind_PAHp)%density ) THEN
       id1 = ind_PAH0
       IF      ( speci(ind_PAHm)%density > speci(ind_PAH0)%density ) THEN
          id1 = ind_PAHm
          id2 = ind_PAH0
       ELSE IF ( speci(ind_PAHm)%density > speci(ind_PAHp)%density ) THEN
          id2 = ind_PAHm
       ELSE
          id2 = ind_PAHp
       ENDIF
    ELSE
       id1 = ind_PAHp
       IF      ( speci(ind_PAHm)%density > speci(ind_PAHp)%density ) THEN
          id1 = ind_PAHm
          id2 = ind_PAHp
       ELSE IF ( speci(ind_PAHm)%density > speci(ind_PAH0)%density ) THEN
          id2 = ind_PAHm
       ELSE
          id2 = ind_PAH0
       ENDIF
    ENDIF
    IF    ( id1 /= charge_icon ) THEN
       PAH_icon = id1
    ELSE
       PAH_icon = id2
    ENDIF

    ! --------------------------------------------
    ! Check that the species selected for 
    ! conservation are not used twice
    ! --------------------------------------------
    DO i = 1, Nelements
       DO ii = 1, Nelements
          IF ( i == ii ) CYCLE
          IF ( (elements(i)%icon == elements(ii)%icon).AND.(elements(i)%icon /= 0) ) THEN
             WRITE(*,*)
             WRITE(*,*) "Error - the same species is used twice for conservation equation"
             WRITE(*,*) "        find another way to select the species"
             WRITE(*,*) "        => code stop"
             WRITE(*,*)
             WRITE(*,*) elements(i)%name, elements(ii)%name
             WRITE(*,*) speci(elements(i)%icon)%name, speci(elements(ii)%icon)%name, speci(elements(ii)%idsp_mn(1))%name
             STOP
          ENDIF
       ENDDO
    ENDDO

    ! --------------------------------------------
    ! Test BG - 2016 - check the dominant carrier
    ! --------------------------------------------
    ! WRITE(*,*)
    ! DO i = 1, Nelements
    !    WRITE(*, '(A8,3X,A8,2X,2(ES10.3,2X),I4)') elements(i)%name, speci(elements(i)%icon)%name, &
    !                                           elements(i)%ab_init, elements(i)%ab_ref, elements(i)%icon
    ! ENDDO
    ! WRITE(*, '(A8,3X,A8,2X,(ES10.3,2X),I4)') "PAH     ", speci(PAH_icon)%name, PAH_abinit, PAH_icon
    ! --------------------------------------------

  END SUBROUTINE INDEX_CONSERVATION


  FUNCTION CHEMICAL_FORMULA (name) RESULT (formula)
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_SPECIES
    ! purpose :
    !    computes for one molecule its composition in elements H, C, ...
    ! subroutine/function needed :
    ! input variables :
    !    name -> molecule name
    ! ouput variables :
    ! results :
    !    formula -> vector containing the number of each element in the molecule
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(len = name_length), INTENT(in)   :: name
    INTEGER(KIND=LONG), DIMENSION(Nelements)   :: formula
    INTEGER(KIND=LONG)                         :: i, j, length
    INTEGER(KIND=LONG)                         :: digit1, digit2, number_element
    INTEGER(KIND=LONG)                         :: howmany, howmany_element, howmany_digit
    CHARACTER(len = 8), DIMENSION(name_length) :: table_name
    CHARACTER(len = 2)                         :: charact

    !-----------------------------------------------------
    ! initialization
    !-----------------------------------------------------
    table_name(:) = ''      ! describes what's in name(i)
    length = LEN_TRIM(name) ! length of name, without blanks
    formula = 0
    number_element = 0
    howmany = 0

    !-----------------------------------------------------
    ! filling the vector table_name
    !-----------------------------------------------------
    DO i = 1, length
       SELECT CASE (name(i:i))
       ! digit
       CASE ('0':'9')
          table_name(i) = "digit"
       ! 'He', 'Si' ...
       CASE ('a':'z')
          table_name(i) = "2letters"
        ! ion or on grain
       CASE ('+', '-', '*')
          table_name(i) = "iongrain"
       ! capital letters
       CASE ('A':'Z')
          table_name(i)="element"
       ! other
       CASE DEFAULT
          STOP "*** WARNING, incorrect character in CHEMICAL_FORMULA"
       END SELECT
    END DO

    !-----------------------------------------------------
    ! chemical formula
    !-----------------------------------------------------
    i = 1
    DO WHILE (i <= length)
       digit1 = 1
       digit2 = 0
       howmany_digit   = 1
       howmany_element = 1 ! 1 by default (if no digit)
       SELECT CASE ( TRIM(table_name(i)) )
       CASE ("element")
          ! determine how many characters to take into acount -> howmany
          IF (i < length) THEN
             IF (table_name(i+1) == "2letters") howmany_element = 2
          ENDIF
          howmany = howmany_element
          ! extracts from name
          charact = name(i:i+howmany-1)
          ! research of the element
          DO j = 1, Nelements
             IF (charact==TRIM(elements(j)%name)) EXIT
          ENDDO
          ! if not found ...
          IF (j == Nelements + 1) then
             PRINT *, 'name=', name
             PRINT *, 'charact=', charact
             STOP "*** WARNING, no chemical element for this molecule"
          endif
          IF ((j == Nelements) .AND. (charact /= TRIM(elements(j)%name))) THEN
             STOP "*** WARNING, no chemical element for this molecule"
          ENDIF
          number_element = j

       CASE ("digit")
          ! determine how many characters to take into acount -> howmany
          IF (i < length) THEN
             IF (table_name(i+1) == "digit") howmany_digit = 2
          ENDIF
          howmany = howmany_digit
          ! research of the number = digit1 or digit1*10+digit2
          digit1 = IACHAR(name(i:i)) - IACHAR('0')
          IF (howmany == 2) digit2 = IACHAR(name(i+1:i+1))-IACHAR('0')
       END SELECT

       ! determination of the number of elements in the molecule
       SELECT CASE (table_name(i))
       CASE('element')
          formula(number_element) = formula(number_element) + 1
       CASE('digit')
          formula(number_element) = formula(number_element) - 1
          IF (howmany_digit == 2) THEN 
             formula(number_element) = formula(number_element)+ digit1 * 10 + digit2
          ELSE
             formula(number_element) = formula(number_element) + digit1
          ENDIF
       CASE('iongrain') ! do nothing
       END SELECT

       ! look at the next characters in the name
       i = i + howmany
    END DO

  END FUNCTION CHEMICAL_FORMULA


  FUNCTION CHARGE_FORMULA (name) RESULT (charge)
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_SPECIES
    ! purpose :
    !    computes for one molecule its charge
    ! subroutine/function needed :
    ! input variables :
    !    name -> molecule name
    ! ouput variables :
    ! results :
    !    charge -> integer containing the number of charge of the species
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    CHARACTER(len = name_length), INTENT(in) :: name
    INTEGER(KIND=LONG)                       :: charge
    INTEGER(KIND=LONG)                       :: i, length

    length = LEN_TRIM(name)

    !-----------------------------------------------------
    ! filling the vector table_name
    !-----------------------------------------------------
    charge = 0
    DO i = length, 1, -1
       SELECT CASE (name(i:i))
       CASE ('+')
          charge = charge + 1
       CASE ('-')
          charge = charge - 1
       END SELECT
    END DO

  END FUNCTION CHARGE_FORMULA


  SUBROUTINE SPEID(espece, especeID, especeHN)
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_SPECIES
    ! purpose :
    !    transforms the name of a species in human readable name
    !    creates and ID without special characters.
    ! subroutine/function needed :
    ! input variables :
    !    name -> molecule name
    ! ouput variables :
    ! results :
    !    especeID and especeHN
    !---------------------------------------------------------------------------
  
    !TODO : ARRANGER DES CAS PARTICULIERS comme :
    !       Moments
    !       Electrons a renommer eventuellement mais attention aux effets de bords

    IMPLICIT NONE
    CHARACTER(len=name_length),   INTENT(IN)    :: espece        ! fortran name
    CHARACTER(len=lenSpeciHname), INTENT(OUT)   :: especeHN      ! human readable name
    CHARACTER(len=lenSpeciID),    INTENT(OUT)   :: especeID      ! ID

    CHARACTER(len=lenSpeciHname)                :: espece_sav
    INTEGER                                     :: pos_m         ! position of "-" in the string. Used for species as a-c3h4
    INTEGER                                     :: length_ini    ! Initial length of the string espece
    INTEGER                                     :: length        ! Length of the string to transform, <= length_ini
    INTEGER                                     :: i             ! Index in the string of espece

    !---------------------------------------------------
    ! Rappel:
    ! 65  --- A; 90  --- Z; 97  --- a; 122 --- z
    ! 43 = +, 45 = -, 35 = #, 58 = :
    !---------------------------------------------------

    especeHN = TRIM(ADJUSTL(espece))

    !--- Create ID (lower case) ------------------------
    ! At this step, we create the human readable names and le IDs

    !--- Look for "-" character to determine if the specy name has the syntax a-c3h4
    ! If so, we will only write in upper case the end of the string
    length_ini = LEN_TRIM(ADJUSTL(especeHN))    ! Initial length of the string
    pos_m = INDEX(especeHN,"-")                 ! Position of the "-" sign

    IF ((pos_m .NE. 0) .AND. (pos_m < length_ini)) THEN
       ! Species is of type a-c3h4
       espece_sav = especeHN                           ! Save the name of the species
       especeHN   = TRIM(ADJUSTL(especeHN(pos_m+1:)))  ! Get the part to transform
    ENDIF

    length = LEN_TRIM(ADJUSTL(especeHN))               ! New length of the string

     especeID = especeHN
     DO i = 1, length
        IF (     (IACHAR(especeID(i:i)) .GE. 65) &
           .AND. (IACHAR(especeID(i:i)) .LE. 90) ) THEN
           especeID(i:i) = ACHAR(IACHAR(especeID(i:i)) + 32 ) ! Transform to upper case
        ELSE IF (especeHN(i:i) .EQ. "+") THEN
           especeID(i:i) = "p"
        ELSE IF (especeHN(i:i) .EQ. "-") THEN
           especeID(i:i) = "m"
        ELSE IF (especeHN(i:i) .EQ. "#") THEN
           especeID(i:i) = "_"
        ELSE IF (especeHN(i:i) .EQ. ":") THEN
           especeID(i:i) = "a"
        ENDIF
     ENDDO
  
     ! If required add the prefix (a- for example) to species names
     IF ((pos_m .NE. 0) .AND. (pos_m < length_ini)) THEN
        especeHN = TRIM(ADJUSTL(espece_sav(1:pos_m)))//TRIM(ADJUSTL(especeHN))
        especeID = TRIM(ADJUSTL(espece_sav(1:pos_m-1)))//"_"//TRIM(ADJUSTL(especeID))
     ENDIF
  
     ! PRINT *,espece,"  ",especeHN,"  ",especeID
  
  END SUBROUTINE SPEID


  REAL(KIND=DP) FUNCTION xgas(ielem,nis)

    IMPLICIT none
    REAL(KIND=DP), DIMENSION(Nspec), INTENT(IN)  :: nis
    REAL(KIND=DP), DIMENSION(:),     ALLOCATABLE :: formulae
    INTEGER                                      :: j
    INTEGER                                      :: ielem
    ALLOCATE ( formulae(1:Nspec) )
    DO j = 1, Nspec
       formulae(j) = DBLE(speci(j)%formula(ielem))
    ENDDO
    xgas = ( sum(formulae(b_neu:e_neu)*nis(b_neu:e_neu))+ &
             sum(formulae(b_neg:e_neg)*nis(b_neg:e_neg))+ &
             sum(formulae(b_ion:e_ion)*nis(b_ion:e_ion)) ) - &
           ( formulae(ind_PAH0)*speci(ind_PAH0)%density + &
             formulae(ind_PAHp)*speci(ind_PAHp)%density + &
             formulae(ind_PAHm)*speci(ind_PAHm)%density  )
    DEALLOCATE(formulae)
    RETURN
  END FUNCTION xgas

  SUBROUTINE whereis(ielem,nis)
    USE module_phys_var, ONLY : nH
    IMPLICIT none
    REAL(KIND=DP), DIMENSION(Nspec), INTENT(IN)  :: nis
    REAL(KIND=DP), DIMENSION(:),     ALLOCATABLE :: formulae
    REAL(KIND=DP), DIMENSION(:),     ALLOCATABLE :: xgas
    REAL(KIND=DP)                                :: xxx
    INTEGER                                      :: i, j
    INTEGER                                      :: ielem

    ALLOCATE ( formulae(0:Nspec), xgas(1:Nspec) )
    formulae(0) = 0d0
    DO j = 1, Nspec
       formulae(j) = dble(speci(j)%formula(ielem))
    ENDDO
    xgas = 0d0
    xgas(b_neu:e_neu)= formulae(b_neu:e_neu)*nis(b_neu:e_neu)
    xgas(b_ion:e_ion)= formulae(b_ion:e_ion)*nis(b_ion:e_ion)
    xgas(b_neg:e_neg)= formulae(b_neg:e_neg)*nis(b_neg:e_neg)
    xgas(ind_PAH0)=0d0
    xgas(ind_PAHm)=0d0
    xgas(ind_PAHp)=0d0

    PRINT *, 'element ', elements(ielem)%name,' is in : (total=',sum(xgas)/nH,')'

    xxx = 0d0
    DO WHILE (any(xgas>0.0))
       i = maxloc(xgas,dim=1)
       xxx = xxx+xgas(i)/nH
       PRINT *, speci(i)%name, '=', xgas(i)/nH, xxx
       xgas(i) = 0.0d0
    ENDDO
    DEALLOCATE(formulae,xgas)
    RETURN
  END SUBROUTINE whereis

  SUBROUTINE info_elements(ifile)
    USE module_phys_var,only:nH
    IMPLICIT none
    INTEGER                                  :: i, j, ifile
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: formulae

    IF ( ifile == 0 ) ifile = 6 ! Stdout
    WRITE(ifile,*) '---- elements distribution ----'
    WRITE(ifile,*) 'element       Xtot         gas          PAH              Mantles      Cores '  
    ALLOCATE(formulae(0:Nspec))
    DO i=1,Nelements
       formulae(0)=0.0d0
       DO j=1,Nspec
          formulae(j)=dble(speci(j)%formula(i))
       ENDDO
       WRITE(ifile, "(a7,3x,5(1pe13.3),1x)") elements(i)%name, &
            sum(formulae(1:Nspec)*speci(1:Nspec)%density,dim=1)/nH, &
            (sum(formulae(b_neu:e_neu)*speci(b_neu:e_neu)%density)+ &
            sum(formulae(b_neg:e_neg)*speci(b_neg:e_neg)%density)+ &
            sum(formulae(b_ion:e_ion)*speci(b_ion:e_ion)%density) - &
            (formulae(ind_PAH0)*speci(ind_PAH0)%density+ &
            formulae(ind_PAHp)*speci(ind_PAHp)%density+ &
            formulae(ind_PAHm)*speci(ind_PAHm)%density) &
            )/nH, &
            (formulae(ind_PAH0)*speci(ind_PAH0)%density+ &
            formulae(ind_PAHp)*speci(ind_PAHp)%density+ &
            formulae(ind_PAHm)*speci(ind_PAHm)%density)/nH, &
            sum(formulae(b_gra:e_gra)*speci(b_gra:e_gra)%density)/nH, &
            sum(formulae(b_cor:e_cor)*speci(b_cor:e_cor)%density)/nH
    ENDDO
    DEALLOCATE(formulae)
  END SUBROUTINE info_elements

END MODULE MODULE_CHEMICAL_SPECIES
