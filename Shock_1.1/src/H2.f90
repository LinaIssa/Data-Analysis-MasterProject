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

MODULE MODULE_H2
  !*****************************************************************************
  !** The module 'MODULE_H2' contains variables and subroutines related to    **
  !** the H2 molecule, except op_H2 which is in MODULE_PHYS_VAR     **
  !**     * levels, ortho/para, reaction coefficients                         **
  !**     * Aij, quadrupole lines                                             **
  !** It contains also variables and subroutines needed to compute H2 cooling **
  !** and ortho/para ratio.                                                   **
  !*****************************************************************************
  USE MODULE_TECHCONFIG
  USE MODULE_CHEMICAL_SPECIES, ONLY : TYPE_UVDATA
  use tex
  IMPLICIT NONE
  INCLUDE "precision.f90"

  CHARACTER(len=lenfilename)    :: name_file_h2 = 'none'

  !-----------------------------------------
  ! energy levels
  ! variables defined in READ_H2_LEVELS
  !-----------------------------------------

  INTEGER                       :: NH2_lev_max ! maximal number of H2 level

  INTEGER(KIND=LONG)            :: Vmin_H2, Vmax_H2 ! min. and max. values for V (vibration)
  INTEGER(KIND=LONG)            :: Jmin_H2, Jmax_H2 ! min. and max. values for J (rotation)

  INTEGER(KIND=LONG)            :: Vmax_H2_eff      ! effective max. values for V
  INTEGER(KIND=LONG)            :: Jmax_H2_eff      ! effective max. values for J

  LOGICAL                       :: op_LTE = .false. ! .TRUE. if o/p H2 has been initialized to LTE value

  !-----------------------------------------
  ! ortho / para levels indices
  !-----------------------------------------
  INTEGER                            :: mlph2 ! number of even levels of H2 effectively used
  INTEGER                            :: mlih2 ! number of odd  levels of H2 effectively used
  INTEGER, ALLOCATABLE, DIMENSION(:) :: levpc ! list of even levels
  INTEGER, ALLOCATABLE, DIMENSION(:) :: levic ! list of odd  levels

  !-----------------------
  ! data type of one level
  !-----------------------
  TYPE TYPE_H2_LEVEL
     CHARACTER(len=lenIDLEV) :: Hname    ! Level human readable name
     CHARACTER(len=lenIDLEV) :: ID       ! Level ID
     INTEGER(KIND=LONG)      :: V, J     ! numbers of vibration and rotation
     REAL(KIND=DP)           :: Weight   ! statistical weight : 2J+1 (para) or 3(2J+1) (ortho)
     REAL(KIND=DP)           :: Energy   ! energy (K)
     REAL(KIND=DP)           :: Density  ! density of the level (cm-3)
     REAL(KIND=DP)           :: Density0 ! initial density of the level (cm-3) 
     REAL(KIND=DP)           :: Dens_old ! density at last call to DRIVE
     REAL(KIND=DP)           :: Col_dens ! column density of the level (cm-2)
     REAL(KIND=DP)           :: Form_gr  ! fraction of H2 formed on grains on level V J
     REAL(KIND=DP)           :: cd_l     ! left  side column density used for self-shielding (Doppler width weighted H2 column density)
     REAL(KIND=DP)           :: cd_l_old ! left  side column density used for self-shielding (Doppler width weighted H2 column density) at last call to DRIVE
     REAL(KIND=DP)           :: cd_r     ! right side column density used for self-shielding (Doppler width weighted H2 column density)
  END TYPE TYPE_H2_LEVEL

  ! initialized in READ_H2_LEVELS
  TYPE(TYPE_H2_LEVEL),DIMENSION(:),         ALLOCATABLE :: H2_lev      ! vector containing the H2 levels
  INTEGER(KIND=LONG), DIMENSION(:,:), SAVE, ALLOCATABLE :: index_VJ_H2 ! index (1..NH2_lev) of 1 level (V,J)

  INTEGER(KIND=LONG), PRIVATE, PARAMETER :: Num_undefined=-1 ! useful to initialize index_VJ_H2

  ! density of para and ortho-H2 (cm-3)
  REAL(KIND=DP) :: Dens_paraH2, Dens_orthoH2

  !--------------------------
  ! data type of one H2 line
  !--------------------------
  TYPE TYPE_H2_LINE
     CHARACTER(len=lenIDLIN) :: Hname     ! Line human readable name
     CHARACTER(len=lenIDLIN) :: ID        ! Line ID
     CHARACTER(LEN=9)        :: name      ! name of the line
     INTEGER(KIND=LONG)      :: Nup, Nlow ! index of the upper and the lower levels
     REAL(KIND=DP)           :: Aij       ! Aij (s-1)
     REAL(KIND=DP)           :: DeltaE    ! energy of the line : Eup-Elow (K)
     REAL(KIND=DP)           :: emiss     ! emiss of the line (erg/cm3/s)
     REAL(KIND=DP)           :: emiss_old ! emissivity at the last call to DRIVE
     REAL(KIND=DP)           :: intensity ! intensity integrated along the shock (erg/cm2/s/sr)
  END TYPE TYPE_H2_LINE

  !--- vector containing the H2 lines (dimension NH2_lines) ---
  INTEGER(KIND=LONG) :: NH2_lines      ! number of H2 lines
  TYPE (TYPE_H2_LINE), DIMENSION(:), SAVE, ALLOCATABLE :: H2_lines ! lines

  !---------------------------------------------
  ! collision rates
  ! read in READ_H2_RATES, used in EVOLUTION_H2
  !---------------------------------------------

  ! their dimension used are (4,NH2_lev,NH2_lev), all values outside are rejected
  REAL(KIND=DP), PRIVATE, DIMENSION(:,:,:), ALLOCATABLE :: rate_H_H2     ! collisions H-H2
  LOGICAL, PRIVATE, DIMENSION(:,:),         ALLOCATABLE :: mask_H_H2     ! collisions H-H2
  REAL(KIND=DP), PRIVATE, DIMENSION(:,:,:), ALLOCATABLE :: rate_He_H2    ! collisions He-H2
  REAL(KIND=DP), PRIVATE, DIMENSION(:,:,:), ALLOCATABLE :: rate_H2_H2    ! collisions H2-H2
  REAL(KIND=DP), PRIVATE, DIMENSION(25,0:16,0:36)       :: r_raw_GR_H2   ! collisions GR-H2 (raw rates)
  REAL(KIND=DP), PRIVATE, DIMENSION(25,0:16,0:36)       :: d2r_raw_GR_H2 ! collisions GR-H2 (raw rates)
  REAL(KIND=DP), PRIVATE, DIMENSION(0:16,0:36)          :: vin_GR_H2     ! corresponding velocity
  REAL(KIND=DP), PRIVATE, DIMENSION(0:16,0:36)          :: pgr0_GR_H2    ! corresponding rate

  !--------------------------------------------
  ! evolution terms for radiative transitions
  ! read in READ_H2_LINES
  !--------------------------------------------
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Aij_H2     ! Aij_H2(Nup,Nlow) : Aij (s-1) of the line 
                                                           ! Nup -> Nlow (N is the index of the level)
  REAL(KIND=DP), DIMENSION(:),   ALLOCATABLE :: sum_Aij_H2 ! sum of the Aij (s-1) of the lines starting 
                                                           ! from the level Nup

  !-------------------------------------------------
  ! Data from Gerlich for ortho-para transfer by H+
  ! with David Flower modification for ortho/para
  ! alternation (3/1).
  !-------------------------------------------------
  REAL(KIND=DP), PRIVATE, PARAMETER, DIMENSION(29) :: Gerlich = &
       (/2.024D-10, 9.575D-10, 2.729D-10, 7.910D-10 &
       , 2.246D-10, 6.223D-10, 1.828D-10, 4.975D-10, 1.484D-10 &
       , 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10 &
       , 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10 &
       , 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10 &
       , 1.000D-10, 3.000D-10, 1.000D-10, 3.000D-10, 1.000D-10 /)

  !--------------------------------------------------
  ! source terms for H2 level population (cm-3.s-1)
  ! and H2 internal energy (erg.cm-3.s-1)
  !--------------------------------------------------
  REAL(KIND=DP), public, parameter         :: H2_dissoc = 4.4781  ! H2 dissociation energy (in eV)
  REAL(KIND=DP), public                    :: H2_int_E            ! H2 internal energy at formation on grains (in K)
  REAL(KIND=DP)                            :: H2_energy=0.0_DP    ! calculated in COMPUTE_H2
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: YN_rovib_H2         ! change in H2 levels number density (cm-3.s-1)

  !--------------------------------------------------------------
  ! source terms of level selective chemial reactions involving H2
  ! 22 juin 2001 - JLB
  ! One type only yet: collisional dissociation of H2
  ! Allocated in : INITIALIZE_ROVIB_H2 (H2.f90)
  ! Used in : DIFFUN & CHEMISTRY
  !---------------------------------------------------------------
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_ch_H2        ! total selective reaction rate
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_ne_H2        ! total selective reaction rate (neutrals only)
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_io_H2        ! total selective reaction rate (ions only)
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_el_H2        ! total selective reaction rate (electrons only)
  REAL(KIND=DP), DIMENSION(:), SAVE, ALLOCATABLE :: Sel_rx_H2        ! partial selective reaction rate
  REAL(KIND=DP) :: Sel_tot_H2       ! Collisional dissociation of H2 (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: Sel_tne_H2       ! Collisional dissociation of H2 by neutrals (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: Sel_tio_H2       ! Collisional dissociation of H2 by ions (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: Sel_tel_H2       ! Collisional dissociation of H2 by electrons (summed over all levels) (cm-3 s-1)
  REAL(KIND=DP) :: For_gr_H2        ! Formation of H2 on grains (cm-3 s-1)

  !-----------------------------------
  ! cooling rate (erg.cm-3.s-1)
  ! due to H2 lines
  ! calculated in COMPUTE_OP_H2
  !-----------------------------------
  REAL(KIND=DP) :: cooling_H2   ! total cooling rate due to H2 lines (erg/cm3/s)
  REAL(KIND=DP) :: cooling_n_H2 ! cooling rate of the neutral fluid due to H2 excitation (erg/cm3/s)
  REAL(KIND=DP) :: cooling_i_H2 ! cooling rate of positively ionized fluid due to H2 excitation (erg/cm3/s)
  REAL(KIND=DP) :: cooling_e_H2 ! cooling rate of negatively ionized fluid due to H2 excitation (erg/cm3/s)

  !------------------------------------
  ! H2 electronic radiative transitions
  !------------------------------------
  CHARACTER(LEN=*), PARAMETER                     :: name_file_uvh2_lyman  = 'input/UVdata/uvh2b29.dat'
  CHARACTER(LEN=*), PARAMETER                     :: name_file_uvh2_werner = 'input/UVdata/uvh2c29.dat'

  TYPE(type_uvdata)                               :: uvh2_bx  ! H2 Lyman  transitions
  TYPE(type_uvdata)                               :: uvh2_cx  ! H2 Werner transitions

  REAL (KIND=dp),   DIMENSION(:,:,:), ALLOCATABLE :: aeinb    ! Prob of emission B-X
  REAL (KIND=dp),   DIMENSION(:,:,:), ALLOCATABLE :: aeinc    ! Prob of emission C-X
  REAL (KIND=dp),   DIMENSION(:,:),   ALLOCATABLE :: aeitb    ! Total prob of emission B-X to bound states
  REAL (KIND=dp),   DIMENSION(:,:,:), ALLOCATABLE :: aeitc    ! Total prob of emission B-X to bound states
  REAL (KIND=dp),   DIMENSION(:,:),   ALLOCATABLE :: tmunh2b  ! inverse lifetime of B levels
  REAL (KIND=dp),   DIMENSION(:,:),   ALLOCATABLE :: tmunh2c  ! inverse lifetime of C levels
  REAL (KIND=dp),   DIMENSION(:,:),   ALLOCATABLE :: cashm    ! H2 cascade coeficients
  REAL (KIND=dp),   DIMENSION(:,:,:), ALLOCATABLE :: bhmbx    ! H2 B-X cascade coeficients
  REAL (KIND=dp),   DIMENSION(:,:,:), ALLOCATABLE :: bhmcx    ! H2 C-X cascade coeficients

  INTEGER                                         :: Vmin_el  ! minimum Vup value for all electronic transitions
  INTEGER                                         :: Jmin_el  ! minimum Jup value for all electronic transitions
  INTEGER                                         :: Vmax_el  ! maximum Vup value for all electronic transitions
  INTEGER                                         :: Jmax_el  ! maximum Jup value for all electronic transitions
  INTEGER                                         :: Vmin_bel ! minimum Vup value for B-X transitions
  INTEGER                                         :: Jmin_bel ! minimum Jup value for B-X transitions
  INTEGER                                         :: Vmax_bel ! maximum Vup value for B-X transitions
  INTEGER                                         :: Jmax_bel ! maximum Jup value for B-X transitions
  INTEGER                                         :: Vmin_cel ! minimum Vup value for C-X transitions
  INTEGER                                         :: Jmin_cel ! minimum Jup value for C-X transitions
  INTEGER                                         :: Vmax_cel ! maximum Vup value for C-X transitions
  INTEGER                                         :: Jmax_cel ! maximum Jup value for C-X transitions

  INTEGER,          PARAMETER                     :: nvbc = 5 ! Highest vibrational level (B and C) for which pop is saved
  INTEGER,          PARAMETER                     :: njbc = 15! Highest rotational level (B and C) for which pop is saved
  REAL (KIND=dp),   DIMENSION(0:nvbc,0:njbc)      :: h2bnua   ! H2 local populations (B level)
  REAL (KIND=dp),   DIMENSION(0:nvbc,0:njbc)      :: h2cnua   ! H2 local populations (C level)
  REAL (KIND=dp)                                  :: h2elec   ! H2 local populations (sum of electronic levels up to nvbc & njbc)

  !----------------------------------------------
  ! evolution terms for pumping of electronic
  ! transitions initialized in READ_H2_ELECTRONIC
  ! and computed in FGKDATA
  !----------------------------------------------
  REAL (KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Pij_H2     ! Pij_H2(i,i) : excitation rate from i to j 
                                                            ! via pumping of an electronic state followed 
                                                            ! by radiative cascade
  REAL (KIND=DP), DIMENSION(:,:), ALLOCATABLE :: Pij_H2_old ! Pij_H2(i,i) at previous iteration
  REAL (KIND=DP), DIMENSION(:),   ALLOCATABLE :: sum_Pij_H2 ! sum of the Aij (s-1) of the lines starting 
                                                            ! from the level Nup
  REAL (KIND=DP), DIMENSION(:),   ALLOCATABLE :: pdh2vJ     ! pdh2vJ(lev) : H2 photodissociation rate from level lev
  REAL (KIND=DP)                              :: probdissH2 ! total H2 dissociation probability

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE READ_H2_LEVELS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads H2 levels : V, J, weight, energy.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    H2_lev, index_VJ_H2
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS,     ONLY  : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY  : Zero
    USE MODULE_PHYS_VAR,  ONLY  : NH2_lev
    IMPLICIT NONE

    CHARACTER(len=1)            :: charact
    INTEGER(KIND=LONG)          :: i, ii
    CHARACTER(LEN=*), PARAMETER :: name_file_H2_lev='input/H2_levels_Evj.in'
    INTEGER                     :: file_H2_lev

    INTEGER                     :: ios, whatever
    
    CHARACTER(LEN=5)            :: car_v, car_j

    ! Read the maximal number of H2 levels
    file_H2_lev = GET_FILE_NUMBER()
    OPEN(file_H2_lev,file=name_file_H2_lev,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    DO i=1,6
       READ(file_H2_lev,*)
    END DO
    NH2_lev_max=0
    DO
       READ (file_H2_lev, *, iostat=ios) whatever
       IF (ios /= 0) EXIT
       NH2_lev_max = NH2_lev_max + 1
    ENDDO
    CLOSE(file_H2_lev)

    ! initialization
    ALLOCATE (H2_lev(NH2_lev_max))
    H2_lev(1:NH2_lev_max)%Hname    = ""
    H2_lev(1:NH2_lev_max)%V        = 0
    H2_lev(1:NH2_lev_max)%J        = 0
    H2_lev(1:NH2_lev_max)%weight   = Zero
    H2_lev(1:NH2_lev_max)%Energy   = Zero
    H2_lev(1:NH2_lev_max)%density  = Zero
    H2_lev(1:NH2_lev_max)%density0 = Zero
    H2_lev(1:NH2_lev_max)%Dens_old = Zero
    H2_lev(1:NH2_lev_max)%Col_dens = Zero
    H2_lev(1:NH2_lev_max)%Form_gr  = Zero

    ! file opening
    file_H2_lev = GET_FILE_NUMBER()
    OPEN(file_H2_lev,file=name_file_H2_lev,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! comments
    DO i=1,6
       READ(file_H2_lev,'(A1)')charact
    END DO

    DO i=1, NH2_lev_max
       READ(file_H2_lev,*) &
            ii, &
            H2_lev(i)%V, &
            H2_lev(i)%J, &
            H2_lev(i)%weight, &
            H2_lev(i)%Energy
            WRITE(car_v,'(I5)') H2_lev(i)%V
            WRITE(car_j,'(I5)') H2_lev(i)%J
            H2_lev(i)%Hname = "v=" // TRIM(ADJUSTL(car_v)) // "," // "J=" // TRIM(ADJUSTL(car_j))
            H2_lev(i)%ID = TRIM(ADJUSTL(car_v)) // "_" // TRIM(ADJUSTL(car_j))
       ! check if all lines have been correctly read, we must have i=ii
       IF (i /= ii) STOP "*** WARNING, error de read des niveaux de H2"
    END DO

    ! min. and max. values of V and J for these levels
    Vmin_H2=MINVAL(H2_lev(1:NH2_lev_max)%V)
    Vmax_H2=MAXVAL(H2_lev(1:NH2_lev_max)%V)
    Jmin_H2=MINVAL(H2_lev(1:NH2_lev_max)%J)
    Jmax_H2=MAXVAL(H2_lev(1:NH2_lev_max)%J)

    ! allocation and filling of the table index_VJ_H2
    ALLOCATE(index_VJ_H2(Vmin_H2:Vmax_H2, Jmin_H2:Jmax_H2))
    index_VJ_H2(Vmin_H2:Vmax_H2,Jmin_H2:Jmax_H2) = num_undefined
    DO i=1, NH2_lev_max
       index_VJ_H2(H2_lev(i)%V,H2_lev(i)%J) = i
    ENDDO

    ! file closure
    CLOSE(file_H2_lev)

    ! Compute the effective maximal value of V
    Vmax_H2_eff = MAXVAL(H2_lev(1:NH2_lev)%V)
    Jmax_H2_eff = MAXVAL(H2_lev(1:NH2_lev)%J)

  END SUBROUTINE READ_H2_LEVELS


  SUBROUTINE LLEV_H2
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    compute the number of para & ortho levels and their respective indices
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    mlph2, levpc, mlih2, levic
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY  : NH2_lev

    IMPLICIT NONE
    INTEGER :: i

    mlph2 = COUNT( MOD(H2_lev(1:NH2_lev)%J,2) == 0)
    ALLOCATE (levpc(mlph2))
    levpc = PACK( (/(i, i=1,NH2_lev)/), mask = (MOD(H2_lev(1:NH2_lev)%J,2) == 0) )
    mlih2 = COUNT( MOD(H2_lev(1:NH2_lev)%J,2) == 1)
    ALLOCATE (levic(mlih2))
    levic = PACK( (/(i, i=1,NH2_lev)/), mask = (MOD(H2_lev(1:NH2_lev)%J,2) == 1) )

  END SUBROUTINE LLEV_H2


  SUBROUTINE INITIALIZE_ROVIB_H2
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    initialize population of the ro-vibrational levels of H2
    !    3 possibilities :
    !          * if input file (h2exfile) -> h2* read. else
    !          * op_H2 = op_H2 read in READ_PARAMETERS or
    !          * op_H2 = op LTE at temperature Tn (if op_H2 > op_H2_LTE)
    !    note : op_H2 is (re)-calculated
    ! subroutine/function needed :
    !    COMPUTE_OP_H2
    ! input variables :
    ! ouput variables :
    ! results :
    !    H2_lev, op_H2
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS,            ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR,         ONLY : Tn, op_H2, op_H2_in, NH2_lev, iforH2
    USE MODULE_CHEMICAL_SPECIES, ONLY : speci, ind_H2
    USE MODULE_CONSTANTS,        ONLY : Zero, kB, EVerg
    USE MODULE_TECHCONFIG,       ONLY : h2exfile

    IMPLICIT NONE

    CHARACTER(len=80)       :: dummy_line
    INTEGER                 :: ndummy
    INTEGER(KIND=LONG)      :: file_h2
    INTEGER                 :: NH2_lev_file

    INTEGER(KIND=LONG)      :: i
    REAL(KIND=DP)           :: Zortho, Zpara, weight_ortho, weight_para
    REAL(KIND=DP),PARAMETER :: op_H2_LTE=999._DP
    REAL(KIND=DP),PARAMETER :: population_limit=1.D-20 ! lower limit for H2 population
    REAL(KIND=DP)           :: int_enrg, T_ini, tform, tf_min, tf_max, res

    if (.not.allocated(YN_rovib_H2)) then
       ALLOCATE (YN_rovib_H2(NH2_lev))
       ALLOCATE (Sel_ch_H2(NH2_lev))
       ALLOCATE (Sel_ne_H2(NH2_lev))
       ALLOCATE (Sel_io_H2(NH2_lev))
       ALLOCATE (Sel_el_H2(NH2_lev))
       ALLOCATE (Sel_rx_H2(NH2_lev))
    endif

    ! ==========================================================================
    ! Initialization of H2 level populations
    ! ==========================================================================

    !-----------------------------------------------------------
    ! Use values stored in an input file
    !-----------------------------------------------------------
    IF (h2exfile /= 'none') THEN

       !--- opening file ---
       file_h2 = GET_FILE_NUMBER()
       name_file_h2 = TRIM(data_dir) // TRIM(h2exfile)
       OPEN(file_h2, file = name_file_h2, status='OLD',&
         access = 'SEQUENTIAL', form = 'FORMATTED', action = 'READ')

       !--- Count number of comment lines ---
       dummy_line = '!'
       ndummy     = -1
       DO WHILE (dummy_line(1:1) == '!')
          ndummy = ndummy+1
          READ(file_h2,'(A80)') dummy_line
       ENDDO
       REWIND(file_h2)

       !--- Skip comments, read nb levels ---
       DO i = 1, ndummy - 2
          READ(file_h2,*)
       END DO
       READ(file_h2,'(15X,I8)') NH2_lev_file
       READ(file_h2,*)
       
       !--- Read h2 level populations ---
       DO i = 1, MIN(NH2_lev,NH2_lev_file)
          READ(file_h2,'(12X,ES11.3E3)') H2_lev(i)%density
       ENDDO

       !--- Set all missing levels to 1e-20 ---
       H2_lev(NH2_lev_file:NH2_lev)%density = population_limit

       !--- debug : check initial populations ---
       !DO i = 1, NH2_lev
       !   WRITE(*,*) H2_lev(i)%density
       !ENDDO
       !STOP

       !--- Close file ---
       CLOSE(file_h2)

    !-----------------------------------------------------------
    ! Or a given O/P ratio
    !-----------------------------------------------------------
    ELSE
       ! which o/p to take ?
       IF (op_H2_in > op_H2_LTE) THEN
          op_LTE=.TRUE.
       ELSE
          op_LTE=.FALSE.
       END IF
   
       ! computes the populations, given the o/p choice
       IF (OP_LTE) THEN
          !-----------------------------------------------------
          ! the o/p has it's LTE value at the temperature Tn ---
          !-----------------------------------------------------
          DO i=1,NH2_lev
             H2_lev(i)%density=H2_lev(i)%weight * &
                  TEXP(-(H2_lev(i)%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/Tn)
          ENDDO
       ELSE
          op_H2=op_H2_in
          !--------------------------------------------------
          ! the o/p has the value read in READ_PARAMETERS ---
          !--------------------------------------------------
          ! partition functions and weights for ortho-H2 and para-H2
          Zortho=Zero
          T_ini = Tn
          DO i=index_VJ_H2(0,1),NH2_lev,2 ! first ortho level : (V=0,J=1)
             Zortho = Zortho + H2_lev(i)%weight * &
                  TEXP(-(H2_lev(i)%energy-H2_lev(index_VJ_H2(0,0))%Energy)/Tn)
          END DO
          weight_ortho=op_H2/(1._DP+op_H2)/Zortho
   
          Zpara=Zero
          DO i=index_VJ_H2(0,0),NH2_lev,2 ! first para level : (V=0,J=0)
             Zpara = Zpara + H2_lev(i)%weight * &
                  TEXP(-(H2_lev(i)%energy-H2_lev(index_VJ_H2(0,0))%Energy)/Tn)
          END DO
          weight_para=1._DP/(1._DP+op_H2)/Zpara
   
          ! populations for ortho levels
          WHERE (MOD(H2_lev(1:NH2_lev)%J,2) > 0)
             H2_lev(1:NH2_lev)%density = weight_ortho * H2_lev(1:NH2_lev)%weight * &
                  EXP( -(H2_lev(1:NH2_lev)%energy - H2_lev(index_VJ_H2(0,0))%energy) / Tn )
          END WHERE
   
          ! populations for para levels
          WHERE (MOD(H2_lev(1:NH2_lev)%J,2) == 0)
             H2_lev(1:NH2_lev)%density = weight_para * H2_lev(1:NH2_lev)%weight * &
                  EXP( -(H2_lev(1:NH2_lev)%energy - H2_lev(index_VJ_H2(0,0))%energy) / Tn )
          END WHERE
       END IF

    ENDIF

    ! avoids too small numbers => set lower limit to the populations
    WHERE (H2_lev(1:NH2_lev)%density < population_limit)
       H2_lev(1:NH2_lev)%density=population_limit
    END WHERE

    H2_lev(1:NH2_lev)%density0 = H2_lev(1:NH2_lev)%density

    !   print *, H2_lev(1:NH2_lev)%density

    !  Formation of H2 on grains :
    !  -1: formation in the v,J = 0,0 and 0,1 levels only
    !   0: 1/3 of 4.4781 eV in internal energy (=> 17249 K) (Allen, 1999)
    !   1: Proportional to Boltzman distribution at 17249 K
    !   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)
    !   3: v = 6, J = 0,1
    !   4: fraction = relative populations at t, initialised as H2_lev(1:NH2_lev)%density
    !                 and changed during integration

    IF (iforH2 == -1) then

      H2_lev(index_VJ_H2(0,0))%Form_gr = 0.25_DP
      H2_lev(index_VJ_H2(0,1))%Form_gr = 0.75_DP

    ELSE IF (iforH2 == 0) then

      int_enrg = H2_dissoc * EVerg / (3.0_DP * kB)
      tform = int_enrg
      tf_min = 0.0_dp
      tf_max = 3.0_DP * int_enrg
      res = 1.0_dp

      ! find temperature for which SUM(H2_lev(1:NH2_lev)%Energy * weight * exp(-H2_lev(1:NH2_lev)%energy / T))
      ! is equal to 1/3 x H2_dissoc
      do while (abs(res) > 1.0e-3_DP)
        res = SUM(H2_lev(1:NH2_lev)%weight * (int_enrg-H2_lev(1:NH2_lev)%Energy) &
            * EXP(-H2_lev(1:NH2_lev)%Energy/tform))
        if (res < 0.0_dp) then
          tf_max = tform
          tform = 0.5_DP * (tf_min + tf_max)
        else
          tf_min = tform
          tform = 0.5_DP * (tf_min + tf_max)
        endif
      end do

      int_enrg = tform

      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%weight * &
           EXP(-(H2_lev(1:NH2_lev)%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/int_enrg)
      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%Form_gr / SUM(H2_lev(1:NH2_lev)%Form_gr)

    ELSE IF (iforH2 == 1) then

      int_enrg = H2_dissoc * EVerg / (3.0_DP * kB)
      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%weight * &
           EXP(-(H2_lev(1:NH2_lev)%Energy-H2_lev(index_VJ_H2(0,0))%Energy)/int_enrg)
      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%Form_gr / SUM(H2_lev(1:NH2_lev)%Form_gr)

    ELSE IF (iforH2 == 2) then

      H2_lev(index_VJ_H2(14,0))%Form_gr = 0.25_DP
      H2_lev(index_VJ_H2(14,1))%Form_gr = 0.75_DP

    ELSE IF (iforH2 == 3) then

      H2_lev(index_VJ_H2(6,0))%Form_gr = 0.25_DP
      H2_lev(index_VJ_H2(6,1))%Form_gr = 0.75_DP

    ELSE IF (iforH2 == 4) then

      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%density
      H2_lev(1:NH2_lev)%Form_gr = H2_lev(1:NH2_lev)%Form_gr / SUM(H2_lev(1:NH2_lev)%Form_gr)

    ENDIF

    H2_int_E = SUM(H2_lev(1:NH2_lev)%Form_gr * H2_lev(1:NH2_lev)%Energy)

    ! normalization at H2 density
    H2_lev(1:NH2_lev)%density = H2_lev(1:NH2_lev)%density &
                              * speci(ind_H2)%density / SUM(DBLE(H2_lev(1:NH2_lev)%density))

    ! op_H2 is (re)-calculated for checking
    CALL COMPUTE_OP_H2
    IF (op_LTE) WRITE(*,'("--- LTE chosen for ortho:para-H2, &
         &o/p = ",ES7.1," ---")')op_H2


  END SUBROUTINE INITIALIZE_ROVIB_H2



  SUBROUTINE READ_H2_RATES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    read collision rates for H-H2, He-H2 and H2-H2.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    rate_H_H2, rate_He_H2, rate_H2_H2
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    USE MODULE_PHYS_VAR, ONLY : NH2_lev, H_H2_flag
    USE NUM_REC, ONLY: spline
    IMPLICIT NONE

    INTEGER(KIND=LONG)          :: Nrate, i, LevU, LevL
    CHARACTER(LEN=25)           :: name_file_H_H2
    INTEGER                     :: file_H_H2
    CHARACTER(LEN=*), PARAMETER :: name_file_He_H2='input/coeff_He_H2.in'
    INTEGER                     :: file_He_H2
    CHARACTER(LEN=*), PARAMETER :: name_file_H2_H2='input/coeff_H2_H2.in'
    INTEGER                     :: file_H2_H2
    ! CHARACTER(LEN=*), PARAMETER :: name_file_GR_H2='input/coeff_GR_H2.in'
    INTEGER                     :: file_GR_H2
    integer :: ii, jj, kk
    integer :: iju, ivu
    REAL(KIND=DP) :: yp1, ypn
    REAL(KIND=DP) :: toto1, toto2, toto3, toto4
    REAL(KIND=DP), DIMENSION(1:25) :: vin
    data yp1, ypn / 1.0d30, 1.0d30 /
    data vin / 2.d0,4.d0,6.d0,8.d0,&
               10.d0,12.d0,14.d0,16.d0,18.d0,&
               20.d0,22.d0,24.d0,26.d0,28.d0,&
               30.d0,32.d0,34.d0,36.d0,38.d0,&
               40.d0,42.d0,44.d0,46.d0,48.d0,50.d0/


    ! initialization
    ALLOCATE (rate_H_H2(4,NH2_lev,NH2_lev))
    ALLOCATE (mask_H_H2(NH2_lev,NH2_lev))
    ALLOCATE (rate_He_H2(4,NH2_lev,NH2_lev))
    ALLOCATE (rate_H2_H2(4,NH2_lev,NH2_lev))
    rate_H_H2(1,:,:)=-50._DP
    rate_H_H2(2:4,:,:)=Zero
    mask_H_H2(:,:)= .FALSE.
    rate_He_H2(1,:,:)=-50._DP
    rate_He_H2(2:4,:,:)=Zero
    rate_H2_H2(1,:,:)=-50._DP
    rate_H2_H2(2:4,:,:)=Zero
    r_raw_GR_H2(:,:,:)=Zero
    d2r_raw_GR_H2(:,:,:)=Zero
    vin_GR_H2(:,:)=0

    !-----------------------
    ! rate H-H2
    !-----------------------

    ! We read either DRF file or MM file or BOTH
    ! If BOTH, then first MM then DRF so that Quantum rates override Semiclassical ones

    if (H_H2_flag == "MM" .OR. H_H2_flag == "BOTH") then
      name_file_H_H2 = 'input/coeff_H_H2_MM.in'
      file_H_H2 = GET_FILE_NUMBER()
      OPEN(file_H_H2,file=name_file_H_H2,status='OLD',&
           access='SEQUENTIAL',form='FORMATTED',action='READ')
      DO i=1,4
         READ(file_H_H2,*)
      ENDDO

      READ(file_H_H2,*)Nrate
      DO i=1, Nrate
         READ(file_H_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
         IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
           rate_H_H2(1,LevU,LevL) = toto1
           rate_H_H2(2,LevU,LevL) = toto2
           rate_H_H2(3,LevU,LevL) = toto3
           rate_H_H2(4,LevU,LevL) = toto4
         ENDIF
      END DO
      CLOSE(file_H_H2)
    endif

    if (H_H2_flag == "DRF" .OR. H_H2_flag == "BOTH") then
      name_file_H_H2 = 'input/coeff_H_H2_DRF.in'
      file_H_H2 = GET_FILE_NUMBER()
      OPEN(file_H_H2,file=name_file_H_H2,status='OLD',&
           access='SEQUENTIAL',form='FORMATTED',action='READ')
      DO i=1,4
         READ(file_H_H2,*)
      ENDDO

      READ(file_H_H2,*)Nrate
      DO i=1, Nrate
         READ(file_H_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
         IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
           rate_H_H2(1,LevU,LevL) = toto1
           rate_H_H2(2,LevU,LevL) = toto2
           rate_H_H2(3,LevU,LevL) = toto3
           rate_H_H2(4,LevU,LevL) = toto4
           mask_H_H2(LevU,LevL) = .TRUE.
         ENDIF
      END DO
      CLOSE(file_H_H2)
    endif

    !-------------
    ! rate He-H2
    !-------------
    ! file opening
    file_He_H2 = GET_FILE_NUMBER()
    OPEN(file_He_H2,file=name_file_He_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! commentss
    DO i=1,4
       READ(file_He_H2,*)
    ENDDO

    ! number of collision rates
    READ(file_He_H2,*)Nrate
    ! read des coefficients ligne par ligne
    DO i=1, Nrate
       READ(file_He_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
       IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
         rate_He_H2(1,LevU,LevL) = toto1
         rate_He_H2(2,LevU,LevL) = toto2
         rate_He_H2(3,LevU,LevL) = toto3
         rate_He_H2(4,LevU,LevL) = toto4
       ENDIF
    END DO
    ! file closure
    CLOSE(file_He_H2)

    !-------------
    ! rate H2-H2
    !-------------
    ! file opening
    file_H2_H2 = GET_FILE_NUMBER()
    OPEN(file_H2_H2,file=name_file_H2_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    ! comments
    DO i=1,4
       READ(file_H2_H2,*)
    ENDDO

    ! number of collision rates
    READ(file_H2_H2,*)Nrate
    ! read des coefficients ligne par ligne
    DO i=1, Nrate
       READ(file_H2_H2,*) LevU, LevL, toto1, toto2, toto3, toto4
       IF ((LevU <= NH2_lev) .AND. (LevL <= NH2_lev)) THEN
         rate_H2_H2(1,LevU,LevL) = toto1
         rate_H2_H2(2,LevU,LevL) = toto2
         rate_H2_H2(3,LevU,LevL) = toto3
         rate_H2_H2(4,LevU,LevL) = toto4
       ENDIF
    END DO
    ! file closure
    CLOSE(file_H2_H2)

    ! Raw rates - file opening Para, then Ortho
    file_GR_H2 = GET_FILE_NUMBER()
    OPEN(file_GR_H2,file="input/PH2GR.DAT",status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    DO ii=1,25
       READ(file_GR_H2,*)
       DO kk=0,36,2
         READ(file_GR_H2,*) iju, (r_raw_GR_H2(ii,jj,kk),jj=0,16)
       END DO
    END DO
    ! file closure
    CLOSE(file_GR_H2)

    file_GR_H2 = GET_FILE_NUMBER()
    OPEN(file_GR_H2,file="input/OH2GR.DAT",status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    DO ii=1,25
       READ(file_GR_H2,*)
       DO kk=1,35,2
         READ(file_GR_H2,*) iju, (r_raw_GR_H2(ii,jj,kk),jj=0,16)
       END DO
    END DO
    ! file closure
    CLOSE(file_GR_H2)

    DO ii=3,NH2_lev
      iju = H2_lev(ii)%J
      ivu = H2_lev(ii)%V
      DO jj=1,25
        if (r_raw_GR_H2(jj,ivu,iju) > 0.0_DP) then
          vin_GR_H2(ivu,iju) = vin(jj)
          pgr0_GR_H2(ivu,iju) = r_raw_GR_H2(jj,ivu,iju)
          exit
        endif
      END DO
      call spline (vin(jj:),r_raw_GR_H2(jj:,ivu,iju),yp1,ypn,d2r_raw_GR_H2(jj:,ivu,iju))
    END DO

  END SUBROUTINE READ_H2_RATES


  SUBROUTINE READ_H2_LINES
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads H2 lines (levels, energies, Aij).
    ! subroutine/function needed :
    !    NAME_H2_LINE
    ! input variables :
    ! ouput variables :
    ! results :
    !    Aij_H2, sum_Aij_H2
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS,    ONLY   : GET_FILE_NUMBER
    USE MODULE_CONSTANTS
    USE MODULE_PHYS_VAR, ONLY   : NH2_lev, NH2_lines_out

    IMPLICIT NONE
    REAL(KIND=DP)               :: Aij_S, Aij_Q, Aij_O
    INTEGER(KIND=LONG)          :: i, Vup, Vlow, Jup, Jlow, Nup, Nlow, line_number
    CHARACTER(LEN=1)            :: line_type
    CHARACTER(LEN=*), PARAMETER :: name_file_Aij_H2='input/Aij_H2.in'
    CHARACTER(LEN=*), PARAMETER :: format_Aij_H2='(3I4,3D15.6)'
    INTEGER                     :: file_Aij_H2

    INTEGER                     :: nblines
    INTEGER                     :: ios, whatever

    !----------------------------------------------------
    ! Read the number of H2 lines
    !----------------------------------------------------
    file_Aij_H2 = GET_FILE_NUMBER()
    OPEN(file_Aij_H2,file=name_file_Aij_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    DO i=1,5
       READ(file_Aij_H2,*)
    END DO
    nblines = 0
    DO
       READ (file_Aij_H2, *, iostat=ios) whatever
       IF (ios /= 0) EXIT
       nblines = nblines + 1
    ENDDO
    CLOSE(file_Aij_H2)

    !----------------------------------------------------
    ! initialization
    !----------------------------------------------------
    ALLOCATE ( Aij_H2(NH2_lev_max,NH2_lev_max) )
    ALLOCATE ( sum_Aij_H2(NH2_lev_max)         )
    Aij_H2(1:NH2_lev_max,1:NH2_lev_max) = Zero
    sum_Aij_H2(1:NH2_lev_max)           = Zero

    !----------------------------------------------------
    ! file opening
    !----------------------------------------------------
    file_Aij_H2 = GET_FILE_NUMBER()
    OPEN(file_Aij_H2,file=name_file_Aij_H2,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    DO i=1,5
       READ(file_Aij_H2,*)
    END DO

    ! read line per line. There is max. 3 Aij per level (lines S, Q, O)
    ! the file is in increasing Vup order
    DO i = 1, nblines
       READ(file_Aij_H2,format_Aij_H2) Vup, Vlow, Jup, Aij_S, Aij_Q, Aij_O
       IF ( (Jup > Jmax_H2).OR.(Vup > Vmax_H2) ) CYCLE
       Nup = index_VJ_H2(Vup,Jup)
       ! ------
       ! S line
       ! ------
       Jlow = Jup - 2
       IF ( (Jlow >= 0) .AND. (Aij_S /= Zero) ) THEN
          Nlow = index_VJ_H2(Vlow,Jlow)
          Aij_H2(Nup,Nlow) = Aij_S
       ENDIF
       ! ------
       ! Q line
       ! ------
        Jlow = Jup
        IF (Aij_Q /= Zero) THEN
           Nlow = index_VJ_H2(Vlow,Jlow)
           Aij_H2(Nup,Nlow) = Aij_Q
        ENDIF
       ! ------
       ! O line
       ! ------
       Jlow = Jup + 2
       IF ( (Jlow <= Jmax_H2) .AND. (Aij_O /= Zero) ) THEN
          Nlow = index_VJ_H2(Vlow,Jlow)
          Aij_H2(Nup,Nlow) = Aij_O
       ENDIF
    END DO

    ! computes the sum of the Aij of the lines
    DO Nup = 1, NH2_lev_max
       sum_Aij_H2(Nup)=SUM(Aij_H2(Nup,1:Nup-1))
    END DO

    ! file closure
    CLOSE(file_Aij_H2)

    ! counts the number of transitions in the model (Aij_H2 > 0)
    ! involving level whose population are computed properly
    NH2_lines = 0
    DO Nup = 1, NH2_lev
       DO Nlow = 1, NH2_lev
          IF ( Aij_H2(Nup,Nlow) > Zero ) THEN
             NH2_lines = NH2_lines + 1
          ENDIF
       ENDDO
    ENDDO
    NH2_lines_out = min(NH2_lines,NH2_lines_out)
    ! allocation and initialization of H2_lines
    ALLOCATE(H2_lines(NH2_lines))
    H2_lines(1:NH2_lines)%name      = ''
    H2_lines(1:NH2_lines)%Nup       = 0
    H2_lines(1:NH2_lines)%Nlow      = 0
    H2_lines(1:NH2_lines)%emiss     = Zero
    H2_lines(1:NH2_lines)%emiss_old = Zero
    H2_lines(1:NH2_lines)%intensity = Zero
    H2_lines(1:NH2_lines)%DeltaE    = Zero
    H2_lines(1:NH2_lines)%Aij       = Zero

    ! vector H2_lines is filled in order of increasing Eup
    line_number = 0
    IF (NH2_lines < 1) RETURN
    DO Nup = index_VJ_H2(0,2), NH2_lev
       DO Nlow = 1, Nup-1
          ! tests if line exists
          IF (Aij_H2(Nup,Nlow) > Zero) THEN
             line_number =line_number+1
             Vup  = H2_lev(Nup)%V
             Jup  = H2_lev(Nup)%J
             Vlow = H2_lev(Nlow)%V
             Jlow = H2_lev(Nlow)%J
             ! finds the type of the line (S, Q, O)
             SELECT CASE(Jup-Jlow)
             CASE (2)
                line_type='S'
             CASE (0)
                line_type='Q'
             CASE (-2)
                line_type='O'
             CASE DEFAULT
             END SELECT
             ! fills H2_lines
             H2_lines(line_number)%name   = NAME_H2_LINE(Vup=Vup,Vlow=Vlow,Jlow=Jlow,line_type=line_type)
             H2_lines(line_number)%Nup    = Nup
             H2_lines(line_number)%Nlow   = Nlow
             H2_lines(line_number)%Aij    = Aij_H2(Nup,Nlow)
             H2_lines(line_number)%DeltaE = H2_lev(Nup)%energy-H2_lev(Nlow)%energy
             H2_lines(line_number)%Hname  = TRIM(ADJUSTL(H2_lev(Nup)%Hname)) &
                                          // "->" // TRIM(ADJUSTL(H2_lev(Nlow)%Hname))
             H2_lines(line_number)%ID     = TRIM(ADJUSTL(H2_lev(Nup)%ID)) &
                                          // "_" // TRIM(ADJUSTL(H2_lev(Nlow)%ID))
             ! WRITE(*,'(4(I4,2X),2(ES14.7,2X))') H2_lev(Nup)%V, H2_lev(Nup)%J, H2_lev(Nlow)%V, H2_lev(Nlow)%J, &
             !                        H2_lines(line_number)%DeltaE, hp*clum / (H2_lines(line_number)%DeltaE * kb) *1e4_dp
          ENDIF
       END DO
    END DO

  END SUBROUTINE READ_H2_LINES


  SUBROUTINE READ_H2_ELECTRONIC
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads H2 electronic radiative transitions files
    !       - uvh2b29.dat: Lyman transitions  (all v) J <= 29
    !       - uvh2c29.dat: Werner transitions (all v) J <= 29 (01 janvier 1989)
    !         data's from January the 1st 1989
    !         completed to J = 29 some times during the 90's
    !         One transition per line.
    ! subroutine/function needed :
    !    CACOEF
    ! input variables :
    ! ouput variables :
    ! results :
    !    uvh2_(B & C)X%vl ..... v lower level
    !    uvh2_(B & C)X%jl ..... J lower level
    !    uvh2_(B & C)X%vu ..... v upper level
    !    uvh2_(B & C)X%dj ..... Delta J (J_up = J_low + Delta J)
    !    uvh2_(B & C)X%iw ..... index of line in wavelength grid
    !    uvh2_(B & C)X%osc .... oscillator strength
    !    uvh2_(B & C)X%wlg .... wavelength (Angstrom)
    !    uvh2_(B & C)X%ilft ... inverse life time (s-1)
    !    uvh2_(B & C)X%pbdi ... Upper level dissociation probability
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY  : pi, qe, me, clum
    USE MODULE_PHYS_VAR, ONLY   : NH2_lev
    USE MODULE_TECHCONFIG, ONLY : fichier
    USE MODULE_TOOLS, ONLY      : GET_FILE_NUMBER

    IMPLICIT NONE
    INTEGER            :: nuv     ! number of radiative transitions in input file
    INTEGER            :: iread   ! file opening number
    INTEGER            :: nvl     ! lower v level
    INTEGER            :: njl     ! lower J level
    INTEGER            :: nvu     ! upper v level
    INTEGER            :: nju     ! upper J level
    INTEGER            :: ndj     ! delta J = J_up - J_lw
    REAL (KIND=dp)     :: uv1     ! oscillator strength
    REAL (KIND=dp)     :: uv2     ! wavelength (Angstrom)
    REAL (KIND=dp)     :: uv3     ! inverse life time (s-1)
    REAL (KIND=dp)     :: uv4     ! Upper level dissociation probability

    REAL (KIND=dp)     :: coefosc ! conversion coefficient
    REAL (KIND=dp)     :: dumy    ! dummy real

    INTEGER            :: lev     ! dummy integer
    INTEGER            :: i       ! dummy integer
    INTEGER            :: j       ! dummy integer
    INTEGER            :: ii      ! dummy integer

    ! Numerical coeficient : 8 pi**2 e**2 / (me * c) x 1e16 (for Angstrom)
    coefosc = 8.0e16_dp * pi*pi * qe**2 / (me * clum)

    !----------------------------------------------------
    ! Prepare H2 electronic transitions
    ! 1 - Determine the number of lines in Lyman & Werner
    !     transitions files
    ! 2 - Find Vmin_bel, Jmin_bel, Vmax_bel, Jmax_bel
    !          Vmin_cel, Jmin_cel, Vmax_cel, Jmax_cel
    !          Vmin_el, Jmin_el, Vmax_el, & Jmax_el
    !----------------------------------------------------
    Vmin_bel    = 0
    Jmin_bel    = 0
    Vmax_bel    = 0
    Jmax_bel    = 0
    uvh2_bx%nuv = 0
    fichier = TRIM(name_file_uvh2_lyman)
    iread = GET_FILE_NUMBER()
    OPEN (iread, FILE = fichier, STATUS = 'old')
    READ (iread,"(i8)") nuv
    DO i = 1, nuv
       READ (iread,*) ii, ii, nvl, njl, nvu, ndj, uv1, uv2, uv3, uv4
       lev = index_VJ_H2(nvl,njl)
       IF (lev > NH2_lev) CYCLE
       uvh2_bx%nuv = uvh2_bx%nuv + 1
       IF(nvu < Vmin_bel)     Vmin_bel = nvu
       IF(njl+ndj < Jmin_bel) Jmin_bel = njl+ndj
       IF(nvu > Vmax_bel)     Vmax_bel = nvu
       IF(njl+ndj > Jmax_bel) Jmax_bel = njl+ndj
    ENDDO
    CLOSE (iread)

    Vmin_cel    = 0
    Jmin_cel    = 0
    Vmax_cel    = 0
    Jmax_cel    = 0
    uvh2_cx%nuv = 0
    fichier = TRIM(name_file_uvh2_werner)
    iread = GET_FILE_NUMBER()
    OPEN (iread, FILE = fichier, STATUS = 'old')
    READ (iread,"(i8)") nuv
    DO i = 1, nuv
       READ (iread,*) ii, ii, nvl, njl, nvu, ndj, uv1, uv2, uv3, uv4
       lev = index_VJ_H2(nvl,njl)
       IF (lev > NH2_lev) CYCLE
       uvh2_cx%nuv = uvh2_cx%nuv + 1
       IF(nvu < Vmin_cel)     Vmin_cel = nvu
       IF(njl+ndj < Jmin_cel) Jmin_cel = njl+ndj
       IF(nvu > Vmax_cel)     Vmax_cel = nvu
       IF(njl+ndj > Jmax_cel) Jmax_cel = njl+ndj
    ENDDO
    CLOSE (iread)

    Vmin_el = MIN(Vmin_bel, Vmin_cel)
    Jmin_el = MIN(Jmin_bel, Jmin_cel)
    Vmax_el = MAX(Vmax_bel, Vmax_cel)
    Jmax_el = MAX(Jmax_bel, Jmax_cel)

    !----------------------------------------------------
    ! Allocate and initialize tables
    !----------------------------------------------------
    ALLOCATE ( uvh2_bx%vl  (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%jl  (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%vu  (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%dj  (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%iw  (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%osc (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%wlg (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%ilft(uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%pbdi(uvh2_bx%nuv) )
    ALLOCATE ( uvh2_bx%fgk (uvh2_bx%nuv) )
    ALLOCATE ( uvh2_cx%vl  (uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%jl  (uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%vu  (uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%dj  (uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%iw  (uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%osc (uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%wlg (uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%ilft(uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%pbdi(uvh2_cx%nuv) )
    ALLOCATE ( uvh2_cx%fgk (uvh2_cx%nuv) )
    ALLOCATE ( aeinb(Vmin_bel:Vmax_bel,Jmin_bel:Jmax_bel,NH2_lev_max) )
    ALLOCATE ( aeinc(Vmin_cel:Vmax_cel,Jmin_cel:Jmax_cel,NH2_lev_max) )
    ALLOCATE ( aeitb(Vmin_bel:Vmax_bel,Jmin_bel:Jmax_bel)             )
    ALLOCATE ( aeitc(Vmin_cel:Vmax_cel,Jmin_cel:Jmax_cel,2)           )
    ALLOCATE ( tmunh2b(Vmin_el:Vmax_el,Jmin_el:Jmax_el)               )
    ALLOCATE ( tmunh2c(Vmin_el:Vmax_el,Jmin_el:Jmax_el)               )
    ALLOCATE ( bhmbx(Vmin_bel:Vmax_bel,Jmin_bel:Jmax_bel,NH2_lev)     )
    ALLOCATE ( bhmcx(Vmin_cel:Vmax_cel,Jmin_cel:Jmax_cel,NH2_lev)     )
    ALLOCATE ( cashm(NH2_lev_max,NH2_lev_max) )

    uvh2_bx%vl   = 0
    uvh2_bx%jl   = 0
    uvh2_bx%vu   = 0
    uvh2_bx%dj   = 0
    uvh2_bx%iw   = 0
    uvh2_bx%osc  = 0.0_dp
    uvh2_bx%wlg  = 0.0_dp
    uvh2_bx%ilft = 0.0_dp
    uvh2_bx%pbdi = 0.0_dp
    uvh2_bx%fgk  = 0.0_dp
    uvh2_cx%vl   = 0
    uvh2_cx%jl   = 0
    uvh2_cx%vu   = 0
    uvh2_cx%dj   = 0
    uvh2_cx%iw   = 0
    uvh2_cx%osc  = 0.0_dp
    uvh2_cx%wlg  = 0.0_dp
    uvh2_cx%ilft = 0.0_dp
    uvh2_cx%pbdi = 0.0_dp
    uvh2_cx%fgk  = 0.0_dp
    aeinb(Vmin_bel:Vmax_bel,Jmin_bel:Jmax_bel,1:NH2_lev_max) = 0.0_dp
    aeinc(Vmin_cel:Vmax_cel,Jmin_cel:Jmax_cel,1:NH2_lev_max) = 0.0_dp
    aeitb(Vmin_bel:Vmax_bel,Jmin_bel:Jmax_bel)               = 0.0_dp
    aeitc(Vmin_cel:Vmax_cel,Jmin_cel:Jmax_cel,1:2)           = 0.0_dp
    tmunh2b(Vmin_el:Vmax_el,Jmin_el:Jmax_el)                 = 0.0_dp
    tmunh2c(Vmin_el:Vmax_el,Jmin_el:Jmax_el)                 = 0.0_dp
    bhmbx(Vmin_bel:Vmax_bel,Jmin_bel:Jmax_bel,1:NH2_lev)     = 0.0_dp
    bhmcx(Vmin_cel:Vmax_cel,Jmin_cel:Jmax_cel,1:NH2_lev)     = 0.0_dp
    cashm(1:NH2_lev_max,1:NH2_lev_max)                       = 0.0_dp

    !----------------------------------------------------
    ! Read H2 Lyman transitions
    !----------------------------------------------------
    fichier = TRIM(name_file_uvh2_lyman)
    iread = GET_FILE_NUMBER()
    OPEN (iread, FILE = fichier, STATUS = 'old')
    READ (iread,"(i8)") nuv
    j = 0
    DO i = 1, nuv
       READ (iread,*) ii, ii, nvl, njl, nvu, ndj, uv1, uv2, uv3, uv4
       lev = index_VJ_H2(nvl,njl)
       nju = njl + ndj
       !-------------------------------------------------
       ! THE FOLLOWING TEST IT TAKEN FROM THE PDR CODE.
       ! ONLY ELECTRONIC TRANSITIONS TO LEVELS < NH2_LEV
       ! ARE TAKENT INTO ACCOUNT. THIS APPROACH PREVENT
       ! TO INCLUDE THE CASCADE THROUGH LEVELS BETWEEN
       ! NH2_LEV AND NH2_LEV_MAX.
       ! IT MAKES THE WHOLE COMPUTATION OF CASCADE COEFF
       ! IN CACOEF USELESS ... ASK JACQUES WHY SUCH A 
       ! TEST HAS BEEN INCLUDED IN THE PDR MODEL
       !-------------------------------------------------
       IF (lev > NH2_lev) CYCLE
       j = j + 1
       dumy =  coefosc * uv1 / uv2**2 * ((2.0_dp*njl+1.0_dp) / (2.0_dp*nju+1.0_dp))
       aeinb(nvu,nju,lev) = dumy
       aeitb(nvu,nju) = aeitb(nvu,nju) + dumy
       uvh2_bx%vl(j)    = nvl
       uvh2_bx%jl(j)    = njl
       uvh2_bx%vu(j)    = nvu
       uvh2_bx%dj(j)    = ndj
       uvh2_bx%osc(j)   = uv1
       uvh2_bx%wlg(j)   = uv2
       uvh2_bx%ilft(j)  = uv3
       uvh2_bx%pbdi(j)  = uv4
       tmunh2b(nvu,nju) = uv3
    ENDDO
    CLOSE (iread)

    !----------------------------------------------------
    ! Read H2 Werner transitions
    !----------------------------------------------------
    fichier = TRIM(name_file_uvh2_werner)
    iread = GET_FILE_NUMBER()
    OPEN (iread, FILE = fichier, STATUS = 'old')
    READ (iread,"(i8)") nuv
    j = 0
    DO i = 1, nuv
       READ (iread,*) ii, ii, nvl, njl, nvu, ndj, uv1, uv2, uv3, uv4
       lev = index_VJ_H2(nvl,njl)
       nju = njl + ndj
       !-------------------------------------------------
       ! SAME REMARK AS ABOVE
       !-------------------------------------------------
       IF (lev > NH2_lev) CYCLE
       j = j + 1
       dumy =  coefosc * uv1 / uv2**2 * ((2.0_dp*njl+1.0_dp) / (2.0_dp*nju+1.0_dp))
       aeinc(nvu,nju,lev) = dumy
       IF (ndj == 0) THEN
          aeitc(nvu,nju,2) = aeitc(nvu,nju,2) + dumy
       ELSE
          aeitc(nvu,nju,1) = aeitc(nvu,nju,1) + dumy
       ENDIF
       uvh2_cx%vl(j)    = nvl
       uvh2_cx%jl(j)    = njl
       uvh2_cx%vu(j)    = nvu
       uvh2_cx%dj(j)    = ndj
       uvh2_cx%osc(j)   = uv1
       uvh2_cx%wlg(j)   = uv2
       uvh2_cx%ilft(j)  = uv3
       uvh2_cx%pbdi(j)  = uv4
       tmunh2c(nvu,nju) = uv3
    ENDDO
    CLOSE (iread)

    !----------------------------------------------------
    ! Compute cascade coefficients (CASHM, BHMBX, BHMCX)
    !----------------------------------------------------
    CALL CACOEF

    !----------------------------------------------------
    ! Allocate and initialize excitation rates via 
    ! electronic pumping and dissociation rates
    !----------------------------------------------------
    ALLOCATE ( Pij_H2(1:NH2_lev,1:NH2_lev)     )
    ALLOCATE ( Pij_H2_old(1:NH2_lev,1:NH2_lev) )
    ALLOCATE ( sum_Pij_H2(1:NH2_lev)           )
    ALLOCATE ( pdh2vJ(1:NH2_lev)               )
    Pij_H2(1:NH2_lev,1:NH2_lev)     = 0.0_dp
    Pij_H2_old(1:NH2_lev,1:NH2_lev) = 0.0_dp
    sum_Pij_H2(1:NH2_lev)           = 0.0_dp
    pdh2vJ(1:NH2_lev)               = 0.0_dp
  END SUBROUTINE READ_H2_ELECTRONIC


  SUBROUTINE CACOEF
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_H2_ELECTRONIC TRANSITIONS
    ! purpose :
    !    compute H2 cascade coefficients
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    bhmbx, bhmcx : cascade coefficients for fluorescence Lyman and Werner
    !    cashm        : cascade coefficients inside the X state
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_TECHCONFIG
    USE MODULE_PHYS_VAR, ONLY   : NH2_lev

    IMPLICIT NONE
    REAL (KIND=dp), DIMENSION (NH2_lev_max,NH2_lev_max) :: cscad  ! probability of spontaneous de-excitation from 1 level to another
    INTEGER                                             :: levu   ! Index of the upper level
    INTEGER                                             :: levl   ! Index of the lower level
    INTEGER                                             :: nv0    ! V value of H2_lev(NH2_lev)
    INTEGER                                             :: nj0    ! J value of H2_lev(NH2_lev)
    INTEGER                                             :: ndj0
    INTEGER                                             :: ndj    ! dummy integer - Delta-J
    INTEGER                                             :: nvu    ! dummy integer - upper V level
    INTEGER                                             :: nju    ! dummy integer - upper J level
    INTEGER                                             :: njl    ! dummy integer - lower J level
    INTEGER                                             :: nj     ! dummy integer - J level
    REAL (KIND=dp)                                      :: total  ! dummy real to compute total de-excitation probability
    REAL (KIND=dp)                                      :: total1 ! dummy real to compute total de-excitation probability (lambda doublet)
    REAL (KIND=dp)                                      :: total2 ! dummy real to compute total de-excitation probability (lambda doublet)
    REAL (KIND=dp)                                      :: dviei  ! dummy real to compute total de-excitation probability
    REAL (KIND=dp)                                      :: dviei1 ! dummy real to compute total de-excitation probability (lambda doublet)
    REAL (KIND=dp)                                      :: dviei2 ! dummy real to compute total de-excitation probability (lambda doublet)
    INTEGER                                             :: k      ! dummy integer

    !----------------------------------------------------
    ! Get V and J quantum numbers of the uppest ro-vib
    ! levels of H2 computed properly in the code
    !----------------------------------------------------
    nv0 = H2_lev(NH2_lev)%V
    nj0 = H2_lev(NH2_lev)%J

    !----------------------------------------------------
    ! Compute cascade coefficients
    ! Viala, Roueff & Abgrall (1988), A&A, 190, 215-236
    !
    ! - NH2_lev_max is the maximum number of ro-vib
    !   levels of H2 in the ground electronic state
    ! - NH2_lev is the number of ro-vib levels of H2 for
    !   which we compute level populations
    ! Level population between NH2_lev and NH2_lev_max
    ! are not directly computed. Because these levels
    ! are mainly populated only by fluorescence and 
    ! chemistry, we can use the cascade formalism to
    ! take into account their contribution to the 
    ! populations of levels <= NH2_lev
    ! CAREFUL : cascade formalism breaks if collisions
    !           become an important excitation pathways
    !           of lev between NH2_lev and  NH2_lev_max
    !----------------------------------------------------

    !----------------------------------------------------
    ! Cascade coefficients for levels upper than NH2_lev
    ! Computed by recurrence starting from upper levels
    ! (Eq. 33)
    ! cashm(levu,levl) = probability that level levu 
    ! de-excite towards levl by any radiative path
    !----------------------------------------------------
    cscad(1:NH2_lev_max,1:NH2_lev_max) = 0.0_dp

    DO levl = NH2_lev_max, NH2_lev+1, -1
       ! By definition (Eq. 31)
       cscad(levl,levl) = 1.0_dp 
       DO levu = levl+1, NH2_lev_max
          total = 0.0_dp
          DO k = levl+1, levu
             IF (sum_Aij_H2(k) == 0.0_dp) CYCLE
             total = total + cscad(levu,k) * Aij_H2(k,levl) / sum_Aij_H2(k)
          ENDDO
          cscad(levu,levl) = total
       ENDDO
    ENDDO
  
    DO levl = 1, NH2_lev
       DO levu = NH2_lev+1, NH2_lev_max
          total = 0.0_dp
          DO k = NH2_lev+1, levu
             IF (sum_Aij_H2(k) == 0.0_dp) CYCLE
             total = total + cscad(levu,k) * Aij_H2(k,levl) / sum_Aij_H2(k)
          ENDDO
          cashm(levu,levl) = total
       ENDDO
    ENDDO

    !----------------------------------------------------
    ! Verification of cashm
    !   Sum over lower levels must be 1
    !   (total number of molecules is constant)
    !----------------------------------------------------
    DO levu = NH2_lev+1, NH2_lev_max
       nvu = H2_lev(levu)%V
       nju = H2_lev(levu)%J
       total = 1.0_dp
       DO levl = 1, NH2_lev
          total = total - cashm(levu,levl)
       ENDDO
       IF (ABS(total) > 1.0e-8_dp .AND.  total /= 1.0_dp) THEN
          PRINT *, 'pb CASHM!:', nvu, nju, total, levu
          STOP
       ENDIF
    ENDDO

    !----------------------------------------------------
    ! Compute deexcitation from excited electronic states
    !   Vmax_bel : highest vib. quantum number of B state
    !   Vmax_cel : highest vib. quantum number of C state
    !   Jmax_bel : highest rot. quantum number of B state
    !   Jmax_cel : highest rot. quantum number of C state
    !
    ! Compute bhmbx(nvu, nju, levl) and bhmcx(nvu, nju, levl)
    !   => proba that level vu, Ju of B (or C) state 
    !      de-excite towards levl by any path
    !----------------------------------------------------

    !----------------------------------------------------
    ! Lyman transitions - compute bhmbx
    !----------------------------------------------------
    ! v of the level of the B state
    DO nvu = Vmin_bel, Vmax_bel
       ! J of the level of the B state
       DO nju = Jmin_bel, Jmax_bel
          ! Inverse lifetime of level nvu, nju of B state (sum of Aij)
          dviei = aeitb(nvu,nju)
          IF (dviei == 0.0_dp) CYCLE
          DO levl = 1, NH2_lev
             total = 0.0_dp
             DO levu = NH2_lev+1, NH2_lev_max
                total = total + aeinb(nvu,nju,levu) * cashm(levu,levl)
             ENDDO
             total = total + aeinb(nvu,nju,levl)
             bhmbx(nvu,nju,levl) = total / dviei
          ENDDO
       ENDDO
    ENDDO

    !----------------------------------------------------
    ! Verification of bhmbx
    !----------------------------------------------------
    DO nvu = Vmin_bel, Vmax_bel
       DO nju = Jmin_bel, Jmax_bel
          dviei = aeitb(nvu,nju)
          IF (dviei == 0.0_dp) CYCLE
          total = 1.0_dp
          DO levl = 1, NH2_lev
             njl = H2_lev(levl)%J
             total = total - bhmbx(nvu,nju,levl)
          ENDDO
          IF (ABS(total) > 1.0e-8_dp) THEN
             PRINT *, 'pb BHMBX!:',nvu,nju,total
             STOP
          ENDIF
       ENDDO
    ENDDO

    !----------------------------------------------------
    ! Werner transitions - compute bhmcx
    !   Beware of Lambda-doubling!
    !   Each sub-level is identified only using Delta-J
    !   parity. This is a bit tricky (and dangerous).
    !   Should we be more systematic?
    !----------------------------------------------------
    DO nvu = Vmin_cel, Vmax_cel
       DO nju = Jmin_cel, Jmax_cel
          dviei1 = aeitc(nvu,nju,1)
          dviei2 = aeitc(nvu,nju,2)
          IF (dviei1 == 0.0_dp .AND. dviei2 == 0.0_dp) CYCLE
          DO levl = 1, NH2_lev
             njl = H2_lev(levl)%J
             ndj0 = ABS(nju - njl)
             IF (MOD(ndj0,2) == 0 .AND. dviei2 /= 0.0_dp) THEN
                total = 0.0_dp
                DO levu = NH2_lev+1, NH2_lev_max
                   nj = H2_lev(levu)%J
                   ndj = ABS(nju - nj)
                   IF (ndj /= 0) CYCLE
                   total = total + aeinc(nvu,nju,levu) * cashm(levu,levl)
                ENDDO
                total = total + aeinc(nvu,nju,levl)
                bhmcx(nvu,nju,levl) = total / dviei2
             ELSE IF (MOD(ndj0,2) == 1 .AND. dviei1 /= 0.0_dp) THEN
                total = 0.0_dp
                DO levu = NH2_lev+1, NH2_lev_max
                   nj = H2_lev(levu)%J
                   ndj = ABS(nju - nj)
                   IF (ndj /= 1) CYCLE
                   total = total + aeinc(nvu,nju,levu) * cashm(levu,levl)
                ENDDO
                total = total + aeinc(nvu,nju,levl)
                bhmcx(nvu,nju,levl) = total / dviei1
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  
    !----------------------------------------------------
    ! Verification of bhmcx
    !----------------------------------------------------
    DO nvu = Vmin_cel, Vmax_cel
       DO nju = Jmin_cel, Jmax_cel
          dviei1 = aeitc(nvu,nju,1)
          dviei2 = aeitc(nvu,nju,2)
          IF (dviei1 == 0.0_dp .AND. dviei2 == 0.0_dp) CYCLE
          IF (dviei1 /= 0.0_dp ) THEN
             total1 = 1.0_dp
             DO levl = 1, NH2_lev
                njl = H2_lev(levl)%J
                ndj = ABS(nju - njl)
                IF (MOD(ndj,2) == 1) THEN
                   total1 = total1 - bhmcx(nvu,nju,levl)
                ENDIF
             ENDDO
          ENDIF
          IF (dviei2 /= 0.0_dp ) THEN
             total2 = 1.0_dp
             DO levl = 1, NH2_lev
                njl = H2_lev(levl)%J
                ndj = ABS(nju - njl)
                IF (MOD(ndj,2) == 0) THEN
                   total2 = total2 - bhmcx(nvu,nju,levl)
                ENDIF
             ENDDO
          ENDIF
          IF ( (ABS(total1) > 1.0e-8_dp .OR. ABS(total2) > 1.0e-8_dp) ) THEN
             PRINT *, 'CX!:',nvu, nju, total1, total2
             STOP
          ENDIF
       ENDDO
    ENDDO

    ! ----------------------------------------------------------------------------------
    ! CHECK THE BHMBX AND BHMCX COEFFICIENTS (Only for debug purposes - BG 08/2016)
    ! ----------------------------------------------------------------------------------
    ! fichier = TRIM(out_dir)//'check_bhmbx_shock.res'
    ! OPEN(iwrtmp,file=TRIM(fichier),status='unknown')
    ! DO i = 1, NH2_lev
    !    DO j = Jmin_bel, Jmax_bel
    !       DO k = Vmin_bel, Vmax_bel
    !          WRITE(iwrtmp,'(3(I2,2X), ES10.3)') i, j, k, bhmbx(k,j,i)
    !       ENDDO
    !    ENDDO
    !    WRITE(iwrtmp,*)
    ! ENDDO
    ! CLOSE(iwrtmp)
    ! fichier = TRIM(out_dir)//'check_bhmcx_shock.res'
    ! OPEN(iwrtmp,file=TRIM(fichier),status='unknown')
    ! DO i = 1, NH2_lev
    !    DO j = Jmin_cel, Jmax_cel
    !       DO k = Vmin_cel, Vmax_cel
    !          WRITE(iwrtmp,'(3(I2,2X), ES10.3)') i, j, k, bhmcx(k,j,i)
    !       ENDDO
    !    ENDDO
    !    WRITE(iwrtmp,*)
    ! ENDDO
    ! CLOSE(iwrtmp)
    ! WRITE(*,*) NH2_lev, (Jmax_bel-Jmin_bel+1), (Vmax_bel-Vmin_bel+1)
    ! WRITE(*,*) NH2_lev, (Jmax_cel-Jmin_cel+1), (Vmax_cel-Vmin_cel+1)
    ! ----------------------------------------------------------------------------------

  END SUBROUTINE CACOEF


  SUBROUTINE COMPUTE_OP_H2
    !---------------------------------------------------------------------------
    ! called by :
    !     INIT_ROVIB_H2
    !     DIFFUN
    ! purpose :
    !     Computes ortho:para H2 ratio from H2_lev%density. Densities of
    !     ortho-H2 and para-H2 are also calculated.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     op_H2, Dens_orthoH2, Dens_paraH2
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : op_H2, NH2_lev
    IMPLICIT NONE

    Dens_paraH2  = SUM(H2_lev(1:NH2_lev)%density,mask=(MOD(H2_lev(1:NH2_lev)%J,2)==0))
    Dens_orthoH2 = SUM(H2_lev(1:NH2_lev)%density,mask=(MOD(H2_lev(1:NH2_lev)%J,2)>0))
    op_H2 = Dens_orthoH2/Dens_paraH2

  END SUBROUTINE COMPUTE_OP_H2


  FUNCTION NAME_H2_LINE(Vup,Vlow,Jlow,line_type) RESULT(name)
    !---------------------------------------------------------------------------
    ! called by :
    !    READ_H2_LINES
    ! purpose :
    !    computes the name of one line of the form 1-0S(1)
    ! subroutine/function needed :
    ! input variables :
    !    * Vup, Jup  -> (V,J) of the upper level  (two digits : < 100)
    !    * Jlow      -> J of th lower level       (two digits : < 100)
    !    * line_type -> 'S', 'O', 'Q'
    ! output variables :
    ! results :
    !    name -> name of the line, without blanks
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG), INTENT(in) :: Vup, Vlow, Jlow
    CHARACTER(LEN=1), INTENT(in)   :: line_type ! 'S', 'O' or 'Q'
    CHARACTER(LEN=2) :: str_Vup, str_Vlow,str_Jlow
    CHARACTER(LEN=9) :: name

    WRITE(str_Vup,'(I2)')Vup
    WRITE(str_Vlow,'(I2)')Vlow
    WRITE(str_Jlow,'(I2)')Jlow

    ! return the name without blanks
    name = TRIM(ADJUSTL(str_Vup)) // '-' // TRIM(ADJUSTL(str_Vlow))  // line_type // &
         '(' // TRIM(ADJUSTL(str_Jlow)) // ')'

  END FUNCTION NAME_H2_LINE



  SUBROUTINE EVOLUTION_H2
    !---------------------------------------------------------------------------
    ! called by :
    !    COMPUTE_H2
    ! purpose :
    !    Calculate the source term (cm-3.s-1) for rovibrational levels of H2,
    !    using Aij_H2 (radiation, s-1) read in READ_H2_LINE and Cij_H2
    !    (collisions, s-1) calculated in this subroutine.
    !    The result, YN_rovib_H2, is used in COMPUTE_H2 and in DIFFUN.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rovib_H2
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR, ONLY : Tn, ABS_DeltaV, NH2_lev, NH2_lev_var, H_H2_flag
    USE MODULE_GRAINS, ONLY : Mgrain, rsq_grm, dens_grc
    USE MODULE_CHEMICAL_SPECIES, ONLY : Dens_H, Dens_H2, Dens_Hplus, Dens_He, mass_H2, mass_hplus
    USE MODULE_CONSTANTS, ONLY: pi, Zero, kB
    USE NUM_REC, ONLY: splint
    USE MODULE_DEBUG_JLB

    IMPLICIT NONE
    REAL(KIND=DP), DIMENSION(NH2_lev,NH2_lev) :: Cij_H2           ! probability of excitation by collision (s-1)
    REAL(KIND=DP), DIMENSION(NH2_lev)         :: sum_Cij_H2       ! sum of Cij_H2 (s-1) starting from a given level
    REAL(KIND=DP), DIMENSION(NH2_lev,NH2_lev) :: Aij_plus_Cij_H2  ! sum of the radiative and collisional terms
    INTEGER(KIND=LONG)                        :: Nup, Jup, Vup    ! upper level
    INTEGER(KIND=LONG)                        :: Nlow, Jlow, Vlow ! lower level
    INTEGER(KIND=LONG)                        :: ABS_deltaJ       ! | delta J | between two levels
    INTEGER(KIND=LONG)                        :: Level1, Level2   ! dummy integers
    INTEGER(KIND=LONG)                        :: i, j             ! dummy integers
    REAL(KIND=DP)                             :: deexcit_to_excit ! deexcit_to_excit is the Boltzmann factor 
    REAL(KIND=DP)                             :: excit_rate       ! excitation rate from one level to another
    REAL(KIND=DP)                             :: deexcit_rate     ! deexcitation rate from one level to another
    REAL(KIND=DP)                             :: decal, Tdd, Tdd2, gam, gam1, gam2, delta_E
    REAL(KIND=DP)                             :: coeff_H_NR, coeff_H_R, coeff_H, coeff_He, coeff_H2, coeff_Hplus
    REAL(KIND=DP)                             :: dv_kms, Vth, prob_excit, excit_grain, Vin_H2
    REAL(KIND=DP), DIMENSION(1:25)            :: vin

    REAL(KIND=DP), DIMENSION(NH2_lev,NH2_lev) :: ec_H, ec_H2, ec_He, ec_Hp, ec_gr
    REAL(KIND=DP)                             :: bth_H, bth_H2, bth_He, bth_Hp, bth_gr
    REAL(KIND=DP)                             :: alp_Hp, alp_gr

    DO i = 1, 25
       vin(i) = 2.0_dp * DBLE(i)
    ENDDO

    !--- initialization ---
    Cij_H2          = Zero ! collisions (2 dimensions)
    sum_Cij_H2      = Zero ! collisions (1 dimension)
    Aij_plus_Cij_H2 = Zero ! collisions + radiation (2 dimensions)
    YN_rovib_H2     = Zero ! source terms for H2-level populations (cm-3.s-1)

    bth_H  = 0.0_DP
    bth_H2 = 0.0_DP
    bth_He = 0.0_DP
    bth_Hp = 0.0_DP
    bth_gr = 0.0_DP

    if (NH2_lev_var<2) then
       return
    endif
    !---------------------------------------------------------
    ! A. Calculate the probability of excitation by collisions
    !    Cij_H2 (s-1)
    !---------------------------------------------------------
    DO Nup = 2, NH2_lev ! loop on upper levels
       Jup = H2_lev(Nup)%J
       Vup = H2_lev(Nup)%V
       DO Nlow = 1, Nup-1 ! loop on lower levels
          Jlow = H2_lev(Nlow)%J
          Vlow = H2_lev(Nlow)%V
          ABS_deltaJ = ABS(Jup - Jlow)
          delta_E = H2_lev(Nup)%energy - H2_lev(Nlow)%energy

          !--- deexcit_to_excit is the Boltzmann factor allowing to ---
          !--- get the excitation rate from the de-excitation rate. ---
!         deexcit_to_excit = TEXP(-MIN(delta_E/Tn, 180._DP)) * &
!              H2_lev(Nup)%weight/ H2_lev(Nlow)%weight
          deexcit_to_excit = TEXP(-delta_E/Tn) * &
               H2_lev(Nup)%weight/ H2_lev(Nlow)%weight

          !--------------------------------------------------------------------------
          ! 1) fit of H2 + H collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          if (H_H2_flag == "DRF") then
            decal = 1.0_DP
          else
            IF (MOD(ABS_deltaJ,2) == 0) THEN
              decal = 1.0_DP
            ELSE
              decal = 0.3_DP
            ENDIF
          endif

!  JLB - 5 V 2003 - pb avec MM a TRES haute temperature
!         Tdd = decal + Tn * 1.0e-3_dp
          Tdd = decal + min(Tn * 1.0e-3_dp, 5.0e1_dp)
          Tdd2 = Tdd*Tdd

          ! 1.1) non-reactive collisions (always possible)
          !      Warning: reactive collision are mixed in Martin & Mandy,
          !               and are thus included in (badly-named) coeff_H_NR

          gam = rate_H_H2(1,Nup,Nlow)        &
              + rate_H_H2(2,Nup,Nlow) / Tdd  &
              + rate_H_H2(3,Nup,Nlow) / Tdd2 &
              + rate_H_H2(4,Nup,Nlow) * Tdd
          coeff_H_NR = tpow(gam)

          coeff_H_R = 0.0_DP
          ! 1.2) reactive collisions are added (David Flower's rule) if H_H2_flag = DRF
          !      and for Delta J = 2 only if H_H2_flag =  BOTH (complement quantum rates)
          IF (H_H2_flag == "DRF" .AND. Vup /= Vlow) THEN
             IF (MOD(ABS_deltaJ,2) == 0) THEN
                coeff_H_R = tpow(gam) * &
                     TEXP(-MAX(0.0_DP,(3900.0_DP-delta_E)/Tn))
             ELSE
                IF (Jlow > 0) THEN
                   Level1 = index_VJ_H2(Vlow,Jlow-1)
                ELSE
                   Level1 = index_VJ_H2(Vlow,Jlow+1)
                ENDIF
                Level2 = index_VJ_H2(Vlow,Jlow+1)
                IF (Level2 > NH2_lev) THEN
                   Level2 = index_VJ_H2(Vlow,Jlow-1)
                ENDIF
                gam1 = rate_H_H2(1,Nup,Level1)         &
                     + rate_H_H2(2,Nup,Level1) / Tdd   &
                     + rate_H_H2(3,Nup,Level1) / Tdd2  &
                     + rate_H_H2(4,Nup,Level1) * Tdd
                gam2 = rate_H_H2(1,Nup,Level2)         &
                     + rate_H_H2(2,Nup,Level2) / Tdd   &
                     + rate_H_H2(3,Nup,Level2) / Tdd2  &
                     + rate_H_H2(4,Nup,Level2) * Tdd
                coeff_H_R = (tpow(gam1) + tpow(gam2)) * 0.5_DP * &
                     TEXP(-MAX(0.0_DP,(3900.0_DP-delta_E)/Tn))
                IF (MOD(Jup,2) == 1) THEN
                   coeff_H_R = coeff_H_R / 3.0_DP
                ENDIF

             ENDIF
          ENDIF
!         IF (H_H2_flag == "BOTH" .AND. Vup /= Vlow  .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow) == .TRUE.) THEN
          IF (H_H2_flag == "BOTH" .AND. Vup /= Vlow  .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow)) THEN
             coeff_H_R = tpow(gam) * &
                     TEXP(-MAX(0.0_DP,(3900.0_DP-delta_E)/Tn))
          ENDIF

          ! 1.3) ortho-para transfer with H from Schofield (only if MM are not included)
          !    + David Flower prescription (only for Delta v=0)
          IF (H_H2_flag == "DRF" &
!       .OR. (H_H2_flag == "BOTH" .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow) == .TRUE.)) THEN
        .OR. (H_H2_flag == "BOTH" .AND. MOD(ABS_deltaJ,2) == 0 .AND. mask_H_H2(Nup,Nlow))) THEN
            IF (Vup == Vlow .AND. ABS_deltaJ < 3) THEN
               IF (MOD(Jup,2) == 1 .AND. MOD(ABS_deltaJ,2) == 1) THEN
                  coeff_H_R = 8.0D-11 * TEXP(-3900.0_DP/Tn) / 3.0_DP
               ELSE
                  coeff_H_R = 8.0D-11 * TEXP(-3900.0_DP/Tn)
               ENDIF
            ENDIF
          ENDIF

          ! coeff_H_R is 0.0 if H_H2_flag = MM

          coeff_H = coeff_H_NR + coeff_H_R

          !--------------------------------------------------------------------------
          ! 2) fit of H2 + He collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          gam = rate_He_H2(1,Nup,Nlow)        &
              + rate_He_H2(2,Nup,Nlow) / Tdd  &
              + rate_He_H2(3,Nup,Nlow) / Tdd2 &
              + rate_He_H2(4,Nup,Nlow) * Tdd
          coeff_He = tpow(gam)

          !--------------------------------------------------------------------------
          ! 3) fit of H2 + para-H2 collision rates.
          !    Tdd is a reduced variable that can be "shifted" with a different value
          !    for reactive collisions (delta(J) odd) and for non-reactive ones.
          !    This allows for different asymptotic forms at low temperature.
          !--------------------------------------------------------------------------
          gam = rate_H2_H2(1,Nup,Nlow)        &
              + rate_H2_H2(2,Nup,Nlow) / Tdd  &
              + rate_H2_H2(3,Nup,Nlow) / Tdd2 &
              + rate_H2_H2(4,Nup,Nlow) * Tdd
          coeff_H2 = tpow(gam)

          !------------------------------------------------
          ! 4) ortho-para transfer with ions (Gerlich data)
          !    only in v=0
          !------------------------------------------------
          IF (Vup == 0 .AND. ABS_deltaJ == 1) THEN
             coeff_Hplus = Gerlich(Jup)
          ELSE
             coeff_Hplus = 0.0_DP
          ENDIF

          !---------------------------
          ! 5) Calculate Cij_H2 (s-1)
          !---------------------------
          deexcit_rate = coeff_H2 * Dens_H2        &
                       + coeff_He * Dens_He        &
                       + coeff_Hplus * Dens_Hplus  &
                       + coeff_H * Dens_H
          excit_rate   = deexcit_rate * deexcit_to_excit

          !----------------------------------
          ! 6) Calculate excitation by grains
          !----------------------------------
          ! Excitation probabilites are fitted by a spline, see f77 , for
          ! DeltaV > VIN (first impact velocity at which excitation
          ! probability becomes non zero), VTH = VIN-2km/s.
          ! Exponential form at DeltaV < VIN, scaled by PRG0 to pass through point at VIN
          ! Raw excitation rates are in r_raw_GR_H2
          ! Second derivatives, as computed by SPLINE are in d2r_raw_GR_H2
          
          ! Note : use of ABS_DeltaV => grains are supposed to be charged, or at least
          ! that their dynamics follows that of the ionized fluid.

          dv_kms = ABS_DeltaV  * 1.0D-5
          Vin_H2 = vin_GR_H2(Vup,Jup)
          Vth = Vin_H2 - 2.0_DP

          excit_grain = 0.0_DP
          prob_excit = 0.0_DP
          !if ((Vlow==0) .and. ((Jlow==0) .or. (Jlow==1)) .and. (mod(Jup-Jlow,2)==0)) then
          if ((Vlow==0) .and. (Jlow <= 7) .and. (mod(Jup-Jlow,2)==0)) then
            if (dv_kms >= Vin_H2) then
               prob_excit = splint(vin,r_raw_GR_H2(:,Vup,Jup),d2r_raw_GR_H2(:,Vup,Jup),dv_kms)
               prob_excit = max (prob_excit,0.0d0)

            else if (dv_kms > 1.0e-20) then
               ! 8.02*dv_kms**2 = 0.1*Ts where 3/2kTs = 1/2mH2 dv**2
               ! extrapolation to velocity difference lower than the first point in the excitation grid
               ! at Vin = vin_GR_H2(Vup,Jup)

               prob_excit = TEXP(delta_E/8.02*(1.0_DP/vin_GR_H2(Vup,Jup)**2 - 1.0_DP/dv_kms**2))
               prob_excit = prob_excit*pgr0_GR_H2(Vup,Jup)

            else

               prob_excit = 0.0_DP

            end if
            excit_grain = prob_excit * dens_grc * pi * rsq_grm * ABS_DeltaV

          else
            excit_grain = 0.0_DP
          end if

          Cij_H2(Nlow,Nup) = Cij_H2(Nlow,Nup) + excit_rate + excit_grain
          Cij_H2(Nup,Nlow) = Cij_H2(Nup,Nlow) + deexcit_rate

          !--------------------------------------------------
          ! 7) Calculate separated matrices for each collider
          !---------------------------------------------------
          ec_H (Nup,Nlow) = coeff_H * Dens_H
          ec_He(Nup,Nlow) = coeff_He * Dens_He
          ec_Hp(Nup,Nlow) = coeff_Hplus * Dens_Hplus
          ec_H2(Nup,Nlow) = coeff_H2 * Dens_H2
          ec_gr(Nup,Nlow) = 0.0_dp
          ec_H (Nlow,Nup) = ec_H (Nup,Nlow) * deexcit_to_excit
          ec_He(Nlow,Nup) = ec_He(Nup,Nlow) * deexcit_to_excit
          ec_Hp(Nlow,Nup) = ec_Hp(Nup,Nlow) * deexcit_to_excit
          ec_H2(Nlow,Nup) = ec_H2(Nup,Nlow) * deexcit_to_excit
          ec_gr(Nlow,Nup) = excit_grain

       END DO ! end of loop on Nlow
    END DO ! end of loop on Nup

    !--- calculate sum_Cij_H2 ---
    DO Nup=1,NH2_lev
       sum_Cij_H2(Nup) = SUM(Cij_H2(Nup,1:NH2_lev))
    END DO

    !---------------------------------------------------------------
    ! B. Add (de)excitation probabilities from collisions (Cij_H2)
    !    and from radiation (Aij_H2)
    !---------------------------------------------------------------
    Aij_plus_Cij_H2(1:NH2_lev,1:NH2_lev) = Aij_H2(1:NH2_lev,1:NH2_lev) &
                                         + Cij_H2(1:NH2_lev,1:NH2_lev) &
                                         + Pij_H2(1:NH2_lev,1:NH2_lev)
    DO i = 1, NH2_lev ! diagonal terms
       Aij_plus_Cij_H2(i,i) = Aij_plus_Cij_H2(i,i) &
                            - sum_Aij_H2(i) &
                            - sum_Cij_H2(i) &
                            - sum_Pij_H2(i)
    END DO

    !----------------------------------------------------------
    ! C. Calculate evolution terms for populations in cm-3.s-1
    !    = matrix multiplication of population density and
    !    (de)excitation probabilities.
    !----------------------------------------------------------
    YN_rovib_H2(1:NH2_lev) = &
         MATMUL(H2_lev(1:NH2_lev)%density, &
         Aij_plus_Cij_H2(1:NH2_lev,1:NH2_lev))

    !----------------------------------------------------------
    ! D. Compute the cooling terms in erg.cm-3.s-1
    !    BG 2016 - separate their contributions for the neutral
    !    positively ionized, and negatively ionized fluids
    !----------------------------------------------------------
    DO i=1,NH2_lev-1
       DO j=i+1,NH2_lev
          delta_E = H2_lev(i)%energy - H2_lev(j)%energy
          bth_H  = bth_H  + (ec_H(i,j)  * H2_lev(i)%density - ec_H(j,i)  * H2_lev(j)%density) * delta_E
          bth_H2 = bth_H2 + (ec_H2(i,j) * H2_lev(i)%density - ec_H2(j,i) * H2_lev(j)%density) * delta_E
          bth_He = bth_He + (ec_He(i,j) * H2_lev(i)%density - ec_He(j,i) * H2_lev(j)%density) * delta_E
          bth_Hp = bth_Hp + (ec_Hp(i,j) * H2_lev(i)%density - ec_Hp(j,i) * H2_lev(j)%density) * delta_E
          bth_gr = bth_gr + (ec_gr(i,j) * H2_lev(i)%density - ec_gr(j,i) * H2_lev(j)%density) * delta_E
       ENDDO
    ENDDO

    alp_Hp = mass_Hplus / (mass_Hplus + mass_H2)
    ! alp_gr = (Mgrain/dens_grc) / ((Mgrain/dens_grc) + mass_H2)
    alp_gr = Mgrain / (Mgrain + dens_grc*mass_H2) ! Tabone 09/18 other way to write the same equation to handle low dens_grc = 0._DP

    cooling_n_H2 = kB * ( bth_H + bth_H2 + bth_He + bth_Hp * alp_Hp + bth_gr * alp_gr )
    cooling_i_H2 = kB * ( bth_Hp * (1.0_DP-alp_Hp) + bth_gr * (1.0_DP-alp_gr) )
    cooling_e_H2 = 0.0_DP

  END SUBROUTINE EVOLUTION_H2



  SUBROUTINE COMPUTE_H2
    !---------------------------------------------------------------------------
    ! called by :
    !    SOURCE
    ! purpose :
    !    Calculate the source terms for H2 levels (cm-3.s-1), H2 internal energy
    !    (erg.cm-3.s-1), the H2 lines (erg.cm-3.s-1) and the corresponding
    !    cooling rate (erg.cm-3.s-1).
    ! subroutine/function needed :
    !    EVOLUTION_H2
    ! input variables :
    ! output variables :
    ! results :
    !    YN_rovib_H2, H2_Energy, H2_lines%emiss, cooling_H2
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : kB
    USE MODULE_PHYS_VAR, ONLY : NH2_lev,NH2_lev_var,Tn,nH,cool_KN
    USE MODULE_DEBUG_JLB
    USE MODULE_CHEMICAL_SPECIES, ONLY : Dens_H,Dens_H2
    use refh2
    IMPLICIT NONE

    !-----------------------------------------
    ! In EVOLUTION_H2, the source terms 
    ! YN_rovib_H2 (cm-3.s-1) is calculated, 
    ! with treatment of radiation & collision
    !-----------------------------------------
    CALL EVOLUTION_H2

    if (cool_KN==2) then
       cooling_H2=0d0
       return
    endif
    !-----------------------------------------
    ! If not enough levels, compute H2 cooling
    ! as an interpolation on a 4D grid
    ! see refH2.f90
    !-----------------------------------------
    if (NH2_lev_var<2) then
       call coolH2(log10(Tn),log10(nH),log10(Dens_H/Dens_H2),3.0_dp,cooling_H2)
       cooling_H2=Dens_H2*10.**cooling_H2
       return
    endif
    !-----------------------------------------
    ! Calculate the emissivity (erg.cm-3.s-1)
    ! of each H2 quadrupolar line.
    ! result in H2_lines%emiss
    !-----------------------------------------
    ! first in K.cm-3.s-1
    H2_lines(1:NH2_lines)%emiss = &
         H2_lines(1:NH2_lines)%DeltaE * H2_lines(1:NH2_lines)%Aij * &
         H2_lev(H2_lines(1:NH2_lines)%Nup)%density
    ! change into erg.cm-3.s-1
    H2_lines%emiss=H2_lines%emiss*kB

    !-------------------------------------------------------------
    ! cooling rate cooling_H2 (erg.cm-3.s-1) = sum of emissivities
    !-------------------------------------------------------------
    ! New integrations of the evolution equations - dont take into
    ! account internal energy of H2 in the energy conservation eq
    ! --- old ---
    ! cooling_H2 = SUM(DBLE(H2_lines(1:NH2_lines)%emiss))
    ! --- new ---
    cooling_H2 = - ( cooling_n_H2 + cooling_i_H2 + cooling_e_H2 )

    !---------------------------------------------------------------
    ! source term for internal energy of H2 H2_energy (erg.cm-3.s-1)
    !---------------------------------------------------------------
    H2_energy = kB * SUM(YN_rovib_H2(1:NH2_lev) * H2_lev(1:NH2_lev)%energy)

  END SUBROUTINE COMPUTE_H2

END MODULE MODULE_H2
