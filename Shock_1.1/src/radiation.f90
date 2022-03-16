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

MODULE MODULE_RADIATION
  !*****************************************************************************
  !** The module 'MODULE_RADIATION' contains variables and subroutines        **
  !** related to the UV (only for now) radiation field :                      **
  !**     * the parameters for the radiation field grid                       **
  !**     * INIT_WAVELENGTH (wavelength initialization subroutine)            **
  !**     * INIT_RADIATION  (radiation field initialization subroutine)       **
  !**     * EXTEND_WLG (extend the wavelength table to a new value)           **
  !**     * LOCATE_LEFT (locate a value in a table with closest left index)   **
  !**     * auxiliary subroutines (integral, ...)                             **
  !** The routines related to the computation of FGK transfer have been taken **
  !** from the Meudon PDR code and adapted for time dependent simulations     **
  !*****************************************************************************
  USE MODULE_TECHCONFIG

  IMPLICIT NONE
  INCLUDE "precision.f90"

  !---------------------------------------------------------------------------------------
  ! input wavelength / position grid files
  !---------------------------------------------------------------------------------------
  CHARACTER(LEN=lenfilename)                   :: name_file_grid

  ! --------------------------------------------------------------------------------------
  ! parameters and variables link to the radiation field grid
  ! --------------------------------------------------------------------------------------
  REAL (KIND=dp)                               :: wlth0              ! start of Wavelength grid
  REAL (KIND=dp)                               :: wlthm              ! end of Wavelength grid
  REAL (KIND=dp)                               :: hiop               ! Hydrogen ionization potential
  REAL (KIND=dp)                               :: wl_lym             ! Lyman limit
  REAL (KIND=dp), PARAMETER                    :: wl_hab  =  2400.0_dp ! Limit of Habing Field
  REAL (KIND=dp), PARAMETER                    :: wl_Hmus = 14700.0_dp ! Ionization threshold for H- Tabone2017  14700.0_dp in Angstrom
  INTEGER,        PARAMETER                    :: nwlg_init = 100    ! initial number of WL points
  INTEGER                                      :: nwlg               ! number of WL points
  REAL (KIND=dp)                               :: pwlg               ! initial WL step (log scale)

  !--------------------------------------------------------------------------------------!
  ! radiation field variables : wavelength, energy density, intensity ...
  !--------------------------------------------------------------------------------------!
  INTEGER,        PARAMETER                    :: nangle = 20        ! number of directions for radiative transfer
  REAL (KIND=dp), DIMENSION(nangle)            :: angles             ! angles used for radiative transfer
  REAL (KIND=dp)                               :: dang               ! angle intervals

  INTEGER                                      :: nwlggrd            ! number of wavelength points in the radiation field grid file
  INTEGER                                      :: nposgrd            ! number of position   points in the radiation field grid file
  REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: wlggrd             ! wavelength points contained in the radiation field grid file
  REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: posgrd             ! position   points contained in the radiation field grid file
  REAL (KIND=dp), DIMENSION(:,:),  ALLOCATABLE :: urfgrd             ! radiation field energy density stored in the radiation field grid file

  REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: wlg                ! table of wavelength values
  REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: irf_free           ! ISRF radiation field specific intensity (erg cm-2 s-1 sr-1 A-1) in free space
  REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: irf_secpho         ! equivalent radiation field specific intensity (erg cm-2 s-1 sr-1 A-1) associated to secondary photons
  REAL (KIND=dp), DIMENSION(:,:),  ALLOCATABLE :: irf                ! radiation field specific intensity (erg cm-2 s-1 sr-1 A-1)
  REAL (KIND=dp), DIMENSION(:,:),  ALLOCATABLE :: irf_old            ! radiation field specific intensity (erg cm-2 s-1 sr-1 A-1) at previous step
  ! REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: irf_0              ! initial radiation field specific intensity (erg cm-2 s-1 sr-1 A-1) at previous step
  REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: urf                ! radiation field energy density (erg cm-3 A-1)
  REAL (KIND=dp), DIMENSION(:),    ALLOCATABLE :: nrf                ! radiation field flux (photon cm-2 s-1 A-1)
  
  REAL (KIND=dp)                               :: J_to_u             ! From mean Intensity (erg cm-2 s-1 Ang-1 sr-1) 
                                                                     ! to energy density (erg cm-3 Ang-1)
  REAL (KIND=dp)                               :: fluphfree          ! ISRF integrated photon flux (photon cm-2 s-1) integrated over 4pi
  REAL (KIND=dp)                               :: fluph              ! integrated photon flux (photon cm-2 s-1)
  REAL (KIND=dp)                               :: fluph0             ! initial value of fluph taking into account G0 and integrated over 2pi
  REAL (KIND=dp)                               :: nrjph              ! photon energy density

  REAL (KIND=dp), DIMENSION (:),   ALLOCATABLE :: qabso_D_rf         ! dust absorption coefficient on rad field grid
  REAL (KIND=dp), DIMENSION (:,:), ALLOCATABLE :: sigma_G_rf         ! gas  absorption coefficient on rad field grid
  REAL (KIND=dp), DIMENSION (:),   ALLOCATABLE :: doptdpth           ! increment of optical depth at all wavelength
  ! REAL (KIND=dp), DIMENSION (:),   ALLOCATABLE :: optdpth            ! total optical depth at all wavelength
  ! REAL (KIND=dp), DIMENSION (:),   ALLOCATABLE :: optdpth_old        ! total optical depth at all wavelength

  !--------------------------------------------------------------------------------------!
  ! Output 
  !--------------------------------------------------------------------------------------!
  CHARACTER (LEN=15)                           :: ISRF_name          !  name of FUV radiation field (used in outputs)

  !--------------------------------------------------------------------------------------!
  ! variables for FGK transfer
  !--------------------------------------------------------------------------------------!
  INTEGER                                      :: mfgk               ! Number of lines (H2 + ...) in FGK transfer
  REAL (KIND=dp), DIMENSION (:), ALLOCATABLE   :: fgkgrd             ! Line shielded radiation field

  ! --------------------------------------------------------------------------------------
  ! variables from "precision.f90" are private
  ! --------------------------------------------------------------------------------------
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE READ_RAD_GRID
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     read the radiation field grid provided by the user
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     nwlggrd, wlggrd
    !     nposgrd, posgrd
    !     urfgrd
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_VAR_VODE,  ONLY : length_max
    USE MODULE_TOOLS,     ONLY : GET_FILE_NUMBER

    IMPLICIT none

    CHARACTER(LEN=100) :: c_line      ! dummy chain of character to rid the grid file values
    LOGICAL            :: file_exists ! boolean - does file exists or not
    REAL(KIND=dp)      :: dum         ! dummy real variable
    REAL(KIND=dp)      :: dumz        ! dummy real variable

    INTEGER            :: i, j        ! dummy integers
    INTEGER            :: io

    ! ---------------------------------------
    ! If no grid file quit the routine
    ! ---------------------------------------
    IF( gridfile == 'none' ) THEN
       nwlggrd   = 0
       nposgrd   = 1
       ALLOCATE (posgrd(0:nposgrd))
       posgrd(0) = 0.0_dp
       posgrd(1) = length_max
       RETURN
    ENDIF

    ! ---------------------------------------
    ! Test the existence of the grid file
    ! ---------------------------------------
    iwrtmp = GET_FILE_NUMBER()
    name_file_grid = TRIM(data_dir) // TRIM(gridfile)
    INQUIRE(FILE=name_file_grid, EXIST=file_exists)
    IF( gridfile /= 'none' .AND. .NOT.file_exists ) THEN
       WRITE(*,*) "The grid file given in input doesnt exist"
       WRITE(*,*) "Please supply the file or write 'none' in"
       WRITE(*,*) "the shock input file to avoid using a grid"
       WRITE(*,*) "based radiation field"
       STOP
    ENDIF

    ! ---------------------------------------
    ! Count the number of wavelengths
    ! contained in the file and read them
    ! ---------------------------------------
    OPEN(iwrtmp, file=name_file_grid, status='OLD', access='SEQUENTIAL',form='FORMATTED',action='READ')
    READ(iwrtmp,*)
    nwlggrd = 0
    READ (iwrtmp, '(A)') c_line 
    DO WHILE (c_line /= ' ')
       nwlggrd = nwlggrd + 1
       READ (iwrtmp, '(A)') c_line 
    ENDDO

    ALLOCATE (wlggrd(nwlggrd))
    REWIND(iwrtmp)
    READ(iwrtmp,*)
    DO i = 1, nwlggrd
       READ(iwrtmp,*) dum, dum, wlggrd(i), dum
    ENDDO

    ! ---------------------------------------
    ! Count the number of positions
    ! contained in the file and read them
    ! along with the radiation field energy
    ! density
    ! ---------------------------------------
    REWIND(iwrtmp)
    nposgrd = 1
    DO
       READ (iwrtmp, '(A)', iostat=io) c_line
       IF (io /= 0) EXIT
       IF (c_line == ' ') THEN
          nposgrd = nposgrd + 1
       ENDIF
    ENDDO

    ALLOCATE (posgrd(0:nposgrd))
    ALLOCATE (urfgrd(nposgrd,nwlggrd))
    REWIND(iwrtmp)
    READ(iwrtmp,*)
    DO j = 1, nposgrd
       DO i = 1, nwlggrd
          READ(iwrtmp,*) dumz, posgrd(j), dum, urfgrd(j,i)
          IF( j == 1 ) THEN
             posgrd(j-1) = dumz
          ENDIF
          IF ( dumz /= posgrd(j-1) ) THEN
             WRITE(*,'(A,A)') "Error in the radiation grid file: ", TRIM(gridfile)
             WRITE(*,'(A)')   "   -> zstart (slab) /= zend (slab-1)"
             WRITE(*,'(A)')   "   -> please correct"
             STOP
          ENDIF
       ENDDO
       IF ( j /= nposgrd .AND. i /= nwlggrd) THEN
          READ(iwrtmp,*)
       ENDIF
    ENDDO

    ! ---------------------------------------
    ! Close input file
    ! ---------------------------------------
    CLOSE(iwrtmp)
    
  END SUBROUTINE READ_RAD_GRID


  SUBROUTINE INIT_WAVELENGTH
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     Initialise model dependant quantities related to radiation field
    !        - Wavelength grid
    !        - Incident ISRF (specific intensity)
    ! subroutine/function needed :
    !     EXTAND_WLG
    ! input variables :
    ! output variables :
    ! results :
    !     wlg
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    ! USE MODULE_H2,             ONLY : uvh2_bx, uvh2_cx, index_VJ_H2
    ! USE MODULE_PHYS_VAR,       ONLY : NH2_lev
    USE MODULE_CHEM_REACT,     ONLY : phdest, Ncross
    USE MODULE_DUST_TREATMENT, ONLY : jwlmin_gr, jwlmax_gr, np_w, wavgrd

    IMPLICIT NONE

    ! INTEGER            :: lev         ! H2 level number
    ! INTEGER            :: njl         ! J quantum number of H2
    ! INTEGER            :: nvl         ! V quantum number of H2
    
    INTEGER            :: jmin, jmax  ! minimal & max indices of photoreactions cross sections
                                      ! to include in the wavelength grid
    INTEGER            :: i, j


    !================================================================================
    ! Initialization of the radiation field grid
    !================================================================================

    ! ---------------------------------------
    ! Set initial grid step - BG 08/2016
    ! from Lyman (wlth0) to Habing (wlthm)
    ! limits !!Tabone can be changed!!
    ! ---------------------------------------
    J_to_u = 4.0_dp * pi / clum
    hiop   = 2.0_dp * pi*pi * qe**4 * me / (hp*hp)
    ! convert hiop in angstrom
    wl_lym = 1.0e8_dp * hp * clum / hiop
    wlth0  = wl_lym
    wlthm  = wl_Hmus  ! wl_Hmus wl_hab MODIF Tabone2017
    nwlg   = nwlg_init

    ! ---------------------------------------
    ! Initialisation of wavelengths
    ! ---------------------------------------
    pwlg = EXP(LOG(wlthm / wlth0) / DBLE(nwlg-1))
    ALLOCATE( wlg(1:nwlg) )
    wlg(1:nwlg) = wlth0 * ( pwlg**(/ (i, i=0,nwlg-1) /) )
    wlg(nwlg)   = wlthm


    !================================================================================
    ! Add grain absorption wavelength in the grid
    !================================================================================

    ! ---------------------------------------
    ! Find range in wavelength grid
    ! ---------------------------------------
    DO i = 2, np_w-1
       IF (wavgrd(i-1) < wlth0 .AND. wavgrd(i) >= wlth0) THEN
          jwlmin_gr = i
       ENDIF
       IF (wavgrd(i) <= wlthm .AND. wavgrd(i+1) > wlthm) THEN
          jwlmax_gr = i
          EXIT
       ENDIF
    ENDDO
    DO i = jwlmin_gr, jwlmax_gr
       CALL EXTAND_WLG (wavgrd(i))
    ENDDO


    !================================================================================
    ! Add photoreaction wavelength in the grid
    !================================================================================
    DO i = 1, Ncross
       ! ---------------------------------------
       ! Find the min index of wavelength to 
       ! store in the grid
       ! All wavelength must be >= wlth0 and
       ! <= wlthm
       ! ---------------------------------------
       jmin = 1
       DO WHILE( (jmin <= phdest(i)%npts) .AND. (phdest(i)%wl(jmin) < wlth0) )
          jmin = jmin + 1
       ENDDO
       jmax = phdest(i)%npts
       DO WHILE( (jmax >= 1) .AND. (phdest(i)%wl(jmax) > wlthm) )
          jmax = jmax - 1
       ENDDO

       DO j = jmin, jmax
          CALL EXTAND_WLG (phdest(i)%wl(j))
       ENDDO
    ENDDO


    !================================================================================
    ! Add line wavelength in the grid
    !================================================================================

    ! ---------------------------------------
    ! H2 lines, only lowest levels
    ! ---------------------------------------
    !!!! DO i = 1, uvh2_bx%nuv
    !!!!    nvl = uvh2_bx%vl(i)
    !!!!    njl = uvh2_bx%jl(i)
    !!!!    lev = index_VJ_H2(nvl,njl)
    !!!!    IF (lev >= NH2_lev) CYCLE
    !!!!    CALL EXTAND_WLG (uvh2_bx%wlg(i))
    !!!! ENDDO

    !!!! DO i = 1, uvh2_cx%nuv
    !!!!    nvl = uvh2_cx%vl(i)
    !!!!    njl = uvh2_cx%jl(i)
    !!!!    lev = index_VJ_H2(nvl,njl)
    !!!!    IF (lev >= NH2_lev) CYCLE
    !!!!    CALL EXTAND_WLG (uvh2_cx%wlg(i))
    !!!! ENDDO


    !================================================================================
    ! Add wavelength from the input radiation grid file
    !================================================================================
    DO i = 1, nwlggrd
       CALL EXTAND_WLG (wlggrd(i))
    ENDDO

    !================================================================================
    ! From here, the set of wavelength is complete and should not be modified
    !================================================================================

    ! ---------------------------------------
    ! Allocate the optical depth table
    ! ---------------------------------------
    ALLOCATE(qabso_D_rf(1:nwlg))
    ALLOCATE(sigma_G_rf(1:Ncross,1:nwlg))
    ALLOCATE(doptdpth(1:nwlg))
    ! ALLOCATE(optdpth(1:nwlg))
    ! ALLOCATE(optdpth_old(1:nwlg))
    qabso_D_rf(1:nwlg)          = 0.0_dp
    sigma_G_rf(1:Ncross,1:nwlg) = 0.0_dp
    doptdpth(1:nwlg)            = 0.0_dp
    ! optdpth(1:nwlg)             = 0.0_dp
    ! optdpth_old(1:nwlg)         = 0.0_dp

  END SUBROUTINE INIT_WAVELENGTH


  SUBROUTINE INIT_RADIATION
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     Initialize the radiation field intensity, sum of different components
    !     compo 1 : FUV radiation field (Mathis or Draine)
    !     UV field is scaled by RAD
    ! subroutine/function needed :
    !     LOCATE_LEFT
    ! input variables :
    ! output variables :
    ! results :
    !     irf_free
    !     irf_secpho
    !     irf
    !     irf_old
    !     urf
    !     nrf
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_PHYS_VAR,  ONLY : F_ISRF, RAD, Zeta
    USE MODULE_TECHCONFIG

    IMPLICIT NONE

    REAL (KIND=dp), PARAMETER :: wlg_D = 1.0e5_dp  ! Draine field cutoff wavelength (extended)
    REAL (KIND=dp), PARAMETER :: wlg_V = 8000.0_dp ! Visible cutoff wavelength
    INTEGER                   :: idx_D             ! Draine field cutoff position in grid
    INTEGER                   :: idx_V             ! index in wavelength Grid of visible cutoff
    INTEGER                   :: idx_L             ! index in wavelength Grid of Lyman cutoff
    INTEGER                   :: idx_Hmus          ! index in wavelength Grid of Hmus cutoff

    INTEGER                   :: i, ii

    !==================================================================================
    ! Initialisation of radiation Field
    !==================================================================================
    ALLOCATE ( irf_free  (1:nwlg) )
    ALLOCATE ( irf_secpho(1:nwlg) )
    ALLOCATE ( irf       (1:nangle,1:nwlg) )
    ALLOCATE ( irf_old   (1:nangle,1:nwlg) )
    ! ALLOCATE ( irf_0     (1:nwlg) )
    ALLOCATE ( urf       (1:nwlg) )
    ALLOCATE ( nrf       (1:nwlg) )
    irf_free  (1:nwlg) = 0.0_dp
    irf_secpho(1:nwlg) = 0.0_dp
    irf       (1:nangle,1:nwlg) = 0.0_dp
    irf_old   (1:nangle,1:nwlg) = 0.0_dp
    urf       (1:nwlg) = 0.0_dp
    nrf       (1:nwlg) = 0.0_dp

    !==================================================================================
    ! CREATION OF THE INCIDENT RADIATION FIELD
    !==================================================================================
    ! The standard interstellar radiation field contains here only two components
    ! - the far UV and the visible component
    !    Draine Formula: Specific intensity
    !    Radiation field in erg cm-2 s-1 A-1 sterad-1
    ! The code need the definition of the standard ISRF for the computation of 
    ! photoreactions and the definition of the secondary photon number
    !==================================================================================

    CALL LOCATE_LEFT (wlg, nwlg, wl_lym, idx_L)
    ! ---------------------------------------
    ! Mathis radiation field
    ! ---------------------------------------
    IF (F_ISRF == 1) THEN
       ISRF_name = "MMP83+JB  "
       !---------- Far UV (fit JLB Jan 2008) ----------
       CALL LOCATE_LEFT (wlg, nwlg, wlg_V, idx_V)
       irf_free(idx_L:idx_V) = ( TANH ( 4.07e-3_dp * wlg(idx_L:idx_V) - 4.5991_dp ) + 1.0_dp ) &
                             * 107.192_dp * wlg(idx_L:idx_V)**(-2.89_dp)
       irf_free(idx_V+1:nwlg) = 2.0_dp * 107.192_dp * wlg(idx_V+1:nwlg)**(-2.89_dp)
       !---------- Tabone 2017 inclue visible part of the Mathis radiation field ----------
       CALL LOCATE_LEFT (wlg, nwlg, wl_Hmus, idx_Hmus)
       DO ii = idx_L, idx_Hmus
          irf_free(ii)  = irf_free(ii) + 2.07e-14_dp * BLKB(wlg(ii),6184.0_dp) &
                                       + 1.09e-14_dp * BLKB(wlg(ii),6123.0_dp) &
                                       + 1.52e-12_dp * BLKB(wlg(ii),2539.0_dp)
       ENDDO
    ! ---------------------------------------
    ! Draine radiation field
    ! ---------------------------------------
    ELSE IF (F_ISRF == 2) THEN
       ISRF_name = "Draine  "
       ! Draine's formula extended to 10000.0 A to give a smooth transition to Visible
       ! 2 Fev 2012 - JLB - Still not enough for Rad = 1.0e5 (required by Leiden Workshop)
       !                    Extend to 1.0e5 A = 10 micron
       CALL LOCATE_LEFT (wlg, nwlg, wlg_D, idx_D)
       ! Draine 1978 radiation field
       ! Expression from M. Kopp PhD thesis page 284, from Sternberg and Dalgarno 1985.
       irf_free(idx_L:idx_D) = 6.36E7_dp    / wlg(idx_L:idx_D)**4 &
                             - 1.0237E11_dp / wlg(idx_L:idx_D)**5 &
                             + 4.0812E13_dp / wlg(idx_L:idx_D)**6
       ! irf_free is now in erg cm-2 s-1 str-1 A-1
       irf_free (1:nwlg) = irf_free (1:nwlg) / (4.0 * pi) 
    ENDIF

    ! ---------------------------------------
    ! Check the non-scaled spectrum
    ! Only for debug purposes (BG 08/2016)
    ! ---------------------------------------
    fichier = TRIM(out_dir) // 'check_isrf.res'
    OPEN (iwrtmp, FILE=TRIM(fichier), STATUS="unknown")
    WRITE(iwrtmp,'(A)') "# wavelength (A)    spec intensity"
    DO ii = 1, nwlg
       WRITE (iwrtmp,'(2(ES16.9,2X))') wlg(ii), irf_free(ii)
    ENDDO
    CLOSE (iwrtmp)

    ! ---------------------------------------
    ! Equivalent radiation field used for
    ! secondary photons - G0 = 1e-5
    ! Careful : the spectrum is set to that
    ! of the standard ISRF even though it is
    ! stupid. It is just a template for now
    ! ---------------------------------------
    irf_secpho(1:nwlg) = 1e-5_dp * Zeta / 1e-17 * irf_free(1:nwlg)

    ! ---------------------------------------
    ! Compute the energy density and photon
    ! flux over 4pi steradians.
    ! This values are the reference for G0=1
    ! and without considering any specific
    ! geometry
    ! ---------------------------------------
    ! Analytical integration
    urf(1:nwlg) = irf_free(1:nwlg) * 4.0_dp * pi / clum
    nrf(1:nwlg) = urf(1:nwlg) * wlg(1:nwlg) * 1e-8_dp / hp
    fluphfree = TAB_INTEGRATION (wlg, nrf, nwlg, wl_lym, wl_hab) ! Tabone to be modified
    WRITE(*,*) "fluphfree = ", fluphfree

    !=================================================================================
    ! Build full spectrum
    !=================================================================================
    ! define the angles - from 0 to PI/2
    dang = pi / 2.0 / DBLE(nangle-1)
    DO i = 1, nangle
      angles(i) = DBLE(i-1) * dang
    ENDDO

    ! ---------------------------------------
    ! scale radiation field with RAD
    ! Compute the energy density and photon
    ! flux - integration over 2pi steradians
    ! The 1D parallel shock is supposed to 
    ! be illuminated from one side only
    ! Integration over the angles should
    ! therefore be performed for theta
    ! varying between 0 and PI/2 (and not PI)
    ! and phi varying between 0 and 2PI
    ! ---------------------------------------
    DO i = 1, nangle
       irf(i,1:nwlg)      = RAD * irf_free(1:nwlg)
       irf_old(i,1:nwlg)  = irf(i,1:nwlg)
    ENDDO
    ! irf_0(1:nwlg) = irf(1,1:nwlg)
    ! ---------------------------------------
    ! Analytic integration over angles
    ! ---------------------------------------
    ! urf(1:nwlg) = irf_0(1:nwlg) * 2 * pi / clum
    ! ---------------------------------------
    ! Numerical integration over angles
    ! ---------------------------------------
    DO ii = 1, nwlg
       urf(ii) = 0
       urf(ii) = urf(ii) + irf(1,ii) / 2.0_dp * sin(angles(1))
       DO i = 1, nangle
          urf(ii) = urf(ii) + irf(i,ii) * sin(angles(i))
       ENDDO
       urf(ii) = urf(ii) + irf(nangle,ii) / 2.0_dp * sin(angles(nangle))
       urf(ii) = urf(ii) * dang * 2*pi / clum
    ENDDO
    ! ---------------------------------------
    ! Compute the photon flux over the entire 
    ! range (911 - 2400 A)
    ! ---------------------------------------
    DO ii = 1, nwlg
       nrf(ii) = urf(ii) * wlg(ii) * 1e-8_dp / hp
    ENDDO
    ! ---------------------------------------
    ! Compute associated integrated fluxes
    ! ---------------------------------------
    fluph = TAB_INTEGRATION (wlg, nrf, nwlg, wl_lym, wl_hab)
    fluph0 = fluph
    nrjph = TAB_INTEGRATION (wlg, urf, nwlg, wl_lym, wl_hab) * clum

  END SUBROUTINE INIT_RADIATION


  SUBROUTINE REINIT_RADIATION

    USE MODULE_CONSTANTS
    USE MODULE_PHYS_VAR,  ONLY : RAD

    IMPLICIT none

    INTEGER :: i
    INTEGER :: ii

    DO i = 1, nangle
       irf(i,1:nwlg)     = RAD * irf_free(1:nwlg)
       irf_old(i,1:nwlg) = irf(i,1:nwlg)
    ENDDO
    ! Numerical integration
    DO ii = 1, nwlg
       urf(ii) = 0
       urf(ii) = urf(ii) + irf(1,ii) / 2.0_dp * sin(angles(1))
       DO i = 1, nangle
          urf(ii) = urf(ii) + irf(i,ii) * sin(angles(i))
       ENDDO
       urf(ii) = urf(ii) + irf(nangle,ii) / 2.0_dp * sin(angles(nangle))
       urf(ii) = urf(ii) * dang * 2*pi / clum
    ENDDO
    ! Compute the photon flux over the entire range (911 - 2400 A)
    DO ii = 1, nwlg
       nrf(ii) = urf(ii) * wlg(ii) * 1e-8_dp / hp
    ENDDO
    fluph = TAB_INTEGRATION (wlg, nrf, nwlg, wl_lym, wl_hab)
    fluph0 = fluph
    nrjph = TAB_INTEGRATION (wlg, urf, nwlg, wl_lym, wl_hab) * clum

  END SUBROUTINE REINIT_RADIATION


  SUBROUTINE INTERP_QABSO(n, qabso, tab_wlg)
     !---------------------------------------------------------------------------
     ! called by :
     !     INITIALIZE
     ! purpose :
     !     Interpolate the dust absorption coefficients over all
     !     the wavelengths included in an input wavelength table
     ! subroutine/function needed :
     ! input variables :
     !     n       -> number of wavelength points
     !     tab_wlg -> wavelength grid
     ! output variables :
     ! results :
     !     qabso -> grain absorption coefficients
     !---------------------------------------------------------------------------
     USE MODULE_DUST_TREATMENT, ONLY : jwlmin_gr, wavgrd, Q_abs

     IMPLICIT none
     
     INTEGER,                       INTENT(in)  :: n
     REAL(KIND=dp), DIMENSION(1:n), INTENT(out) :: qabso
     REAL(KIND=dp), DIMENSION(1:n), INTENT(in)  :: tab_wlg
     INTEGER                                    :: i, j

     j = jwlmin_gr-1
     DO i = 1, n
        qabso(i) = 0.0_dp
        DO
           IF( (wavgrd(j) <= tab_wlg(i)).and.(wavgrd(j+1) > tab_wlg(i)) ) THEN
              EXIT
           ELSE
              j = j + 1
           ENDIF
        ENDDO
        qabso(i) = LOG(Q_abs(j)) + LOG( Q_abs(j+1) / Q_abs(j) ) &
                 * LOG( tab_wlg(i) / wavgrd(j) ) / LOG( wavgrd(j+1) / wavgrd(j) )
        qabso(i) = EXP(qabso(i))
     ENDDO
  END SUBROUTINE INTERP_QABSO


  SUBROUTINE INTERP_SIGMA_G_RF
     !---------------------------------------------------------------------------
     ! called by :
     !     INITIALIZE
     ! purpose :
     !     Interpolate the gas absorption cross section over the wavelength
     !     grid for all known photoionization or photodissociation cross sec.
     ! subroutine/function needed :
     ! input variables :
     ! output variables :
     ! results :
     !     sigma_G_rf
     !---------------------------------------------------------------------------
     USE MODULE_CHEM_REACT, ONLY : CROSS_SEC, phdest, Ncross

     IMPLICIT none

     INTEGER        :: i, j, k
     INTEGER        :: npts
     INTEGER        :: i1, i2
     REAL (KIND=dp) :: wl1, wl2, sig1, sig2

     sigma_G_rf(:,:) = 0.0_dp
     DO k = 1, Ncross
        npts = phdest(k)%npts
        DO j = 1, npts-1
           wl1  = phdest(k)%wl(j)
           wl2  = phdest(k)%wl(j+1)
           sig1 = phdest(k)%sigma(j)
           sig2 = phdest(k)%sigma(j+1)
           CALL LOCATE_LEFT (wlg, nwlg, wl1, i1)
           IF(i1 < 1) i1 = 1
           CALL LOCATE_LEFT (wlg, nwlg, wl2, i2)
           IF(i2 < 1) i2 = 1
           IF(i1 == i2) CYCLE
           DO i = i1, i2
              sigma_G_rf(k,i) = 1.0_dp / ( wl2 - wl1 ) &
                              * ( sig1 * ABS(wlg(i) - wl2) &
                                + sig2 * ABS(wlg(i) - wl1) )
           ENDDO
        ENDDO
     ENDDO
  END SUBROUTINE INTERP_SIGMA_G_RF


  SUBROUTINE COMPUTE_PHOTOELEC
    !---------------------------------------------------------------------------
    ! called by :
    !    DIFFUN
    ! purpose :
    !    compute the photoelectric effect rates
    ! subroutine/function needed :
    !    VECT_INT
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_GRAINS
    USE MODULE_CHEM_REACT, ONLY : phele_reac, Nphele

    IMPLICIT none

    REAL (KIND=dp), PARAMETER         :: x_la = 1.0e-6_dp ! Photon attenuation length (Bakes & Tielens)
    REAL (KIND=dp), PARAMETER         :: x_le = 1.0e-7_dp ! Photon attenuation length (Bakes & Tielens)
    REAL (KIND=dp)                    :: x_nc, x_nc_t
    REAL (KIND=dp)                    :: aalp
    REAL (KIND=dp)                    :: zetph
    REAL (KIND=dp)                    :: fy
    REAL (KIND=dp)                    :: bt_fy
    REAL (KIND=dp)                    :: cte_J
    REAL (KIND=dp)                    :: x_ip, wl_cut_IP, cte_IP
    INTEGER                           :: iz
    INTEGER                           :: i

    REAL (KIND=dp), DIMENSION(1:nwlg) :: vector
    REAL (KIND=dp)                    :: intg1, intg2


    x_nc   = (r_grm-d_ads)**3.0_dp * rho_grc + (r_grm**3.0_dp-(r_grm-d_ads)**3.0_dp) * rho_grm
    x_nc   = x_nc * 4.0_dp * pi / (3.0_dp * amu * 12.0_dp)
    x_nc_t = x_nc**(1.0_dp/3.0_dp)

    aalp = r_grm / x_la + r_grm / x_le
    zetph = r_grm / x_la
    fy = aalp*aalp - 2.0_dp * aalp + 2.0_dp - 2.0_dp * EXP(-aalp)
    fy = fy / (zetph*zetph - 2.0_dp * zetph + 2.0_dp - 2.0_dp * EXP(-zetph))
    bt_fy = fy * (zetph / aalp)**2

    cte_J = pi * rsq_grm * 0.14_dp * bt_fy

    DO i = 1, Nphele
       iz = phele_reac(i)%charge
       IF (iz < 0) THEN
          x_ip = 0.6_dp + (iz + 0.5_dp) * 11.1_dp / x_nc_t
       ELSE
          x_ip = 4.4_dp + (iz + 0.5_dp) * 11.1_dp / x_nc_t
       ENDIF
       IF (x_ip > 0.0_dp) THEN
          wl_cut_IP = 1.0e8_dp*hp*clum/(EVerg*x_ip)
       ELSE
          wl_cut_IP = wlg(nwlg)
       ENDIF
       cte_IP = 1.0_dp / wl_cut_IP

       vector(1:nwlg) = qabso_D_rf(1:nwlg)
       intg1          = VECT_INT (nwlg, vector, wlg, nrf, wl_cut_IP)
       vector(1:nwlg) = qabso_D_rf(1:nwlg) * wlg(1:nwlg)
       intg2          = VECT_INT (nwlg, vector, wlg, nrf, wl_cut_IP)

       phele_reac(i)%rate = cte_J * (intg1 - cte_IP * intg2)
    ENDDO

  END SUBROUTINE COMPUTE_PHOTOELEC


  SUBROUTINE COMPUTE_GRAIN_TEMP
    !---------------------------------------------------------------------------
    ! called by :
    !    DIFFUN
    ! purpose :
    !    compute the photoelectric effect rates
    ! subroutine/function needed :
    !    FUNC_TGRAIN
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_GRAINS, ONLY : Tgrain

    IMPLICIT none

    REAL(KIND=DP)             :: x, y, dx, f, df
    ! REAL(KIND=DP), PARAMETER  :: eps = 1.0e-06
    INTEGER                   :: IMAX = 10
    INTEGER                   :: i

    x = Tgrain

    DO i = 1, IMAX
       CALL FUNC_TGRAIN(x, f, df)
       dx = f / df
       y = x - dx

       !--------------------------------------------------
       ! Pierre's tip : when a Newton Raphson is called
       ! inside a Newton Raphson, we must provide smooth
       ! quantities. In particular, computing the numerical
       ! jacobian in dvode requires to always perform NR
       ! steps when computing Tgrain.
       ! For those reasons, it is better to always apply
       ! the same number of steps in the NR to find Tgrain
       ! and not use IF statements to exit NR
       ! Hence comment all tests below
       !--------------------------------------------------

       !!!!! ! Convergence
       !!!!! IF ( abs(dx / x) < eps )  THEN
       !!!!!    EXIT
       !!!!! ENDIF

       !!!!! ! Test - prevent exploration of negative values
       !!!!! IF (y <= 0.0_DP) THEN
       !!!!!    x = x / 10.0_DP
       !!!!! ELSE
       !!!!!    x = y
       !!!!! ENDIF

       x = y
    ENDDO

    !-----------------------------------------------------
    ! Put a lower limit to the grain temperature
    ! To avoid setting an infrared radiation field
    !-----------------------------------------------------
    ! Tgrain = x
    Tgrain = MAX(x,15.0_dp)

  END SUBROUTINE COMPUTE_GRAIN_TEMP


  SUBROUTINE FUNC_TGRAIN(x, f, df)
     !---------------------------------------------------------------------------
     ! called by :
     !    COMPUTE_GRAIN_TEMP
     ! purpose :
     !    f is a function whose root gives the grain temperature
     !    FUNC_TGRAIN computes the value of f and of its derivative df at x
     ! subroutine/function needed :
     ! input variables :
     ! ouput variables :
     ! results :
     !    f(x) and df/dx(x)
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS, ONLY : clum, pi
     USE MODULE_PHYS_VAR,  ONLY : nH, Tn
     USE MODULE_GRAINS,    ONLY : dens_grc, rsq_grm

     IMPLICIT none

     REAL(KIND=DP),                   INTENT(IN)  :: x
     REAL(KIND=DP),                   INTENT(OUT) :: f, df
     INTEGER,       PARAMETER                     :: NptB = 100
     REAL(KIND=DP), PARAMETER                     :: wl1 = 3.0e4_dp !     3 micron (grain at 70 K -> lmax =  42 microns)
     REAL(KIND=DP), PARAMETER                     :: wl2 = 1.0e8_dp ! 10000 micron (grain at  5 K -> lmax = 600 microns
     REAL(KIND=DP)                                :: pwlbb
     REAL(KIND=DP), DIMENSION(1:NptB)             :: Qabso
     REAL(KIND=DP), DIMENSION(1:NptB)             :: wlg_bb
     REAL(KIND=DP), DIMENSION(1:NptB)             :: irf_bb
     REAL(KIND=DP), DIMENSION(1:NptB)             :: dirf_bb
     INTEGER                                      :: i
     REAL(KIND=DP)                                :: intg1, intg2, intg3, intg4
     REAL(KIND=DP)                                :: therm, dtherm

     pwlbb          = EXP(LOG(wl2 / wl1) / DBLE(NptB-1))
     wlg_bb(1:NptB) = wl1 * ( pwlbb**(/ (i, i=0,NptB-1) /) )
     wlg_bb(NptB)   = wl2

     DO i = 1, NptB
        irf_bb(i)  = BLKB(wlg_bb(i),x)
        dirf_bb(i) = DERBLKB(wlg_bb(i),x)
     ENDDO

     CALL INTERP_QABSO(NptB, Qabso, wlg_bb)

     !----------------------------------------------------
     ! Integrals of absorption over the radiation field
     ! - external radiation field
     ! - guessed secondary photon radiation field
     !----------------------------------------------------
     intg1 = clum   * VECT_INT (nwlg, qabso_D_rf, wlg, urf, wlthm)
     intg2 = 4 * pi * VECT_INT (nwlg, qabso_D_rf, wlg, irf_secpho, wlthm)
     !----------------------------------------------------
     ! Intergrals of emission at a given temperature
     ! - emission
     ! - derivative dE/dT
     !----------------------------------------------------
     intg3 = 4 * pi * VECT_INT (NptB, qabso, wlg_bb, irf_bb,  wl2)
     intg4 = 4 * pi * VECT_INT (NptB, qabso, wlg_bb, dirf_bb, wl2)
     !----------------------------------------------------
     ! Grain thermalization by collision with gas
     ! - energy gain
     ! - derivative dE/dT
     !----------------------------------------------------
     therm  =  3.5e-34_dp * SQRT(Tn) * (Tn - x) * nH**2.0_DP / dens_grc / (pi * rsq_grm)
     dtherm = -3.5e-34_dp * SQRT(Tn) * nH**2.0_DP / dens_grc / (pi * rsq_grm)

     !----------------------------------------------------
     ! add up the different terms in the F=0 equation (NR)
     !----------------------------------------------------
     f  = intg1 + intg2 + therm - intg3
     df = dtherm - intg4

  END SUBROUTINE FUNC_TGRAIN


  REAL(KIND=dp) FUNCTION  BLKB(x,T)
    !---------------------------------------------------------------------------
    ! called by :
    !     FUNC_TGRAIN
    ! purpose :
    !     compute the black body specific intensity in erg cm-2 s-1 A-1 sr-1
    !     BLKB(x,T) at a wavelength x (in angstrom) and for a temperature T
    ! subroutine/function needed :
    ! input variables :
    !     x -> wavelength (in A)
    !     T -> black body temperature (in K)
    ! output variables :
    ! results :
    !     BLKB = specific intensity
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : hcsurk, hp, clum, l_hu

    IMPLICIT NONE

    REAL (KIND = dp), INTENT (IN)      :: x          ! Wavelength (Angstrom)
    REAL (KIND = dp), INTENT (IN)      :: T          ! Temperature (K)

    IF ( hcsurk * 1.0e8_dp / (x*T) >= l_hu ) THEN
       BLKB = 0.0_dp
    ELSE
       BLKB = 2.0_dp * hp * clum**2 * 1.0e32_dp / x**5 / ( DEXP(hcsurk*1.0e8_dp / (x*T)) - 1.0_dp )
    ENDIF

  END FUNCTION BLKB


  REAL(KIND=dp) FUNCTION  DERBLKB(x,T)
    !---------------------------------------------------------------------------
    ! called by :
    !     FUNC_TGRAIN
    ! purpose :
    !     compute the derivative of black body specific intensity
    !     BLKB(x,T) (fot Newton-Raphson determination of dust temperatures)
    ! subroutine/function needed :
    ! input variables :
    !     x -> wavelength (in A)
    !     T -> black body temperature (in K)
    ! output variables :
    ! results :
    !     DERBLKB = BLKB derivative
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : hcsurk, hp, clum, l_hu

     IMPLICIT NONE

     REAL (KIND = dp), INTENT (IN)      :: x          ! Wavelength (Angstrom)
     REAL (KIND = dp), INTENT (IN)      :: T          ! Temperature (K)
     REAL (KIND = dp)                   :: exp

     IF ( hcsurk*1.0e8_dp / (x*T) >= l_hu ) THEN
        DERBLKB = 0.0_dp
     ELSE
        exp = DEXP(hcsurk*1.0e8_dp / (x*T))
        DERBLKB = 2.0_dp * hp * clum**2 * 1.0e40_dp * hcsurk / (T**2 * x**6) * exp / (exp - 1.0_dp)**2
     ENDIF

  END FUNCTION DERBLKB


  SUBROUTINE SIMPLE_TRANSFER
     !---------------------------------------------------------------------------
     ! called by :
     !     INITIALIZE
     !     MHD_VODE
     !     DIFFUN
     ! purpose :
     !     Compute the UV continuum radiative transfer assuming only absorption
     !     by gas and dust, and using a single size and absorption coefficient
     !     for dust
     ! subroutine/function needed :
     !     TAB_INTEGRATION
     ! input variables :
     ! output variables :
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS
     ! USE NUM_REC,         ONLY : gamminc0

     IMPLICIT none

     ! REAL(KIND=DP) :: x
     ! REAL(KIND=DP) :: g
     INTEGER       :: i, ii

     ! ---------------------------------------
     ! compute the specific intensity at all
     ! angle
     ! ---------------------------------------
     DO i = 1, nangle
       irf(i,1:nwlg) = irf(i,1:nwlg) * EXP( - doptdpth(1:nwlg) / cos(angles(i)) )
     ENDDO

     ! ---------------------------------------
     ! compute energy density and photon flux 
     ! with analytical integration over angles
     ! ---------------------------------------
     ! DO ii = 1, nwlg
     !    x = optdpth(ii)
     !    CALL gamminc0(x, g, .FALSE.)
     !    urf(ii) = irf_0(ii) * 2 * pi / clum * ( EXP(-x) - x * g )
     !    nrf(ii) = urf(ii) * wlg(ii) * 1e-8_dp / hp
     ! ENDDO
     ! ---------------------------------------
     ! compute energy density and photon flux 
     ! with numerical integration over angles
     ! ---------------------------------------
     DO ii = 1, nwlg
        urf(ii) = 0
        urf(ii) = urf(ii) + irf(1,ii) / 2.0_dp * sin(angles(1))
        DO i = 2, nangle-1
           urf(ii) = urf(ii) + irf(i,ii) * sin(angles(i))
        ENDDO
        urf(ii) = urf(ii) + irf(nangle,ii) / 2.0_dp * sin(angles(nangle))
        urf(ii) = urf(ii) * dang * 2*pi / clum
        nrf(ii) = urf(ii) * wlg(ii) * 1e-8_dp / hp
     ENDDO

     ! ---------------------------------------
     ! Compute integrated quantities
     ! ---------------------------------------
     fluph = TAB_INTEGRATION (wlg, nrf, nwlg, wl_lym, wl_hab)
     nrjph = TAB_INTEGRATION (wlg, urf, nwlg, wl_lym, wl_hab) * clum

  END SUBROUTINE SIMPLE_TRANSFER


  SUBROUTINE SET_RAD_GRID(jpos)
     !---------------------------------------------------------------------------
     ! called by :
     !     INITIALIZE
     !     MHD_VODE
     !     DIFFUN
     ! purpose :
     !     set urf and nrf with the radiation field provided in the grid file 
     !     Compute the associated radiation fluxes
     ! subroutine/function needed :
     !     LOCATE_LEFT
     !     TAB_INTEGRATION
     ! input variables :
     !     index 
     ! output variables :
     ! results :
     !---------------------------------------------------------------------------
     USE MODULE_CONSTANTS

     IMPLICIT none

     INTEGER, INTENT(in) :: jpos ! position of the radiation field slab
     REAL(KIND=dp)       :: dum
     INTEGER             :: i
     INTEGER             :: ii

     urf(1:nwlg) = 0.0_dp

     ! ---------------------------------------
     ! Add the radiation field provided in 
     ! the grid file 
     ! - first position only
     ! ---------------------------------------
     DO ii = 1, nwlg
        IF( wlg(ii) <  wlggrd(1) .OR. wlg(ii) >  wlggrd(nwlggrd) ) CYCLE
        CALL LOCATE_LEFT (wlggrd, nwlggrd, wlg(ii), i)
        IF( i == nwlggrd ) THEN
           urf(ii) = urfgrd(jpos,i)
        ELSE
           dum     = ( wlg(ii) - wlggrd(i) ) / ( wlggrd(i+1) - wlggrd(i) )
           urf(ii) = urfgrd(jpos,i) + ( urfgrd(jpos,i+1) - urfgrd(jpos,i) ) * dum
        ENDIF
     ENDDO

     ! ---------------------------------------
     ! Compute the photon flux over the entire
     ! range (911 - 2400 A)
     ! ---------------------------------------
     DO ii = 1, nwlg
        nrf(ii) = urf(ii) * wlg(ii) * 1e-8_dp / hp
     ENDDO

     fluph = TAB_INTEGRATION (wlg, nrf, nwlg, wl_lym, wl_hab)
     nrjph = TAB_INTEGRATION (wlg, urf, nwlg, wl_lym, wl_hab) * clum

     ! WRITE(*,*) fluph, nrjph

  END SUBROUTINE SET_RAD_GRID


  SUBROUTINE EXTAND_WLG (x)
    !---------------------------------------------------------------------------
    ! called by :
    !     INIT_WAVELENGTH
    ! purpose :
    !     extand the list of wavelength with an additional value
    ! subroutine/function needed :
    !     LOCATE_TAB
    ! input variables :
    !     x : value to add in the table
    ! output variables :
    ! results :
    !     new ordered wlg table
    !     nwlg -> nwlg + 1
    !---------------------------------------------------------------------------

    IMPLICIT NONE

    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: tabtmp
    REAL(KIND=dp),               INTENT(IN)   :: x
    INTEGER                                   :: i

    ! ---------------------------------------
    ! Locate x
    ! if x is outside wlg, or x is
    ! already in wlg no need to extand
    ! ---------------------------------------
    CALL LOCATE_LEFT (wlg, nwlg, x, i)
    IF ( (i == -1) .OR. (i == nwlg) ) THEN
       RETURN
    END IF
    IF ( wlg(i) == x ) THEN
       RETURN
    ENDIF

    ! ---------------------------------------
    ! Copy wlg in temporary table
    ! add the new point in wlg
    ! ---------------------------------------
    ALLOCATE(tabtmp(1:nwlg))
    tabtmp(1:nwlg) = wlg(1:nwlg)

    DEALLOCATE(wlg)
    ALLOCATE(wlg(1:nwlg+1))
    wlg(1:i)        = tabtmp(1:i)
    wlg(i+1)        = x
    wlg(i+2:nwlg+1) = tabtmp(i+1:nwlg)
    nwlg               = nwlg + 1

    DEALLOCATE(tabtmp)

  END SUBROUTINE EXTAND_WLG


  SUBROUTINE LOCATE_LEFT (tab, n, x, j)
    !---------------------------------------------------------------------------
    ! called by :
    !     EXTAND_WLG
    ! purpose :
    !     locate index j of tab such that tab(j) <= x < tab(j+1)
    !     NEED : tab should be sorted in ascending order
    !     Note : if x is outside tab, return 0 or n
    ! subroutine/function needed :
    ! input variables :
    !     tab : a real table
    !     x : value to locate in the table
    ! output variables :
    ! results :
    !     j : index such that tab(j) <= x < tab(j+1)
    !---------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,                        INTENT (IN)  :: n
    REAL (KIND=dp), DIMENSION(1:n), INTENT (IN)  :: tab
    REAL (KIND=dp),                 INTENT (IN)  :: x
    INTEGER,                        INTENT (OUT) :: j
    INTEGER                                      :: jl
    INTEGER                                      :: ju
    INTEGER                                      :: jm

    IF (x >= tab(n)) THEN
      j = n
    ELSE IF (x < tab(1)) THEN
      j = 0
    ELSE
      jl = 1
      ju = n
      DO WHILE ((ju-jl) > 1)
        jm = (ju + jl) / 2
        IF (x < tab(jm)) THEN
          ju = jm
        ELSE
          jl = jm
        ENDIF
      ENDDO
      j = jl
    ENDIF

  END SUBROUTINE LOCATE_LEFT


  REAL(KIND=dp) FUNCTION  TAB_INTEGRATION (tabx, taby, n, xmin, xmax)
     !--------------------------------------------------------------------------
     ! Integrate a function using the trapezium rules.
     ! With uneven step
     ! *** Tabone 12/2017 integration range included to spot
     !     only 91-250nm range to compute FUV field.
     !--------------------------------------------------------------------------

     IMPLICIT NONE

     INTEGER,                        INTENT (IN) :: n      ! number of points
     REAL (KIND=dp), DIMENSION(1:n), INTENT (IN) :: tabx
     REAL (KIND=dp), DIMENSION(1:n), INTENT (IN) :: taby
     REAL (KIND=dp),                 INTENT (IN) :: xmin
     REAL (KIND=dp),                 INTENT (IN) :: xmax
     INTEGER                                     :: l
     INTEGER                                     :: u
     INTEGER                                     :: i

     CALL LOCATE_LEFT (tabx,n,xmin,l)
     CALL LOCATE_LEFT (tabx,n,xmax,u)

     TAB_INTEGRATION = 0.5_dp * SUM( (/ ( (tabx(i+1) - tabx(i) ) &
                     * (taby(i+1) + taby(i) ), i = l, u-1) /) )

  END FUNCTION TAB_INTEGRATION


  REAL(KIND=dp) FUNCTION VECT_INT (n, vect, tabx, taby, wlim)
    !---------------------------------------------------------------------------
    ! called by :
    !     COMPUTE_PHOTOELEC
    !     COMPUTE_FUNC_GRAIN
    ! purpose :
    !     integrate a vector of real over an rf structure
    !     up to a wavelength wlim
    ! subroutine/function needed :
    ! input variables :
    !     n    -> size of the vector
    !     vect -> the vector to integrate
    !     tabx -> abcissa of integration
    !     taby -> integration flux
    ! output variables :
    !     VECT_INT = result of the integral
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER,                        INTENT (IN) :: n
    REAL (KIND=dp), DIMENSION(1:n), INTENT (IN) :: vect  ! Vector to integrate
    REAL (KIND=dp), DIMENSION(1:n), INTENT (IN) :: tabx ! Abscissa
    REAL (KIND=dp), DIMENSION(1:n), INTENT (IN) :: taby ! Y-axis
    REAL (KIND=dp),                 INTENT (IN) :: wlim
    REAL (KIND=dp)                              :: dum1, dum2
    REAL (KIND=dp)                              :: a, b
    REAL (KIND=dp)                              :: sum
    INTEGER                                     :: i

    sum  = 0.0_dp
    DO i = 1, n-1
       IF      (tabx(i) > wlim) THEN
          EXIT
       ELSE IF (tabx(i+1) > wlim) THEN
          a    = vect(i) + (vect(i+1)-vect(i)) &
               * (wlim - tabx(i)) / (tabx(i+1) - tabx(i))
          b    = taby(i) + (taby(i+1)-taby(i)) &
               * (wlim - tabx(i)) / (tabx(i+1) - tabx(i))
          dum1 = vect(i) * taby(i)
          dum2 = a * b
          sum = sum + (dum1 + dum2) * (wlim - tabx(i)) / 2.0_dp
          
       ELSE
          dum1 = vect(i)   * taby(i)
          dum2 = vect(i+1) * taby(i+1)
          sum = sum + (dum1 + dum2) * (tabx(i+1) - tabx(i)) / 2.0_dp
       ENDIF
    ENDDO

    VECT_INT = sum

  END FUNCTION VECT_INT


  REAL(KIND=dp) FUNCTION  SECT_INT (sect)
    !---------------------------------------------------------------------------
    ! called by :
    !     DIFFUN
    ! purpose :
    !     integrate cross section with format CROSS_SEC over photon flux
    !     nrf from wl1 to wl2 where wl1 and wl2 are the wavelength limits 
    !     stored in the CROSS_SEC data.
    ! subroutine/function needed :
    ! input variables :
    !     sect -> the cross section to integrate
    ! output variables :
    ! results :
    !     SECT_INT = result of the integral
    !---------------------------------------------------------------------------
    USE MODULE_CHEM_REACT,       ONLY : CROSS_SEC

    IMPLICIT NONE

    TYPE(CROSS_SEC), INTENT (IN) :: sect             ! Section to integrate
 
    INTEGER                      :: npts
    INTEGER                      :: i_wl1, i_wl2
    REAL (KIND = dp)             :: wl1, wl2         ! wavelength limits for integration
    REAL (KIND = dp)             :: wla, wlb         ! wavelength steps for integration
    REAL (KIND = dp)             :: cross_a, cross_b ! cross sections at wla & wlb
    REAL (KIND = dp)             :: sum
    INTEGER                      :: i, j

    npts = sect%npts
    wl1 = sect%wl(1)
    wl2 = sect%wl(npts)
    IF(wl1 < wlth0) wl1 = wlth0
    IF(wl2 > wlthm) wl2 = wlthm
    CALL LOCATE_LEFT (wlg, nwlg, wl1, i_wl1)
    CALL LOCATE_LEFT (wlg, nwlg, wl2, i_wl2)

    sum = 0.0_dp
    j = 1
    DO i = i_wl1, i_wl2-1

       wla = wlg(i)
       wlb = wlg(i+1)

       IF(wla < sect%wl(j)) CYCLE
       ! find j which verifies 
       ! sect%wl(j) <= wla & sect%wl(j+1) > wla
       DO WHILE (sect%wl(j+1) < wla)
          j = j + 1
       ENDDO
       IF(sect%wl(j) > wla) THEN
          WRITE(*,*) "Error in the integration of a cross section"
          WRITE(*,*) "-> code stops"
          STOP
       ENDIF

       ! Compute the cross section at wla and wlb
       cross_a = sect%sigma(j) + (sect%sigma(j+1)-sect%sigma(j)) / (sect%wl(j+1) - sect%wl(j)) * (wla - sect%wl(j))
       cross_b = sect%sigma(j) + (sect%sigma(j+1)-sect%sigma(j)) / (sect%wl(j+1) - sect%wl(j)) * (wlb - sect%wl(j))

       ! perform integration
       sum = sum + (nrf(i) * cross_a + nrf(i+1) * cross_b) / 2.0_dp * (wlb - wla)
    ENDDO

    SECT_INT = sum

  END FUNCTION SECT_INT


  SUBROUTINE FGKCOEF
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     Computes line absorption using FGK (Federman, Glassgold, & Kwan 1979)
    !     approximation. fgkgrd contains radiation field attenuated by H2 lines
    !     assuming voigt profiles with a fixed turbulent velocity
    !     So far the routine computes FGK coefficient taking only into account
    !     absorption by the H2 gas buffer and not by the shock itself. Same goes
    !     for the radiation field which is attenuated by the dust buffer and not
    !     by the dust in the shock (see INIT_RADIATION)
    !     Note : including shock in the calculation of the FGK coefficients is
    !            tricky because shock line profile is not a Voigt profile with
    !            a turbulent velocity -> FGK approximation probably not adapted
    !                                 -> need to integrate line prof numerically
    !     Radiation field from the right is supposed to be 0
    !     optical depth on the right is supposed to be infinite
    ! subroutine/function needed :
    !     AJDF
    ! input variables :
    ! output variables :
    ! results :
    !     fgkgrd
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS,        ONLY : pi, sqpi, kB, qe, me, clum
    USE MODULE_TOOLS,            ONLY : MINIMUM, MAXIMUM
    USE MODULE_TECHCONFIG
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES, ONLY : speci, ind_H2, ind_CO
    USE MODULE_H2
    USE MODULE_CO

    IMPLICIT NONE

    !----------------------------------------------------
    ! Useful constants
    !----------------------------------------------------
    REAL (KIND=dp)                             :: ct1        ! absorption cross section (cm2)
    REAL (KIND=dp)                             :: vt2        ! vturb^2 (cm2 s-2)
    REAL (KIND=dp)                             :: bdop_h2    ! Doppler broadening (in A s-1)
    REAL (KIND=dp)                             :: bt_s_bd_h2 ! ratio of turbulent Doppler to total Doppler widths
    REAL (KIND=dp)                             :: bdop_co    ! Doppler broadening (in A s-1)
    REAL (KIND=dp)                             :: bt_s_bd_co ! ratio of turbulent Doppler to total Doppler widths

    !----------------------------------------------------
    ! Calculation of H2 and CO column densities
    !----------------------------------------------------
    REAL (KIND=dp), DIMENSION (1:NH2_lev)      :: cdh2_l  ! H2 column density (left  side)
    REAL (KIND=dp), DIMENSION (1:NH2_lev)      :: cdh2_r  ! H2 column density (right side)
    REAL (KIND=dp), DIMENSION (1:NCO_lev)      :: cdco_l  ! CO column density (left  side)
    REAL (KIND=dp), DIMENSION (1:NCO_lev)      :: cdco_r  ! CO column density (right side)

    !----------------------------------------------------
    ! line properties
    !----------------------------------------------------
    REAL (KIND=dp)                             :: xlam0   ! transition wavelength (A)
    INTEGER                                    :: nvl     ! vibrational quantum number V
    INTEGER                                    :: njl     ! rotational quantum number J
    INTEGER                                    :: lev     ! index for H2 levels
    INTEGER                                    :: jwl     ! index of the line in the wlg table
    INTEGER                                    :: i       ! dumy index

    !----------------------------------------------------
    ! Calculation of FGK coefficients
    !----------------------------------------------------
    REAL (KIND=dp), PARAMETER                  :: ff1 = 3.02e+00_dp ! constant of Eq. (A6) of FGK (1979)
    REAL (KIND=dp), PARAMETER                  :: ff2 = 1.00e+03_dp ! constant of Eq. (A6) of FGK (1979)
    REAL (KIND=dp), PARAMETER                  :: ff3 = 6.40e-02_dp ! constant of Eq. (A6) of FGK (1979)
    REAL (KIND=dp)                             :: gama       ! inverse lifetime / 4pi
    REAL (KIND=dp)                             :: betas      ! intermediary variable (same notation as FGK 1979)
    REAL (KIND=dp)                             :: betas0     ! intermediary variable (same notation as FGK 1979)
    REAL (KIND=dp)                             :: r          ! intermediary variable (same notation as FGK 1979)
    REAL (KIND=dp)                             :: t1         ! intermediary variable (same notation as FGK 1979)
    REAL (KIND=dp)                             :: opt0       ! intermediary variable
    REAL (KIND=dp)                             :: taud_l     ! tauD (same notation as FGK 1979) (left side)
    REAL (KIND=dp)                             :: taud_r     ! tauD (same notation as FGK 1979) (right side)
    REAL (KIND=dp)                             :: u1_l       ! intermediary variable (same notation as FGK 1979) (left side)
    REAL (KIND=dp)                             :: u1_r       ! intermediary variable (same notation as FGK 1979) (right side)
    REAL (KIND=dp)                             :: ajr_l      ! JR (same notation as FGK 1979) (left side)
    REAL (KIND=dp)                             :: ajr_r      ! JR (same notation as FGK 1979) (right side)
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: fgk_h2bx_l ! Line absorption (left side)
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: fgk_h2bx_r ! Line absorption (right side)
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: fgk_h2cx_l ! Line absorption (left side)
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: fgk_h2cx_r ! Line absorption (right side)
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: fgk_co_l   ! Line absorption (left side)
    REAL (KIND=dp), DIMENSION (:), ALLOCATABLE :: fgk_co_r   ! Line absorption (right side)

    !----------------------------------------------------
    ! radiation field
    !----------------------------------------------------
    REAL (KIND=dp), DIMENSION(1:nwlg)          :: rfu_l   ! left  side radiation field energy density (erg cm-3 A-1)
    REAL (KIND=dp), DIMENSION(1:nwlg)          :: rfu_r   ! right side radiation field energy density (erg cm-3 A-1)
    REAL (KIND=dp)                             :: rapwl   ! ratio of wavelength
    REAL (KIND=dp)                             :: aux1    ! dumy variable
    REAL (KIND=dp)                             :: aux2    ! dumy variable

    !===========================================================================
    ! 1 - Allocation and initialization
    !===========================================================================

    !----------------------------------------------------
    ! find the index of each transition in the wlg table
    !----------------------------------------------------
    DO i = 1, uvh2_bx%nuv
       xlam0 = uvh2_bx%wlg(i)
       CALL LOCATE_LEFT (wlg, nwlg, xlam0, jwl)
       uvh2_bx%iw(i) = jwl
    ENDDO
    DO i = 1, uvh2_cx%nuv
       xlam0 = uvh2_cx%wlg(i)
       CALL LOCATE_LEFT (wlg, nwlg, xlam0, jwl)
       uvh2_cx%iw(i) = jwl
    ENDDO
    DO i = 1, uvco%nuv
       xlam0 = uvco%wlg(i)
       CALL LOCATE_LEFT (wlg, nwlg, xlam0, jwl)
       uvco%iw(i) = jwl
    ENDDO

    !----------------------------------------------------
    ! allocate and initialize tables
    ! reset all fgk grid to their default values
    !----------------------------------------------------
    ALLOCATE ( fgk_h2bx_l(1:uvh2_bx%nuv) )
    ALLOCATE ( fgk_h2bx_r(1:uvh2_bx%nuv) )
    ALLOCATE ( fgk_h2cx_l(1:uvh2_cx%nuv) )
    ALLOCATE ( fgk_h2cx_r(1:uvh2_cx%nuv) )
    ALLOCATE ( fgk_co_l  (1:uvco%nuv   ) )
    ALLOCATE ( fgk_co_r  (1:uvco%nuv   ) )
    uvh2_bx%fgk(1:uvh2_bx%nuv) = 0.0_dp
    uvh2_cx%fgk(1:uvh2_cx%nuv) = 0.0_dp
    uvco%fgk   (1:uvco%nuv   ) = 0.0_dp
    fgk_h2bx_l (1:uvh2_bx%nuv) = 1.0_dp
    fgk_h2bx_r (1:uvh2_bx%nuv) = 1.0_dp
    fgk_h2cx_l (1:uvh2_cx%nuv) = 1.0_dp
    fgk_h2cx_r (1:uvh2_cx%nuv) = 1.0_dp
    fgk_co_l   (1:uvco%nuv   ) = 1.0_dp
    fgk_co_r   (1:uvco%nuv   ) = 1.0_dp

    !===========================================================================
    ! 2 - Definition of constants and computation of column dens (left & right)
    !===========================================================================

    !----------------------------------------------------
    ! Uselful constants
    !   - vt2 = vturb^2 square of turbulent velocity
    !   - bdop = Doppler width (converted in A)
    !   - bt_s_bd = ratio of turb to total Doppler widths
    !   - ct1 = pi e^2/mc = 2.654057e-2_dp (abs cross 
    !     section integrated over frequencies in cm2)
    ! BE CAREFUL : NOT PERFECT YET
    ! VTURB NEEDS TO BE REPLACE BY A PRESCRIPTION FOR THE
    ! SHOCK VELOCITY
    !----------------------------------------------------
    ct1     = pi * qe*qe / me / clum
    vt2     = vturb * vturb

    bdop_h2    = SQRT(vt2 + 2.0_dp * kB * Tn / speci(ind_H2)%mass)
    bt_s_bd_h2 = vturb / bdop_h2
    bdop_h2    = 1.0e8_dp * bdop_h2
    bdop_co    = SQRT(vt2 + 2.0_dp * kB * Tn / speci(ind_CO)%mass)
    bt_s_bd_co = vturb / bdop_co
    bdop_co    = 1.0e8_dp * bdop_co

    !----------------------------------------------------
    ! Computes Doppler width weighted H2 column density 
    ! for each level. We suppose that H2 column density 
    ! on the right is 10^23 (100 magnitudes of pure H2 
    ! extinction) -> shock only illuminated form behind.
    !----------------------------------------------------
    ! If S, H2 & CO level column dens in buffer computed
    ! - with a constant population distribution given by 
    !   the input o/p ratio if no H2 entry file and
    !   with a Boltzmann distribution for CO
    ! - with a constant population distribution given in
    !   the H2 entry file (ex. H2levels.in)
    ! If P, C, or J, H2 & CO level column dens computed
    ! in the buffer + in the structure taking into
    ! account variation of temperature, level population
    ! and microturbulence (but not velocity structure)
    !----------------------------------------------------
    DO lev = 1, NH2_lev
       cdh2_l(lev) = H2_lev(lev)%cd_l
       cdh2_r(lev) = H2_lev(lev)%cd_r
    ENDDO
    DO lev = 1, NCO_lev
       cdco_l(lev) = CO_lev(lev)%cd_l
       cdco_r(lev) = CO_lev(lev)%cd_r
    ENDDO

    !===========================================================================
    ! 3 - Calculation of FGK coefficients according to Eqs. in FGK (1979)
    !
    !     If there is a radiation field grid, FGK coefficient are set to 1 
    !      -> no need to compute opacities of H2 because the radiation field 
    !         is given by an external program
    !===========================================================================
    IF (F_COUP_RAD == 1) THEN

       DO i = 1, uvh2_bx%nuv
          !-------------------------------------------------
          ! Line properties
          !-------------------------------------------------
          nvl   = uvh2_bx%vl(i)
          njl   = uvh2_bx%jl(i)
          lev   = index_VJ_H2(nvl,njl)
          xlam0 = uvh2_bx%wlg(i)
          gama  = uvh2_bx%ilft(i) / (4.0_dp * pi)

          !-------------------------------------------------
          ! intermediary variables of FGK (1979)
          !-------------------------------------------------
          betas  = bdop_h2 / xlam0
          betas0 = 1.0e8_dp * vturb / xlam0
          r      = gama / (sqpi * betas)
          t1     = ff1 * (ff2 * r) ** (-ff3)
          opt0   = ct1 * uvh2_bx%osc(i) / (sqpi * betas0)

          !-------------------------------------------------
          ! Computation of tauD, JR, and JD (left and right)
          !-------------------------------------------------
          taud_l = opt0 * cdh2_l(lev)
          u1_l   = SQRT(taud_l * r) / t1
          ajr_l  = (r / t1) / SQRT(u1_l * u1_l + pi / 4.0_dp)
          taud_r = opt0 * cdh2_r(lev)
          u1_r   = SQRT(taud_r * r) / t1
          ajr_r  = (r / t1) / SQRT(u1_r * u1_r + pi / 4.0_dp)

          fgk_h2bx_l(i) = MINIMUM(AJDF(taud_l) + ajr_l, 1.0_dp)
          fgk_h2bx_r(i) = MINIMUM(AJDF(taud_r) + ajr_r, 1.0_dp)
       ENDDO

       DO i = 1, uvh2_cx%nuv
          !-------------------------------------------------
          ! Line properties
          !-------------------------------------------------
          nvl   = uvh2_cx%vl(i)
          njl   = uvh2_cx%jl(i)
          lev   = index_VJ_H2(nvl,njl)
          xlam0 = uvh2_cx%wlg(i)
          gama  = uvh2_cx%ilft(i) / (4.0_dp * pi)

          !-------------------------------------------------
          ! intermediary variables of FGK (1979)
          !-------------------------------------------------
          betas  = bdop_h2 / xlam0
          betas0 = 1.0e8_dp * vturb / xlam0
          r      = gama / (sqpi * betas)
          t1     = ff1 * (ff2 * r) ** (-ff3)
          opt0   = ct1 * uvh2_cx%osc(i) / (sqpi * betas0)

          !-------------------------------------------------
          ! Computation of tauD, JR, and JD (left and right)
          !-------------------------------------------------
          taud_l = opt0 * cdh2_l(lev)
          u1_l   = SQRT(taud_l * r) / t1
          ajr_l  = (r / t1) / SQRT(u1_l * u1_l + pi / 4.0_dp)
          taud_r = opt0 * cdh2_r(lev)
          u1_r   = SQRT(taud_r * r) / t1
          ajr_r  = (r / t1) / SQRT(u1_r * u1_r + pi / 4.0_dp)

          fgk_h2cx_l(i) = MINIMUM(AJDF(taud_l) + ajr_l, 1.0_dp)
          fgk_h2cx_r(i) = MINIMUM(AJDF(taud_r) + ajr_r, 1.0_dp)
       ENDDO

       DO i = 1, uvco%nuv
          !-------------------------------------------------
          ! Line properties
          !-------------------------------------------------
          nvl   = uvco%vl(i)
          njl   = uvco%jl(i)
          lev   = index_VJ_CO(nvl,njl)
          xlam0 = uvco%wlg(i)
          gama  = uvco%ilft(i) / (4.0_dp * pi)

          !-------------------------------------------------
          ! intermediary variables of FGK (1979)
          !-------------------------------------------------
          betas  = bdop_co / xlam0
          betas0 = 1.0e8_dp * vturb / xlam0
          r      = gama / (sqpi * betas)
          t1     = ff1 * (ff2 * r) ** (-ff3)
          opt0   = ct1 * uvco%osc(i) / (sqpi * betas0)

          !-------------------------------------------------
          ! Computation of tauD, JR, and JD (left and right)
          !-------------------------------------------------
          taud_l = opt0 * cdco_l(lev)
          u1_l   = SQRT(taud_l * r) / t1
          ajr_l  = (r / t1) / SQRT(u1_l * u1_l + pi / 4.0_dp)
          taud_r = opt0 * cdco_r(lev)
          u1_r   = SQRT(taud_r * r) / t1
          ajr_r  = (r / t1) / SQRT(u1_r * u1_r + pi / 4.0_dp)

          fgk_co_l(i) = MINIMUM(AJDF(taud_l) + ajr_l, 1.0_dp)
          fgk_co_r(i) = MINIMUM(AJDF(taud_r) + ajr_r, 1.0_dp)
       ENDDO

    ELSE IF (F_COUP_RAD == 2) THEN
       fgk_h2bx_l (1:uvh2_bx%nuv) = 1.0_dp
       fgk_h2bx_r (1:uvh2_bx%nuv) = 1.0_dp
       fgk_h2cx_l (1:uvh2_cx%nuv) = 1.0_dp
       fgk_h2cx_r (1:uvh2_cx%nuv) = 1.0_dp
       fgk_co_l   (1:uvco%nuv   ) = 1.0_dp
       fgk_co_r   (1:uvco%nuv   ) = 1.0_dp

    ENDIF

    ! ----------------------------------------------------------------------------------
    ! CHECK THE fgk_h2_l COEFFICIENTS (Only for debug purposes - BG 08/2016)
    ! ----------------------------------------------------------------------------------
    ! fichier = TRIM(out_dir)//'check_fgkrl_shock.res'
    ! OPEN(iwrtmp,file=TRIM(fichier),status='unknown')
    ! WRITE(iwrtmp,'(A8,ES10.3)') "vturb = ", vturb
    ! WRITE(iwrtmp,'(A8,ES10.3)') "Tn    = ", Tn
    ! WRITE(iwrtmp,'(A8,ES10.3)') "bdop  = ", bdop
    ! WRITE(iwrtmp,'(A8,ES10.3)') "NH2   = ", N_H2_0
    ! WRITE(iwrtmp,'(A8)')        "NH2*  = "
    ! DO i = 1, NH2_lev
    !    WRITE(iwrtmp,'(8X,ES10.3)') cdh2_l(i)
    ! ENDDO
    ! WRITE(*,*)
    ! DO i = 1, uvh2_bx%nuv
    !    WRITE(iwrtmp,'(2I4, ES10.3)') i, index_VJ_H2(uvh2_bx%vl(i),uvh2_bx%jl(i)), fgk_h2bx_l(i)
    ! ENDDO
    ! DO i = 1, uvh2_cx%nuv
    !    WRITE(iwrtmp,'(2I4, ES10.3)') i, index_VJ_H2(uvh2_cx%vl(i),uvh2_cx%jl(i)), fgk_h2cx_l(i)
    ! ENDDO
    ! CLOSE(iwrtmp)
    ! STOP
    ! ----------------------------------------------------------------------------------

    !===========================================================================
    ! 4 - Calculation of the radiation field in each line
    !===========================================================================

    !----------------------------------------------------
    ! Radiation field attenuation between the edge of the
    ! cloud and the present position
    !----------------------------------------------------
    rfu_l(1:nwlg) = urf(1:nwlg)
    rfu_r(1:nwlg) = 0.0_dp

    !----------------------------------------------------
    ! Radiation field in each transition
    !----------------------------------------------------
    DO i = 1, uvh2_bx%nuv
       xlam0 = uvh2_bx%wlg(i)
       jwl   = uvh2_bx%iw(i)
       IF (jwl < 1) THEN
          uvh2_bx%fgk(i) = 0.0_dp
       ELSE
          rapwl = (xlam0 - wlg(jwl)) / (wlg(jwl+1) - wlg(jwl))
          aux1 = rfu_l(jwl) + (rfu_l(jwl+1) - rfu_l(jwl)) * rapwl
          aux2 = rfu_r(jwl) + (rfu_r(jwl+1) - rfu_r(jwl)) * rapwl
          uvh2_bx%fgk(i) = aux1 * fgk_h2bx_l(i) + aux2 * fgk_h2bx_r(i)
       ENDIF
    ENDDO

    DO i = 1, uvh2_cx%nuv
       xlam0 = uvh2_cx%wlg(i)
       jwl   = uvh2_cx%iw(i)
       IF (jwl < 1) THEN
          uvh2_cx%fgk(i) = 0.0_dp
       ELSE
          rapwl = (xlam0 - wlg(jwl)) / (wlg(jwl+1) - wlg(jwl))
          aux1 = rfu_l(jwl) + (rfu_l(jwl+1) - rfu_l(jwl)) * rapwl
          aux2 = rfu_r(jwl) + (rfu_r(jwl+1) - rfu_r(jwl)) * rapwl
          uvh2_cx%fgk(i) = aux1 * fgk_h2cx_l(i) + aux2 * fgk_h2cx_r(i)
       ENDIF
    ENDDO

    DO i = 1, uvco%nuv
       xlam0 = uvco%wlg(i)
       jwl   = uvco%iw(i)
       IF (jwl < 1) THEN
          uvco%fgk(i) = 0.0_dp
       ELSE
          rapwl = (xlam0 - wlg(jwl)) / (wlg(jwl+1) - wlg(jwl))
          aux1 = rfu_l(jwl) + (rfu_l(jwl+1) - rfu_l(jwl)) * rapwl
          aux2 = rfu_r(jwl) + (rfu_r(jwl+1) - rfu_r(jwl)) * rapwl
          uvco%fgk(i) = aux1 * fgk_co_l(i) + aux2 * fgk_co_r(i)
       ENDIF
    ENDDO

    DEALLOCATE ( fgk_h2bx_l )
    DEALLOCATE ( fgk_h2bx_r )
    DEALLOCATE ( fgk_h2cx_l )
    DEALLOCATE ( fgk_h2cx_r )
    DEALLOCATE ( fgk_co_l   )
    DEALLOCATE ( fgk_co_r   )

  END SUBROUTINE FGKCOEF


  REAL (KIND=dp) FUNCTION AJDF(taud)
    !---------------------------------------------------------------------------
    ! called by :
    !     FGKCOEF
    ! purpose :
    ! subroutine/function needed :
    ! input variables :
    !     taud (opacity)
    ! output variables :
    ! results :
    !     AJDF -> JD at a given taud (FGK 1979, A8)
    !---------------------------------------------------------------------------
    IMPLICIT none

    REAL (KIND=dp), INTENT (IN) ::  taud

    ! IF      (taud < 2.0e+0_dp) THEN
    !    AJDF = EXP(-2.0_dp * taud / 3.0_dp)
    ! ELSE IF (taud < 1.0e+1_dp) THEN
    !    AJDF = 0.638_dp * taud**(-1.25_dp)
    ! ELSE IF (taud < 1.0e+2_dp) THEN
    !    AJDF = 0.505_dp * taud**(-1.15_dp)
    ! ELSE
    !    AJDF = 0.344_dp * taud**(-1.0667_dp)
    ! ENDIF

    ! IF      (taud < 2.0e+0_dp) THEN
    !    AJDF = EXP(-2.0_dp * taud / 3.0_dp)
    ! ELSE
    !    AJDF = EXP(-4.0_dp / 3.0_dp) * 2.0_dp**1.16_dp / taud**1.16_dp
    ! ENDIF

    ! AJDF = 1.0_dp/(1.0_dp+(taud/0.65_dp)**1.16_dp)
    AJDF = 1.0_dp/(1.0_dp+(taud/0.65_dp)**1.16_dp) + 0.15_dp * exp(-(log(taud)-log(0.8_dp))**2.0_dp)
  END FUNCTION AJDF


  SUBROUTINE FGKDATA
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    !     DIFFUN
    ! purpose :
    !     Prepares data for H2 detailed balance
    !     Computes H2 photodissociation rates (NOT YET BUT POSSIBLE OUTCOME)
    !     Computes populations of electronic levels
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     Pij_H2, sum_Pij_H2, h2bnua, h2cnua, probdissH2, probdissCO
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : pi, hp, qe, me, clum
    USE MODULE_PHYS_VAR
    USE MODULE_H2
    USE MODULE_CO

    IMPLICIT NONE

    !----------------------------------------------------
    ! Useful constants
    !----------------------------------------------------
    REAL (KIND=dp)                             :: coefb   ! constant used for conversion from f to Bij
    REAL (KIND=dp)                             :: vt2     ! vturb^2 (cm2 s-2)
    REAL (KIND=dp)                             :: sumh2   ! population of all H2 levels (ground elec state) (in cm-3)
    REAL (KIND=dp)                             :: sumco   ! population of all CO levels (ground elec state) (in cm-3)

    !----------------------------------------------------
    ! line properties
    !----------------------------------------------------
    REAL (KIND=dp)                             :: xlam0   ! transition wavelength (A)
    REAL (KIND=dp)                             :: xreh2   ! population of the lower level of H2 (in cm-3)
    REAL (KIND=dp)                             :: xreco   ! population of the lower level of co (in cm-3)
    INTEGER                                    :: nvu     ! vibrational quantum number V (upper level)
    INTEGER                                    :: nju     ! rotational  quantum number J (upper level)
    INTEGER                                    :: nvl     ! vibrational quantum number V (lower level)
    INTEGER                                    :: njl     ! rotational  quantum number J (lower level)
    INTEGER                                    :: ndj     ! delta J of the transitio
    INTEGER                                    :: lev     ! index for H2 levels
    INTEGER                                    :: i       ! dumy index
    INTEGER                                    :: j       ! dumy index

    !----------------------------------------------------
    ! einstein coefficients, dissociation rates
    !----------------------------------------------------
    REAL (KIND=dp)                             :: probdi  ! line dissociation probability
    REAL (KIND=dp)                             :: beinst  ! B Einstein coefficient
    REAL (KIND=dp)                             :: absorb  ! photon absorption rate
    REAL (KIND=dp)                             :: disici  ! line dissociation rate
    REAL (KIND=dp)                             :: dumh2   ! total H2 dissociation rate
    REAL (KIND=dp)                             :: dumco   ! total CO dissociation rate
    REAL (KIND=dp)                             :: pumping ! pumping rate (without dissociation)

    !===========================================================================
    ! 1 - Initialization
    !===========================================================================

    !----------------------------------------------------
    ! Numerical coefficient : 
    ! pi e**2 / (me * h * c**2) * 1.0e-16 (for Angstrom)
    !----------------------------------------------------
    coefb = 1.0e-16_dp * pi * qe*qe / (me * hp * clum * clum)
    vt2   = vturb * vturb
    sumh2 = SUM(H2_lev(1:NH2_lev)%density)
    sumco = SUM(CO_lev(1:NCO_lev)%density)

    !----------------------------------------------------
    ! Probabilities and electronic level populations
    !----------------------------------------------------
    Pij_H2(1:NH2_lev,1:NH2_lev) = 0.0_dp
    h2bnua(0:nvbc,0:njbc)       = 0.0_dp
    h2cnua(0:nvbc,0:njbc)       = 0.0_dp

    !----------------------------------------------------
    !  H2 dissociation probability
    !----------------------------------------------------
    dumh2             = 0.0_dp
    pdh2vJ(1:NH2_lev) = 0.0_dp

    !===========================================================================
    ! 2 - Treatment of H2 Lyman lines
    !===========================================================================

    DO i = 1, uvh2_bx%nuv
       !-------------------------------------------------
       ! Line properties
       !-------------------------------------------------
       nvl   = uvh2_bx%vl(i)
       njl   = uvh2_bx%jl(i)
       nvu   = uvh2_bx%vu(i)
       ndj   = uvh2_bx%dj(i)
       nju   = njl + ndj
       lev   = index_VJ_H2(nvl,njl)
       xlam0 = uvh2_bx%wlg(i)
       xreh2 = H2_lev(lev)%density / sumh2

       !-------------------------------------------------
       ! dissociation rates (s-1)
       !    - of each line ...... disici
       !    - from each level ... pdh2vJ(lev)
       !    - total ............. dumh2
       ! einstein coefficients, photon absorption rate
       !-------------------------------------------------
       probdi = uvh2_bx%pbdi(i)
       beinst = coefb * uvh2_bx%osc(i) * (xlam0**3)
       absorb = beinst * uvh2_bx%fgk(i)
       disici = absorb * probdi
       dumh2  = dumh2 + disici * xreh2
       pdh2vJ(lev) = pdh2vJ(lev) + disici

       !-------------------------------------------------
       ! Analysis - find the maximum contribution to
       !            H2 photodissociation
       !            dissociation probabilities are larger
       !            for large vu-vl transitions
       !-------------------------------------------------
       ! WRITE(*,*) "coucou", maxx, absorb * probdi * xreh2
       ! IF (absorb * probdi * xreh2 > maxx) THEN
       !    maxx = absorb * probdi * xreh2
       !    WRITE(*,*) nvl, njl, nvu, nju, maxx, absorb, probdi, xreh2
       ! ENDIF
       !-------------------------------------------------

       !-------------------------------------------------
       ! Calculation of band B populations (step 1)
       ! and of the radiative pumping matrix
       !    Suppress dissociation from equilibrium matrix
       !    (now taken care of in chemistry)
       !-------------------------------------------------
       IF(pumpH2 == 1) THEN
          IF ( (nvu <= nvbc) .AND. (nju <= njbc) ) THEN
             h2bnua(nvu,nju) = h2bnua(nvu,nju) + absorb * H2_lev(lev)%density
          ENDIF
          pumping = absorb * (1.0_dp - probdi)
          IF (MOD(njl,2) == 0) THEN
             Pij_H2(lev,levpc) = Pij_H2(lev,levpc) + pumping * bhmbx(nvu, nju, levpc)
          ELSE
             Pij_H2(lev,levic) = Pij_H2(lev,levic) + pumping * bhmbx(nvu, nju, levic)
          ENDIF
       ENDIF
    ENDDO

    !----------------------------------------------------
    ! Calculation of band B populations (final step)
    !----------------------------------------------------
    DO nvu = 0, nvbc
       DO nju = 0, njbc
          IF (tmunh2b(nvu,nju) > 1.0e-30_dp) THEN
             h2bnua(nvu,nju) = h2bnua(nvu,nju) / tmunh2b(nvu,nju)
          ELSE
             h2bnua(nvu,nju) = -1.0_dp
          ENDIF
       ENDDO
    ENDDO

    !===========================================================================
    ! 3 - Treatment of H2 Werner lines
    !===========================================================================

    DO i = 1, uvh2_cx%nuv
       !-------------------------------------------------
       ! Line properties
       !-------------------------------------------------
       nvl   = uvh2_cx%vl(i)
       njl   = uvh2_cx%jl(i)
       nvu   = uvh2_cx%vu(i)
       ndj   = uvh2_cx%dj(i)
       nju   = njl + ndj
       lev   = index_VJ_H2(nvl,njl)
       xlam0 = uvh2_cx%wlg(i)
       xreh2 = H2_lev(lev)%density / sumh2  

       !-------------------------------------------------
       ! dissociation rates (s-1)
       !    - of each line ...... disici
       !    - from each level ... pdh2vJ(lev)
       !    - total ............. dumh2
       ! einstein coefficients, photon absorption rate
       !-------------------------------------------------
       probdi = uvh2_cx%pbdi(i)
       beinst = coefb * uvh2_cx%osc(i) * (xlam0**3)
       absorb = beinst * uvh2_cx%fgk(i)
       disici = absorb * probdi
       dumh2  = dumh2 + disici * xreh2
       pdh2vJ(lev) = pdh2vJ(lev) + disici
  
       !-------------------------------------------------
       ! Calculation of band C populations (step 1)
       ! and of the radiative pumping matrix
       !    Suppress dissociation from equilibrium matrix
       !    (now taken care of in chemistry)
       !-------------------------------------------------
       IF(pumpH2 == 1) THEN
          IF ( (nvu <= nvbc) .AND. (nju <= njbc) ) THEN
             h2cnua(nvu,nju) = h2cnua(nvu,nju) + absorb * H2_lev(lev)%density
          ENDIF
          pumping = absorb * (1.0_dp - probdi)
          IF (MOD(njl,2) == 0) THEN
             Pij_H2(lev,levpc) = Pij_H2(lev,levpc) + pumping * bhmcx(nvu, nju, levpc)
          ELSE
             Pij_H2(lev,levic) = Pij_H2(lev,levic) + pumping * bhmcx(nvu, nju, levic)
          ENDIF
       ENDIF
    ENDDO

    !----------------------------------------------------
    ! Calculation of band C populations (final step)
    !----------------------------------------------------
    DO nvu = 0, nvbc
       DO nju = 0, njbc
          IF (tmunh2c(nvu,nju) > 1.0e-30_dp) THEN
             h2cnua(nvu,nju) = h2cnua(nvu,nju) / tmunh2c(nvu,nju)
          ELSE
             h2cnua(nvu,nju) = -1.0_dp
          ENDIF
       ENDDO
    ENDDO

    !----------------------------------------------------
    ! Levels deexcitation by pumping sum_Pij_H2
    !----------------------------------------------------
    DO j = 1, NH2_lev
       sum_Pij_H2(j) = SUM(Pij_H2(j,1:NH2_lev))
    END DO

    !----------------------------------------------------
    ! Calculation of elec level populations (in cm-3)
    !----------------------------------------------------
    DO nvu = 0, nvbc
       DO nju = 0, njbc
          IF (tmunh2b(nvu,nju) > 1.0e-30_dp) THEN
             h2elec = h2elec + h2bnua(nvu,nju)
          ENDIF
          IF (tmunh2c(nvu,nju) > 1.0e-30_dp) THEN
             h2elec = h2elec + h2cnua(nvu,nju)
          ENDIF
       ENDDO
    ENDDO

    !----------------------------------------------------
    ! Check coefficients and level populations
    ! - for debug purposes only
    !----------------------------------------------------
    ! DO j = 1, NH2_lev
    !    WRITE(*,'(I3, 2X, ES10.3)') j, sum_Pij_H2(j)
    ! END DO
    ! DO nvu = 0, nvbc
    !    DO nju = 0, njbc
    !       WRITE(*,'(I3, 2X, I3, ES10.3)') nvu, nju, h2cnua(nvu,nju)
    !    ENDDO
    ! ENDDO
    !----------------------------------------------------

    !----------------------------------------------------
    ! set the photodissociation probab of H2
    !----------------------------------------------------
    IF( F_COUP_RAD == 1 .OR. F_COUP_RAD == 2 ) THEN
       probdissH2 = dumh2
    ENDIF

    !===========================================================================
    ! 4 - Treatment of CO lines
    !===========================================================================
    dumco = 0.0_dp

    DO i = 1, uvco%nuv
       !-------------------------------------------------
       ! Line properties
       !-------------------------------------------------
       nvl   = uvco%vl(i)
       njl   = uvco%jl(i)
       lev   = index_VJ_CO(nvl,njl)
       xlam0 = uvco%wlg(i)
       xreco = CO_lev(lev)%density / sumco
       

       !-------------------------------------------------
       ! dissociation rates (s-1)
       !    - of each line ...... disici
       !    - total ............. dumco
       ! einstein coefficients, photon absorption rate
       !-------------------------------------------------
       probdi = uvco%pbdi(i)
       beinst = coefb * uvco%osc(i) * (xlam0**3)
       absorb = beinst * uvco%fgk(i)
       disici = absorb * probdi
       dumco  = dumco + disici * xreco
    ENDDO

    !----------------------------------------------------
    ! set the photodissociation probab of CO
    !----------------------------------------------------
    IF( F_COUP_RAD == 1 .OR. F_COUP_RAD == 2 ) THEN
       probdissCO = dumco
    ENDIF

  END SUBROUTINE FGKDATA



  REAL(KIND=dp) FUNCTION  HEATING_PHOTOCHEM_SECT_INT (sect)
    !---------------------------------------------------------------------------
    !****Created by Tabone 01/2018************
    ! called by :
    !     DIFFUN
    ! purpose :
    !     integrate cross section with format CROSS_SEC multipled by the difference
    !     btw energy of the ionizing/discociating photon and the ionization threshold
    !     given by the the first wt in the cross section file
    ! subroutine/function needed :
    ! input variables :
    !     sect -> the cross section to integrate
    ! output variables :
    ! results :
    !     SECT_INT = result of the integral
    !---------------------------------------------------------------------------
    USE MODULE_CHEM_REACT,       ONLY : CROSS_SEC
    USE MODULE_CONSTANTS

    IMPLICIT NONE

    TYPE(CROSS_SEC), INTENT (IN) :: sect             ! Section to integrate
 
    INTEGER                      :: npts
    INTEGER                      :: i_wl1, i_wl2
    REAL (KIND = dp)             :: wl1, wl2         ! wavelength limits for integration
    REAL (KIND = dp)             :: wla, wlb         ! wavelength steps for integration
    REAL (KIND = dp)             :: wlthphoto        ! wavelength of the ionization/photodestruction threshold Tabone 01/18
    REAL (KIND = dp)             :: cross_a, cross_b ! cross sections at wla & wlb
    REAL (KIND = dp)             :: sum
    INTEGER                      :: i, j

    npts      = sect%npts
    wl1       = sect%wl(1)      ! first wlth of the cross section
    wl2       = sect%wl(npts)   ! last wlth of the cross section
    wlthphoto = wl2             ! ionization/photodestruction threshold assumed to be the last element of the simga file (Tabone 01/18)
    IF(wl1 < wlth0) wl1 = wlth0
    IF(wl2 > wlthm) wl2 = wlthm
    CALL LOCATE_LEFT (wlg, nwlg, wl1, i_wl1) ! define index i_wl1 in the rf table of wl1
    CALL LOCATE_LEFT (wlg, nwlg, wl2, i_wl2) ! define index i_wl2 in the rf table of wl2

    sum = 0.0_dp
    j = 1
    DO i = i_wl1, i_wl2-1

       wla = wlg(i)
       wlb = wlg(i+1)

       IF(wla < sect%wl(j)) CYCLE
       ! find j which verifies 
       ! sect%wl(j) <= wla & sect%wl(j+1) > wla
       DO WHILE (sect%wl(j+1) < wla)
          j = j + 1
       ENDDO
       IF(sect%wl(j) > wla) THEN
          WRITE(*,*) "Error in the integration of a cross section"
          WRITE(*,*) "-> code stops"
          STOP
       ENDIF

       ! Compute the cross section at wla and wlb
       cross_a = sect%sigma(j) + (sect%sigma(j+1)-sect%sigma(j)) / (sect%wl(j+1) - sect%wl(j)) * (wla - sect%wl(j))
       cross_b = sect%sigma(j) + (sect%sigma(j+1)-sect%sigma(j)) / (sect%wl(j+1) - sect%wl(j)) * (wlb - sect%wl(j))

       ! perform integration of the photon flux multiplied by the difference in energy between ionizatin th
       ! and the energy of the photon
       sum = sum +  &
            (nrf(i)   * cross_a * (wlthphoto - wla) / wla + &
             nrf(i+1) * cross_b * (wlthphoto - wlb) / wlb) &
            / 2.0_dp * (wlb - wla)
       
    ENDDO
    sum = (clum*hp/(1.e-8*wlthphoto))*sum   ! in erg.s-1
    HEATING_PHOTOCHEM_SECT_INT = sum       

  END FUNCTION HEATING_PHOTOCHEM_SECT_INT


END MODULE MODULE_RADIATION
