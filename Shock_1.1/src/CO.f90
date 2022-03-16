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

MODULE MODULE_CO
  !*****************************************************************************
  !** The module 'MODULE_CO' contains variables and subroutines related to    **
  !** the CO molecule, except op_H2 which is in MODULE_PHYS_VAR               **
  !**     * levels and electronic transitions                                 **
  !** Note that this module is used to compute the photodissociation rate of  **
  !** CO and not the time dependent evolution of CO level populations         **
  !*****************************************************************************
  USE MODULE_TECHCONFIG
  USE MODULE_CHEMICAL_SPECIES, ONLY : TYPE_UVDATA
  use tex
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !-----------------------------------------
  ! data type of one level
  !-----------------------------------------
  TYPE TYPE_CO_LEVEL
     CHARACTER(len=lenIDLEV) :: Hname    ! Level human readable name
     CHARACTER(len=lenIDLEV) :: ID       ! Level ID
     INTEGER(KIND=LONG)      :: V, J     ! numbers of vibration and rotation
     REAL(KIND=DP)           :: Weight   ! statistical weight : 2J+1 (para) or 3(2J+1) (ortho)
     REAL(KIND=DP)           :: Energy   ! energy (K)
     REAL(KIND=DP)           :: Density  ! density of the level (cm-3)
     REAL(KIND=DP)           :: Density0 ! initial density of the level (cm-3) 
     REAL(KIND=DP)           :: Dens_old ! density at last call to DRIVE
     REAL(KIND=DP)           :: Col_dens ! column density of the level (cm-2)
     REAL(KIND=DP)           :: cd_l     ! left  side column density used for self-shielding (Doppler width weighted CO column density)
     REAL(KIND=DP)           :: cd_l_old ! left  side column density used for self-shielding (Doppler width weighted CO column density) at last call to DRIVE
     REAL(KIND=DP)           :: cd_r     ! right side column density used for self-shielding (Doppler width weighted CO column density)
  END TYPE TYPE_CO_LEVEL

  !-----------------------------------------
  ! energy levels
  ! variables defined in READ_CO_LEVELS
  !-----------------------------------------
  INTEGER                                          :: NCO_lev_max      ! maximal number of CO level
  INTEGER(KIND=LONG)                               :: Vmin_CO, Vmax_CO ! min. and max. values for V (vibration)
  INTEGER(KIND=LONG)                               :: Jmin_CO, Jmax_CO ! min. and max. values for J (rotation)

  TYPE(TYPE_CO_LEVEL), DIMENSION(:),   ALLOCATABLE :: CO_lev           ! vector containing the CO levels
  INTEGER(KIND=LONG),  DIMENSION(:,:), ALLOCATABLE :: index_VJ_CO      ! index (1..NCO_lev) of 1 level (V,J)

  INTEGER(KIND=LONG),  PARAMETER                   :: Num_undefined=-1 ! useful to initialize index_VJ_CO

  TYPE(type_uvdata)                                :: uvco             ! CO transitions

  REAL (KIND=DP)                                   :: probdissCO       ! total CO dissociation probability

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE READ_CO_LEVELS
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads CO levels : J, weight, energy.
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    CO_lev
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS,     ONLY : GET_FILE_NUMBER
    USE MODULE_CONSTANTS, ONLY : Zero
    IMPLICIT NONE

    CHARACTER(len=1)            :: charact
    INTEGER(KIND=LONG)          :: i, ii
    CHARACTER(LEN=*), PARAMETER :: name_file_CO_lev='input/CO_levels_Evj.in'
    INTEGER                     :: file_CO_lev

    INTEGER                     :: ios, whatever

    CHARACTER(LEN=5)            :: car_v, car_j

    ! Read the maximal number of CO levels
    file_CO_lev = GET_FILE_NUMBER()
    OPEN(file_CO_lev,file=name_file_CO_lev,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')
    DO i=1,6
       READ(file_CO_lev,*)
    END DO
    NCO_lev_max=0
    DO
       READ (file_CO_lev, *, iostat=ios) whatever
       IF (ios /= 0) EXIT
       NCO_lev_max = NCO_lev_max + 1
    ENDDO
    CLOSE(file_CO_lev)

    ! initialization
    ALLOCATE (CO_lev(NCO_lev_max))
    CO_lev(1:NCO_lev_max)%Hname    = ""
    CO_lev(1:NCO_lev_max)%V        = 0
    CO_lev(1:NCO_lev_max)%J        = 0
    CO_lev(1:NCO_lev_max)%weight   = Zero
    CO_lev(1:NCO_lev_max)%Energy   = Zero
    CO_lev(1:NCO_lev_max)%density  = Zero
    CO_lev(1:NCO_lev_max)%density0 = Zero
    CO_lev(1:NCO_lev_max)%Dens_old = Zero
    CO_lev(1:NCO_lev_max)%Col_dens = Zero

    ! file opening
    file_CO_lev = GET_FILE_NUMBER()
    OPEN(file_CO_lev,file=name_file_CO_lev,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    ! comments
    DO i=1,6
       READ(file_CO_lev,'(A1)')charact
    END DO

    DO i=1, NCO_lev_max
       READ(file_CO_lev,*) &
            ii, &
            CO_lev(i)%V, &
            CO_lev(i)%J, &
            CO_lev(i)%weight, &
            CO_lev(i)%Energy
            WRITE(car_v,'(I5)') CO_lev(i)%V
            WRITE(car_j,'(I5)') CO_lev(i)%J
            CO_lev(i)%Hname = "v=" // TRIM(ADJUSTL(car_v)) // "," // "J=" // TRIM(ADJUSTL(car_j))
            CO_lev(i)%ID = TRIM(ADJUSTL(car_v)) // "_" // TRIM(ADJUSTL(car_j))
       ! check if all lines have been correctly read, we must have i=ii
       IF (i /= ii) STOP "*** WARNING, error de read des niveaux de CO"
    END DO

    ! min. and max. values of V and J for these levels
    Vmin_CO=MINVAL(CO_lev(1:NCO_lev_max)%V)
    Vmax_CO=MAXVAL(CO_lev(1:NCO_lev_max)%V)
    Jmin_CO=MINVAL(CO_lev(1:NCO_lev_max)%J)
    Jmax_CO=MAXVAL(CO_lev(1:NCO_lev_max)%J)

    ! allocation and filling of the table index_VJ_CO
    ALLOCATE(index_VJ_CO(Vmin_CO:Vmax_CO, Jmin_CO:Jmax_CO))
    index_VJ_CO(Vmin_CO:Vmax_CO,Jmin_CO:Jmax_CO) = num_undefined
    DO i=1, NCO_lev_max
       index_VJ_CO(CO_lev(i)%V,CO_lev(i)%J) = i
    ENDDO

    ! file closure
    CLOSE(file_CO_lev)

  END SUBROUTINE READ_CO_LEVELS


  SUBROUTINE READ_CO_ELECTRONIC
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    ! purpose :
    !    reads CO electronic radiative transitions files
    !       - T112C16O.txt
    !         One transition per line
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    uvco%vl ..... v lower level
    !    uvco%jl ..... J lower level
    !    uvco%vu ..... v upper level (not used)
    !    uvco%dj ..... Delta J (J_up = J_low + Delta J) (not used)
    !    uvco%iw ..... index of line in wavelength grid
    !    uvco%osc .... oscillator strength
    !    uvco%wlg .... wavelength (Angstrom)
    !    uvco%ilft ... inverse life time (s-1)
    !    uvco%pbdi ... Upper level dissociation probability
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR,         ONLY : NCO_lev
    USE MODULE_TECHCONFIG,       ONLY : fichier
    USE MODULE_TOOLS,            ONLY : GET_FILE_NUMBER

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: name_file_uvco = 'input/UVdata/T112C16O.txt'

    INTEGER                     :: iread   ! file opening number
    INTEGER                     :: nvl     ! lower v level
    INTEGER                     :: njl     ! lower J level
    REAL (KIND=dp)              :: uv1     ! wavelength (in nm)
    REAL (KIND=dp)              :: uv2     ! wavenumber (in cm-1)
    REAL (KIND=dp)              :: uv3     ! oscillator strength
    REAL (KIND=dp)              :: uv4     ! upper state reciprocal lifetime (in s-1)
    REAL (KIND=dp)              :: uv5     ! upper state branching ratio towards dissociation

    INTEGER                     :: nuv     ! dummy integer
    INTEGER                     :: lev     ! dummy integer
    INTEGER                     :: i       ! dummy integer
    INTEGER                     :: j       ! dummy integer

    !----------------------------------------------------
    ! Prepare CO electronic transitions
    ! 1 - Determine the number of lines in file
    ! 2 - Find Vmin_bel, Jmin_bel, Vmax_bel, Jmax_bel
    !          Vmin_cel, Jmin_cel, Vmax_cel, Jmax_cel
    !          Vmin_el, Jmin_el, Vmax_el, & Jmax_el
    !----------------------------------------------------
    uvco%nuv = 0
    fichier = TRIM(name_file_uvco)
    iread   = GET_FILE_NUMBER()
    OPEN (iread, FILE = fichier, STATUS = 'old')
    READ (iread,'(/////////////////)')
    READ (iread,'(13X,I5)') nuv
    DO i = 1, nuv
       READ (iread,'(32X)',ADVANCE='no')
       READ (iread,*) nvl, njl, uv1, uv2, uv3, uv4, uv5
       lev = index_VJ_CO(nvl,njl)
       IF (lev > NCO_lev) CYCLE
       uvco%nuv = uvco%nuv + 1
    ENDDO
    CLOSE (iread)

    !----------------------------------------------------
    ! Allocate and initialize tables
    !----------------------------------------------------
    ALLOCATE ( uvco%vl  (uvco%nuv) )
    ALLOCATE ( uvco%jl  (uvco%nuv) )
    ALLOCATE ( uvco%vu  (uvco%nuv) )
    ALLOCATE ( uvco%dj  (uvco%nuv) )
    ALLOCATE ( uvco%iw  (uvco%nuv) )
    ALLOCATE ( uvco%osc (uvco%nuv) )
    ALLOCATE ( uvco%wlg (uvco%nuv) )
    ALLOCATE ( uvco%ilft(uvco%nuv) )
    ALLOCATE ( uvco%pbdi(uvco%nuv) )
    ALLOCATE ( uvco%fgk (uvco%nuv) )
    uvco%vl   = 0
    uvco%jl   = 0
    uvco%vu   = 0
    uvco%dj   = 0
    uvco%iw   = 0
    uvco%osc  = 0.0_dp
    uvco%wlg  = 0.0_dp
    uvco%ilft = 0.0_dp
    uvco%pbdi = 0.0_dp
    uvco%fgk  = 0.0_dp

    !----------------------------------------------------
    ! Read CO transitions
    !----------------------------------------------------
    fichier = TRIM(name_file_uvco)
    iread   = GET_FILE_NUMBER()
    OPEN (iread, FILE = fichier, STATUS = 'old')
    READ (iread,'(/////////////////)')
    READ (iread,'(13X,I5)') nuv
    j = 0
    DO i = 1, nuv
       READ (iread,'(32X)',ADVANCE='no')
       READ (iread,*) nvl, njl, uv1, uv2, uv3, uv4, uv5
       lev = index_VJ_CO(nvl,njl)
       IF (lev > NCO_lev) CYCLE
       j = j + 1
       uvco%vl(j)   = nvl
       uvco%jl(j)   = njl
       uvco%osc(j)  = uv3
       uvco%wlg(j)  = uv1 * 10.0_dp
       uvco%ilft(j) = uv4
       uvco%pbdi(j) = uv5
    ENDDO
    CLOSE (iread)
 
  END SUBROUTINE READ_CO_ELECTRONIC


  SUBROUTINE BOLTZMANN_ROVIB_CO
    !---------------------------------------------------------------------------
    ! called by :
    !    INITIALIZE
    !    DIFFUN
    ! purpose :
    !    compute the populations of the rovibrational levels of CO 
    !    assuming a Boltzmann distribution at the kinetic temperature 
    !    of the neutrals Tn
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !    CO_lev(*)%density
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR,         ONLY : NCO_lev, Tn
    USE MODULE_CHEMICAL_SPECIES, ONLY : speci, ind_co

    IMPLICIT NONE

    REAL (KIND=dp) :: dum ! dummy real
    INTEGER        :: i   ! dummy integer

    dum = 0.0
    DO i = 1, NCO_lev
       CO_lev(i)%density = CO_lev(i)%weight * TEXP(-(CO_lev(i)%Energy-CO_lev(index_VJ_CO(0,0))%Energy) / Tn)
       dum = dum + CO_lev(i)%density
    ENDDO
    CO_lev(1:NCO_lev)%density = CO_lev(1:NCO_lev)%density / dum * speci(ind_CO)%density

  END SUBROUTINE BOLTZMANN_ROVIB_CO

END MODULE MODULE_CO
