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

MODULE MODULE_SWITCH_CJ

  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


CONTAINS


  SUBROUTINE FIND_JUMP_RANGE(t_min_jump, t_max_jump, i_min_jump, i_max_jump, cs_active)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     find the range of potential jump conditions
    ! subroutine/function needed :
    ! input variables :
    !     counter -> the last iteration before reaching the sonic speed
    ! output variables :
    !     t_jump_min -> minimum distance for jump
    !     t_jump_max -> maximum distance for jump
    !     i_jump_min -> minimum index for jump
    !     i_jump_max -> maximum index for jump
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_TECHCONFIG
    USE MODULE_TOOLS,            ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR
    USE MODULE_PROFIL_TABLES
    USE MODULE_EVOLUTION

    IMPLICIT none

    REAL(KIND=DP), INTENT(out)   :: t_min_jump, t_max_jump
    INTEGER,       INTENT(out)   :: i_min_jump, i_max_jump
    LOGICAL,       INTENT(inout) :: cs_active
    LOGICAL                      :: allowed
    INTEGER                      :: i_jump
    REAL(KIND=DP)                :: dVn_old
    REAL(KIND=DP)                :: dVi_old
    REAL(KIND=DP)                :: Vn_old
    REAL(KIND=DP)                :: t_jump
    CHARACTER(LEN=lenfilename)   :: name_jump = 'jump_conditions.txt'
    INTEGER(KIND=LONG)           :: file_jump
    INTEGER                      :: i

    cs_active = .false.

    ! ----------------------------------------
    ! Find the maximal time to jump
    ! ----------------------------------------
    i_jump = counter - 1
    allowed = .false.
    DO WHILE ( .NOT.allowed .AND. i_jump > 2 )
       !--- reinitialize physical variables
       CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
       dVn_old = dVn
       dVi_old = dVi
       Vn_old  = Vn
       t_jump  = distance

       ! ------------------------------------
       ! Test the values pre and post jump
       ! Just for debug - not use otherwise
       ! ------------------------------------
       ! WRITE(*,'(I4,6(ES12.3,2X))',advance='no') i_jump, Vn, Vi, Tn, Ti, rhoN, rhoI
       ! ------------------------------------

       !--- Apply Rankine-Hugoniot relations for an adiabatic jump
       !--- hydrodynamical shock - on neutral fluid only
       CALL RH_JUMP_HYDRO

       ! ------------------------------------
       ! Test the values pre and post jump
       ! Just for debug - not use otherwise
       ! ------------------------------------
       ! WRITE(*,'(I4,6(ES12.3,2X))') i_jump, Vn, Vi, Tn, Ti, rhoN, rhoI
       ! ------------------------------------

       !--- call diffun for initial estimations
       distance_old = 0.0_dp
       CALL diffun(d_v_var,t_jump,v_lvariab,v_dvariab)
       dVn = v_dvariab(iv_Vn)
       dVi = v_dvariab(iv_Vi)

       ! The jump is too late
       IF( (dVn*dVn_old < 0.0_dp).OR.(Vn > Vn_old) ) THEN
          allowed = .false.
          t_max_jump = distance
          i_jump = i_jump - 1
       ELSE
          allowed = .true.
       ENDIF

       ! Identify potential C* shocks
       IF ( Vn > Vn_old ) THEN
          cs_active = .true.
       ENDIF

    ENDDO
    i_max_jump = i_jump

    ! ----------------------------------------
    ! Find the minimal time to jump
    ! ----------------------------------------
    i_jump = 2
    allowed = .false.
    DO WHILE ( .NOT.allowed .AND. i_jump < counter - 1 )
       !--- reinitialize physical variables
       CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
       dVn_old = dVn
       dVi_old = dVi
       t_jump  = distance

       ! ------------------------------------
       ! Test the values pre and post jump
       ! Just for debug - not use otherwise
       ! ------------------------------------
       ! WRITE(*,'(I4,6(ES12.3,2X))',advance='no') i_jump, Vn, Vi, Tn, Ti, rhoN, rhoI
       ! ------------------------------------

       !--- Apply Rankine-Hugoniot relations for an adiabatic jump
       !--- hydrodynamical shock - on neutral fluid only
       CALL RH_JUMP_HYDRO

       ! ------------------------------------
       ! Test the values pre and post jump
       ! Just for debug - not use otherwise
       ! ------------------------------------
       ! WRITE(*,'(I4,6(ES12.3,2X))') i_jump, Vn, Vi, Tn, Ti, rhoN, rhoI
       ! ------------------------------------

       !--- call diffun for initial estimations
       distance_old = 0.0_dp
       CALL diffun(d_v_var,t_jump,v_lvariab,v_dvariab)
       dVn = v_dvariab(iv_Vn)
       dVi = v_dvariab(iv_Vi)

       ! The jump is too soon
       IF ( ( Vn < Vi ).OR.( dVi*dVi_old < 0.0_dp ) ) THEN
          allowed = .false.
          t_min_jump = distance
          i_jump = i_jump + 1
       ELSE
          allowed = .true.
       ENDIF
    ENDDO
    i_min_jump = i_jump

    ! -------------------------------------------
    ! Print info on screen
    ! -------------------------------------------
    WRITE(*,*)
    WRITE(*,*) "range of potential jumping indices and time"
    WRITE(*,*) i_min_jump, i_max_jump, t_min_jump, t_max_jump
    WRITE(*,*)

    ! -------------------------------------------
    ! Adjust i_max for C* shocks if necessary,
    ! i.e. if the jump conditions are connected 
    ! to the the sonic point
    ! -------------------------------------------
    shock_type_fin  = "CJ"
    IF ( cs_active ) THEN
       CALL REINIT_PHYS_VARIABLES(traj_main(i_max_jump+1))
       Vn_old  = Vn
       CALL RH_JUMP_HYDRO
       IF ( Vn > Vn_old ) THEN
          shock_type_fin  = "C*"
          i_max_jump = i_max_jump + 1
          WRITE(*,*) "adjusted potential jumping indices and time"
          WRITE(*,*) i_min_jump, i_max_jump, t_min_jump, t_max_jump
          WRITE(*,*)
       ELSE
          cs_active = .false.
       ENDIF
    ENDIF

    ! -------------------------------------------
    ! Print potential jump conditions in file
    ! -------------------------------------------
    file_jump = GET_FILE_NUMBER()
    name_jump = TRIM(out_dir) // TRIM(name_jump)
    OPEN(file_jump,file=name_jump,status='UNKNOWN',access='SEQUENTIAL',form='FORMATTED',action='WRITE')
    WRITE(file_jump,'(A16,2X)',advance='no') '#     Vn (<jump)'
    WRITE(file_jump,'(A16,2X)',advance='no') '      Vi (<jump)'
    WRITE(file_jump,'(A16,2X)',advance='no') '      Vn (>jump)'
    WRITE(file_jump,'(A16,2X)',advance='no') '      Vi (>jump)'
    WRITE(file_jump,*)
    DO i = i_min_jump, i_max_jump
       CALL REINIT_PHYS_VARIABLES(traj_main(i))
       WRITE(file_jump,'(2(ES16.6,2X))',advance='no') Vn, Vi
       CALL RH_JUMP_HYDRO
       WRITE(file_jump,'(2(ES16.6,2X))') Vn, Vi
    ENDDO
    CLOSE(file_jump)

  END SUBROUTINE FIND_JUMP_RANGE


  SUBROUTINE FIND_RECOUPLING_JUMP(i_recoup_jump)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     find the index along the trajectory at which the neutrals
    !     jump exactly on the ion velocity
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    !     i_recoup_jump -> index of a recoupling jump
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_TECHCONFIG
    USE MODULE_PHYS_VAR
    USE MODULE_PROFIL_TABLES

    IMPLICIT none

    INTEGER,       INTENT(out)   :: i_recoup_jump
    INTEGER                      :: i_jump
    REAL(KIND=DP)                :: min

    min = Vs_cm
    DO i_jump = 2, counter-1
       CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
       CALL RH_JUMP_HYDRO
       IF( abs(Vn-Vi) < min ) THEN
           min = abs(Vn-Vi)
           i_recoup_jump = i_jump
       ENDIF
    ENDDO
    WRITE(*,*) "index recoupling jump = ", i_recoup_jump

  END SUBROUTINE FIND_RECOUPLING_JUMP


  SUBROUTINE FIND_CS_CROSS(i)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     find the index above which sound speed is crossed
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    !     i -> index along main trajectory above which sound speed is crossed
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR
    USE MODULE_PROFIL_TABLES

    IMPLICIT none

    INTEGER,       INTENT(out)   :: i

    ! ----------------------------------------
    ! Find the maximal time to jump
    ! ----------------------------------------
    i = 2
    CALL REINIT_PHYS_VARIABLES(traj_main(i))
    DO WHILE( Vn > Vsound )
       i = i + 1
       CALL REINIT_PHYS_VARIABLES(traj_main(i))
    ENDDO

  END SUBROUTINE FIND_CS_CROSS


  SUBROUTINE RH_JUMP_HYDRO
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    !     FIND_JUMP_RANGE
    ! purpose :
    !     apply the Rankine Hugoniot relations for
    !     an hydrodynamical adiabatic jump
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------

    USE MODULE_PHYS_VAR
    USE MODULE_GAMMA
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_H2
    USE MODULE_CO
    USE MODULE_CONSTANTS

    IMPLICIT none

    REAL(KIND=DP) :: mach
    REAL(KIND=DP) :: N1, N2
    REAL(KIND=DP) :: P1, P2
    REAL(KIND=DP) :: V1, V2
    REAL(KIND=DP) :: T1, T2
    REAL(KIND=DP) :: RHO1, RHO2
    REAL(KIND=DP) :: MU1, MU2
    INTEGER       :: i

    !---------------------------------------------------------------------------
    ! Initial jump conditions
    !---------------------------------------------------------------------------
    RHO1 = RhoN
    P1   = DensityN * kB * Tn
    V1   = Vn
    T1   = Tn
    MU1  = muN
    mach = Vn / Vsound
    N1   = DensityN

    !---------------------------------------------------------------------------
    ! Jump relations (p75 Guillet's thesis)
    !---------------------------------------------------------------------------
    P2   = ( 2_dp * GAMMA * mach**2.0_dp - (GAMMA - 1.0_dp) ) / (GAMMA + 1.0_dp) * P1
    RHO2 = (GAMMA + 1.0_dp) * mach**2.0_dp / ( (GAMMA - 1.0_dp) * mach**2.0_dp + 2.0_dp ) * RHO1
    V2   = V1 * RHO1 / RHO2
    MU2  = muN
    T2   = MU2 / MU1 * ( 2_dp * GAMMA * mach**2.0_dp - (GAMMA - 1.0_dp) ) &
         * ( (GAMMA - 1.0_dp) * mach**2.0_dp + 2.0_dp ) &
         / ( (GAMMA + 1.0_dp)**2.0_dp * mach**2.0_dp ) * T1
    N2   = P2 / (kB * T2)

    !---------------------------------------------------------------------------
    ! Implement results in neutral fluid
    !---------------------------------------------------------------------------

    ! --- hydrodynamics
    rhoN       = RHO2
    Tn         = T2
    DensityN   = N2
    muN        = MU2
    Vn         = V2
    compr_n    = compr_n * N2 / N1

    ! --- densities & H2 excitation
    DO i = b_neu, e_neu
       speci(i)%density = speci(i)%density * RHO2 / RHO1
    ENDDO
    nH = speci(ind_H)%density + 2._DP * speci(ind_H2)%density + speci(ind_Hplus)%density

    DO i = 1, NH2_lev
       H2_lev(i)%density = H2_lev(i)%density * RHO2 / RHO1
    ENDDO

    DO i = 1, NCO_lev
       CO_lev(i)%density = CO_lev(i)%density * RHO2 / RHO1
    ENDDO

    ! --- dvode tables
    v_variab(iv_Vn        ) = Vn
    v_variab(iv_vCH       ) = Vn
    v_variab(iv_vS        ) = Vn
    v_variab(iv_vSH       ) = Vn
    v_variab(iv_RhoN      ) = RhoN
    v_variab(iv_Tn        ) = Tn
    v_variab(iv_DensityN  ) = DensityN
    v_variab(iv_compr_n  )  = compr_n
    v_variab(bv_neu:ev_neu) = speci(b_neu:e_neu)%density
    v_variab(bv_H2_lev:ev_H2_lev) = H2_lev(1:NH2_lev_var)%density
    WHERE (v_variab > 0)
       v_lvariab = LOG(v_variab)
    ELSEWHERE
       v_lvariab = minus_infinity
    END WHERE

  END SUBROUTINE RH_JUMP_HYDRO


  SUBROUTINE INTERP_JUMP_CONDITIONS(t_jump,i_jump)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     interpolate jump conditions at an input time (distance) between
    !     two points previously calculated
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_RADIATION
    USE MODULE_H2
    USE MODULE_CO
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_PHYS_VAR
    USE MODULE_GRAINS
    USE MODULE_PROFIL_TABLES

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(in)                     :: t_jump
    INTEGER,       INTENT(in)                     :: i_jump

    REAL(KIND=DP)                                 :: distance_1   , distance_2
    REAL(KIND=DP)                                 :: Av_1         , Av_2
    REAL(KIND=DP)                                 :: timeN_1      , timeN_2
    REAL(KIND=DP)                                 :: timeI_1      , timeI_2
    REAL(KIND=DP)                                 :: nH_1         , nH_2
    REAL(KIND=DP)                                 :: rhoN_1       , rhoN_2
    REAL(KIND=DP)                                 :: rhoI_1       , rhoI_2
    REAL(KIND=DP)                                 :: rhoA_1       , rhoA_2
    REAL(KIND=DP)                                 :: rhoNeg_1     , rhoNeg_2
    REAL(KIND=DP)                                 :: DensityN_1   , DensityN_2
    REAL(KIND=DP)                                 :: DensityI_1   , DensityI_2
    REAL(KIND=DP)                                 :: DensityA_1   , DensityA_2
    REAL(KIND=DP)                                 :: DensityNeg_1 , DensityNeg_2
    REAL(KIND=DP)                                 :: muN_1        , muN_2
    REAL(KIND=DP)                                 :: muI_1        , muI_2
    REAL(KIND=DP)                                 :: muA_1        , muA_2
    REAL(KIND=DP)                                 :: muNeg_1      , muNeg_2
    REAL(KIND=DP)                                 :: Vsound_1     , Vsound_2
    REAL(KIND=DP)                                 :: Vmagnet_1    , Vmagnet_2
    REAL(KIND=DP)                                 :: Valfven_1    , Valfven_2
    REAL(KIND=DP)                                 :: Vn_1         , Vn_2
    REAL(KIND=DP)                                 :: Vi_1         , Vi_2
    REAL(KIND=DP)                                 :: dVn_1        , dVn_2
    REAL(KIND=DP)                                 :: dVi_1        , dVi_2
    REAL(KIND=DP)                                 :: Vgrad_1      , Vgrad_2
    REAL(KIND=DP)                                 :: grad_V_1     , grad_V_2
    REAL(KIND=DP)                                 :: v_CH_1       , v_CH_2
    REAL(KIND=DP)                                 :: v_S_1        , v_S_2
    REAL(KIND=DP)                                 :: v_SH_1       , v_SH_2
    REAL(KIND=DP)                                 :: Tn_1         , Tn_2
    REAL(KIND=DP)                                 :: Ti_1         , Ti_2
    REAL(KIND=DP)                                 :: Te_1         , Te_2
    REAL(KIND=DP)                                 :: Tgrain_1     , Tgrain_2
    REAL(KIND=DP)                                 :: Teff_grain_1 , Teff_grain_2
    REAL(KIND=DP)                                 :: op_H2_1      , op_H2_2
    REAL(KIND=DP)                                 :: compr_n_1    , compr_n_2
    REAL(KIND=DP)                                 :: compr_i_1    , compr_i_2
    REAL(KIND=DP)                                 :: fluph_1      , fluph_2
    REAL(KIND=DP)                                 :: coldens_h_1  , coldens_h_2
    REAL(KIND=DP)                                 :: coldens_h2_1 , coldens_h2_2
    REAL(KIND=DP)                                 :: coldens_co_1 , coldens_co_2
    REAL(KIND=DP)                                 :: Col_Dens_nH_1, Col_Dens_nH_2
    REAL(KIND=DP),      DIMENSION(nangle,1:nwlg)  :: irf_1        , irf_2
    ! REAL(KIND=DP),      DIMENSION(1:nwlg)         :: optdpth_1    , optdpth_2
    TYPE(TYPE_H2_LEVEL),DIMENSION(:), ALLOCATABLE :: H2_lev_1     , H2_lev_2
    TYPE(TYPE_H2_LEVEL),DIMENSION(:), ALLOCATABLE :: CO_lev_1     , CO_lev_2
    TYPE(type_specy),   DIMENSION(:), ALLOCATABLE :: speci_1      , speci_2

    REAL(KIND=DP)                                 :: dum
    INTEGER                                       :: i

    ! -------------------------------------------
    ! Allocate necessary tables
    ! -------------------------------------------
    ALLOCATE (H2_lev_1(NH2_lev))
    ALLOCATE (H2_lev_2(NH2_lev))
    ALLOCATE (CO_lev_1(NCO_lev))
    ALLOCATE (CO_lev_2(NCO_lev))
    ALLOCATE (speci_1(0:Nspec_plus))
    ALLOCATE (speci_2(0:Nspec_plus))

    ! -------------------------------------------
    ! Conditions at first point
    ! -------------------------------------------
    CALL REINIT_PHYS_VARIABLES(traj_main(i_jump))
    distance_1   = distance
    Av_1         = Av
    timeN_1      = timeN
    timeI_1      = timeI
    nH_1         = nH
    rhoN_1       = rhoN
    rhoI_1       = rhoI
    rhoA_1       = rhoA
    rhoNeg_1     = rhoNeg
    DensityN_1   = densityN
    DensityI_1   = densityI
    DensityA_1   = densityA
    DensityNeg_1 = densityNeg
    muN_1        = muN
    muI_1        = muI
    muA_1        = muA
    muNeg_1      = muNeg
    Vsound_1     = Vsound
    Vmagnet_1    = Vmagnet
    Valfven_1    = Valfven
    Vn_1         = Vn
    Vi_1         = Vi
    dVn_1        = dVn
    dVi_1        = dVi
    Vgrad_1      = Vgrad
    grad_V_1     = grad_V
    v_CH_1       = v_CH
    v_S_1        = V_S
    v_SH_1       = V_SH
    Tn_1         = Tn
    Ti_1         = Ti
    Te_1         = Te
    Tgrain_1     = Tgrain
    Teff_grain_1 = Teff_grain
    op_H2_1      = op_H2
    compr_n_1    = compr_n
    compr_i_1    = compr_i
    fluph_1      = fluph
    coldens_h_1  = coldens_h
    coldens_h2_1 = coldens_h2
    coldens_co_1 = coldens_co
    Col_Dens_nH_1= Col_Dens_nH
    IF (F_AV /= 0) THEN
       DO i = 1, nangle
          irf_1(i,1:nwlg) = irf(i,1:nwlg)
       ENDDO
       ! optdpth_1(1:nwlg) = optdpth(1:nwlg)
    ENDIF
    DO i = 1, NH2_lev
       H2_lev_1(i)%cd_l = H2_lev(i)%cd_l
    ENDDO
    DO i = 1, NCO_lev
       CO_lev_1(i)%cd_l = CO_lev(i)%cd_l
    ENDDO
    DO i = 1, Nspec
       speci_1(i)%density  = speci(i)%density 
       speci_1(i)%Col_dens = speci(i)%Col_dens
    ENDDO
    DO i = 1, NH2_lev
       H2_lev_1(i)%density  = H2_lev(i)%density
       H2_lev_1(i)%Col_dens = H2_lev(i)%Col_dens
    ENDDO
    DO i = 1, NCO_lev
       CO_lev_1(i)%density  = CO_lev(i)%density
       CO_lev_1(i)%Col_dens = CO_lev(i)%Col_dens
    ENDDO

    ! -------------------------------------------
    ! Conditions at second point
    ! -------------------------------------------
    CALL REINIT_PHYS_VARIABLES(traj_main(i_jump+1))
    distance_2   = distance
    Av_2         = Av
    timeN_2      = timeN
    timeI_2      = timeI
    nH_2         = nH
    rhoN_2       = rhoN
    rhoI_2       = rhoI
    rhoA_2       = rhoA
    rhoNeg_2     = rhoNeg
    DensityN_2   = densityN
    DensityI_2   = densityI
    DensityA_2   = densityA
    DensityNeg_2 = densityNeg
    muN_2        = muN
    muI_2        = muI
    muA_2        = muA
    muNeg_2      = muNeg
    Vsound_2     = Vsound
    Vmagnet_2    = Vmagnet
    Valfven_2    = Valfven
    Vn_2         = Vn
    Vi_2         = Vi
    dVn_2        = dVn
    dVi_2        = dVi
    Vgrad_2      = Vgrad
    grad_V_2     = grad_V
    v_CH_2       = v_CH
    v_S_2        = V_S
    v_SH_2       = V_SH
    Tn_2         = Tn
    Ti_2         = Ti
    Te_2         = Te
    Tgrain_2     = Tgrain
    Teff_grain_2 = Teff_grain
    op_H2_2      = op_H2
    compr_n_2    = compr_n
    compr_i_2    = compr_i
    fluph_2      = fluph
    coldens_h_2  = coldens_h
    coldens_h2_2 = coldens_h2
    coldens_co_2 = coldens_co
    Col_Dens_nH_2= Col_Dens_nH
    IF (F_AV /= 0) THEN
       DO i = 1, nangle
          irf_2(i,1:nwlg) = irf(i,1:nwlg)
       ENDDO
       ! optdpth_2(1:nwlg) = optdpth(1:nwlg)
    ENDIF
    DO i = 1, NH2_lev
       H2_lev_2(i)%cd_l = H2_lev(i)%cd_l
    ENDDO
    DO i = 1, NCO_lev
       CO_lev_2(i)%cd_l = CO_lev(i)%cd_l
    ENDDO
    DO i = 1, Nspec
       speci_2(i)%density  = speci(i)%density 
       speci_2(i)%Col_dens = speci(i)%Col_dens
    ENDDO
    DO i = 1, NH2_lev
       H2_lev_2(i)%density  = H2_lev(i)%density
       H2_lev_2(i)%Col_dens = H2_lev(i)%Col_dens
    ENDDO
    DO i = 1, NCO_lev
       CO_lev_2(i)%density  = CO_lev(i)%density
       CO_lev_2(i)%Col_dens = CO_lev(i)%Col_dens
    ENDDO

    ! -------------------------------------------
    ! Logarithmic interpolation
    ! -------------------------------------------
    dum = ( LOG(t_jump) - LOG(distance_1) ) / ( LOG(distance_2) - LOG(distance_1) )

    distance    = LOG(ABS(distance_1   )) + ( LOG(ABS(distance_2  )) - LOG(ABS(distance_1   )) ) * dum
    ! Av          = LOG(ABS(Av_1         )) + ( LOG(ABS(Av_2        )) - LOG(ABS(Av_1         )) ) * dum
    ! timeN       = LOG(ABS(timeN_1      )) + ( LOG(ABS(timeN_2     )) - LOG(ABS(timeN_1      )) ) * dum
    ! timeI       = LOG(ABS(timeI_1      )) + ( LOG(ABS(timeI_2     )) - LOG(ABS(timeI_1      )) ) * dum
    ! nH          = LOG(ABS(nH_1         )) + ( LOG(ABS(nH_2        )) - LOG(ABS(nH_1         )) ) * dum
    rhoN        = LOG(ABS(rhoN_1       )) + ( LOG(ABS(rhoN_2      )) - LOG(ABS(rhoN_1       )) ) * dum
    ! rhoI        = LOG(ABS(rhoI_1       )) + ( LOG(ABS(rhoI_2      )) - LOG(ABS(rhoI_1       )) ) * dum
    ! rhoA        = LOG(ABS(rhoA_1       )) + ( LOG(ABS(rhoA_2      )) - LOG(ABS(rhoA_1       )) ) * dum
    ! rhoNeg      = LOG(ABS(rhoNeg_1     )) + ( LOG(ABS(rhoNeg_2    )) - LOG(ABS(rhoNeg_1     )) ) * dum
    ! densityN    = LOG(ABS(DensityN_1   )) + ( LOG(ABS(DensityN_2  )) - LOG(ABS(DensityN_1   )) ) * dum
    ! densityI    = LOG(ABS(DensityI_1   )) + ( LOG(ABS(DensityI_2  )) - LOG(ABS(DensityI_1   )) ) * dum
    ! densityA    = LOG(ABS(DensityA_1   )) + ( LOG(ABS(DensityA_2  )) - LOG(ABS(DensityA_1   )) ) * dum
    ! densityNeg  = LOG(ABS(DensityNeg_1 )) + ( LOG(ABS(DensityNeg_2)) - LOG(ABS(DensityNeg_1 )) ) * dum
    ! muN         = LOG(ABS(muN_1        )) + ( LOG(ABS(muN_2       )) - LOG(ABS(muN_1        )) ) * dum
    ! muI         = LOG(ABS(muI_1        )) + ( LOG(ABS(muI_2       )) - LOG(ABS(muI_1        )) ) * dum
    ! muA         = LOG(ABS(muA_1        )) + ( LOG(ABS(muA_2       )) - LOG(ABS(muA_1        )) ) * dum
    ! muNeg       = LOG(ABS(muNeg_1      )) + ( LOG(ABS(muNeg_2     )) - LOG(ABS(muNeg_1      )) ) * dum
    Vsound      = LOG(ABS(Vsound_1     )) + ( LOG(ABS(Vsound_2    )) - LOG(ABS(Vsound_1     )) ) * dum
    ! Vmagnet     = LOG(ABS(Vmagnet_1    )) + ( LOG(ABS(Vmagnet_2   )) - LOG(ABS(Vmagnet_1    )) ) * dum
    ! Valfven     = LOG(ABS(Valfven_1    )) + ( LOG(ABS(Valfven_2   )) - LOG(ABS(Valfven_1    )) ) * dum
    Vn          = LOG(ABS(Vn_1         )) + ( LOG(ABS(Vn_2        )) - LOG(ABS(Vn_1         )) ) * dum
    Vi          = LOG(ABS(Vi_1         )) + ( LOG(ABS(Vi_2        )) - LOG(ABS(Vi_1         )) ) * dum
    ! dVn         = LOG(ABS(dVn_1        )) + ( LOG(ABS(dVn_2       )) - LOG(ABS(dVn_1        )) ) * dum
    ! dVi         = LOG(ABS(dVi_1        )) + ( LOG(ABS(dVi_2       )) - LOG(ABS(dVi_1        )) ) * dum
    ! Vgrad       = LOG(ABS(Vgrad_1      )) + ( LOG(ABS(Vgrad_2     )) - LOG(ABS(Vgrad_1      )) ) * dum
    ! grad_V      = LOG(ABS(grad_V_1     )) + ( LOG(ABS(grad_V_2    )) - LOG(ABS(grad_V_1     )) ) * dum
    ! v_CH        = LOG(ABS(v_CH_1       )) + ( LOG(ABS(v_CH_2      )) - LOG(ABS(v_CH_1       )) ) * dum
    ! V_S         = LOG(ABS(v_S_1        )) + ( LOG(ABS(v_S_2       )) - LOG(ABS(v_S_1        )) ) * dum
    ! V_SH        = LOG(ABS(v_SH_1       )) + ( LOG(ABS(v_SH_2      )) - LOG(ABS(v_SH_1       )) ) * dum
    Tn          = LOG(ABS(Tn_1         )) + ( LOG(ABS(Tn_2        )) - LOG(ABS(Tn_1         )) ) * dum
    Ti          = LOG(ABS(Ti_1         )) + ( LOG(ABS(Ti_2        )) - LOG(ABS(Ti_1         )) ) * dum
    Te          = LOG(ABS(Te_1         )) + ( LOG(ABS(Te_2        )) - LOG(ABS(Te_1         )) ) * dum
    ! Tgrain      = LOG(ABS(Tgrain_1     )) + ( LOG(ABS(Tgrain_2     )) - LOG(ABS(Tgrain_1     )) ) * dum
    ! Teff_grain  = LOG(ABS(Teff_grain_1 )) + ( LOG(ABS(Teff_grain_2 )) - LOG(ABS(Teff_grain_1 )) ) * dum
    ! op_H2       = LOG(ABS(op_H2_1      )) + ( LOG(ABS(op_H2_2      )) - LOG(ABS(op_H2_1      )) ) * dum
    ! compr_n     = LOG(ABS(compr_n_1    )) + ( LOG(ABS(compr_n_2   )) - LOG(ABS(compr_n_1    )) ) * dum
    ! compr_i     = LOG(ABS(compr_i_1    )) + ( LOG(ABS(compr_i_2   )) - LOG(ABS(compr_i_1    )) ) * dum
    ! fluph       = LOG(ABS(fluph_1      )) + ( LOG(ABS(fluph_2      )) - LOG(ABS(fluph_1      )) ) * dum
    ! coldens_h   = LOG(ABS(coldens_h_1  )) + ( LOG(ABS(coldens_h_2  )) - LOG(ABS(coldens_h_1  )) ) * dum
    ! coldens_h2  = LOG(ABS(coldens_h2_1 )) + ( LOG(ABS(coldens_h2_2 )) - LOG(ABS(coldens_h2_1 )) ) * dum
    ! coldens_co  = LOG(ABS(coldens_co_1 )) + ( LOG(ABS(coldens_co_2 )) - LOG(ABS(coldens_co_1 )) ) * dum
    ! Col_Dens_nH = LOG(ABS(Col_Dens_nH_1)) + ( LOG(ABS(Col_Dens_nH_2)) - LOG(ABS(Col_Dens_nH_1)) ) * dum
    ! IF (F_AV /= 0) THEN
    !    DO i = 1, nangle
    !       irf(i,1:nwlg)    = LOG(irf_1(i,1:nwlg)     ) + ( LOG(irf_2(i,1:nwlg)    ) - LOG(irf_1(i,1:nwlg)    ) ) * dum
    !    ENDDO
    ! ENDIF
    ! DO i = 1, NH2_lev
    !    H2_lev(i)%cd_l      = LOG(H2_lev_1(i)%cd_l    ) + ( LOG(H2_lev_2(i)%cd_l   ) - LOG(H2_lev_1(i)%cd_l   ) ) * dum
    ! ENDDO
    ! DO i = 1, Nspec
    !    speci(i)%density    = LOG(speci_1(i)%density  ) + ( LOG(speci_2(i)%density ) - LOG(speci_1(i)%density ) ) * dum
    !    speci(i)%Col_dens   = LOG(speci_1(i)%Col_dens ) + ( LOG(speci_2(i)%Col_dens) - LOG(speci_1(i)%Col_dens) ) * dum
    ! ENDDO
    ! DO i = 1, NH2_lev
    !    H2_lev(i)%density   = LOG(H2_lev_1(i)%density ) + ( LOG(H2_lev_2(i)%density ) - LOG(H2_lev_1(i)%density ) ) * dum
    !    H2_lev(i)%Col_dens  = LOG(H2_lev_1(i)%Col_dens) + ( LOG(H2_lev_2(i)%Col_dens) - LOG(H2_lev_1(i)%Col_dens) ) * dum
    ! ENDDO
    
    ! -------------------------------------------
    ! Exponential development
    ! -------------------------------------------
    distance    = SIGN( EXP(distance   ) , distance_1   )
    ! Av          = SIGN( EXP(Av         ) , Av_1         )
    ! timeN       = SIGN( EXP(timeN      ) , timeN_1      )
    ! timeI       = SIGN( EXP(timeI      ) , timeI_1      )
    ! nH          = SIGN( EXP(nH         ) , nH_1         )
    rhoN        = SIGN( EXP(rhoN       ) , rhoN_1       )
    ! rhoI        = SIGN( EXP(rhoI       ) , rhoI_1       )
    ! rhoA        = SIGN( EXP(rhoA       ) , rhoA_1       )
    ! rhoNeg      = SIGN( EXP(rhoNeg     ) , rhoNeg_1     )
    ! densityN    = SIGN( EXP(densityN   ) , DensityN_1   )
    ! densityI    = SIGN( EXP(densityI   ) , DensityI_1   )
    ! densityA    = SIGN( EXP(densityA   ) , DensityA_1   )
    ! densityNeg  = SIGN( EXP(densityNeg ) , DensityNeg_1 )
    ! muN         = SIGN( EXP(muN        ) , muN_1        )
    ! muI         = SIGN( EXP(muI        ) , muI_1        )
    ! muA         = SIGN( EXP(muA        ) , muA_1        )
    ! muNeg       = SIGN( EXP(muNeg      ) , muNeg_1      )
    Vsound      = SIGN( EXP(Vsound     ) , Vsound_1     )
    ! Vmagnet     = SIGN( EXP(Vmagnet    ) , Vmagnet_1    )
    ! Valfven     = SIGN( EXP(Valfven    ) , Valfven_1    )
    Vn          = SIGN( EXP(Vn         ) , Vn_1         )
    Vi          = SIGN( EXP(Vi         ) , Vi_1         )
    ! dVn         = SIGN( EXP(dVn        ) , dVn_1        )
    ! dVi         = SIGN( EXP(dVi        ) , dVi_1        )
    ! Vgrad       = SIGN( EXP(Vgrad      ) , Vgrad_1      )
    ! grad_V      = SIGN( EXP(grad_V     ) , grad_V_1     )
    ! v_CH        = SIGN( EXP(v_CH       ) , v_CH_1       )
    ! V_S         = SIGN( EXP(V_S        ) , v_S_1        )
    ! V_SH        = SIGN( EXP(V_SH       ) , v_SH_1       )
    Tn          = SIGN( EXP(Tn         ) , Tn_1         )
    Ti          = SIGN( EXP(Ti         ) , Ti_1         )
    Te          = SIGN( EXP(Te         ) , Te_1         )
    ! Tgrain      = SIGN( EXP(Tgrain     ) , Tgrain_1     )
    ! Teff_grain  = SIGN( EXP(Teff_grain ) , Teff_grain_1 )
    ! op_H2       = SIGN( EXP(op_H2      ) , op_H2_1      )
    ! compr_n     = SIGN( EXP(compr_n    ) , compr_n_1    )
    ! compr_i     = SIGN( EXP(compr_i    ) , compr_i_1    )
    ! fluph       = SIGN( EXP(fluph      ) , fluph_1      )
    ! coldens_h   = SIGN( EXP(coldens_h  ) , coldens_h_1  )
    ! coldens_h2  = SIGN( EXP(coldens_h2 ) , coldens_h2_1 )
    ! coldens_co  = SIGN( EXP(coldens_co ) , coldens_co_1 )
    ! Col_Dens_nH = SIGN( EXP(Col_Dens_nH) , Col_Dens_nH_1)
    ! IF (F_AV /= 0) THEN
    !    DO i = 1, nangle
    !       irf(i,1:nwlg)    = EXP(irf(i,1:nwlg)     )
    !    ENDDO
    ! ENDIF
    ! DO i = 1, NH2_lev
    !    H2_lev(i)%cd_l      = EXP(H2_lev(i)%cd_l    )
    ! ENDDO
    ! DO i = 1, Nspec
    !    speci(i)%density    = EXP(speci(i)%density  )
    !    speci(i)%Col_dens   = EXP(speci(i)%Col_dens )
    ! ENDDO
    ! DO i = 1, NH2_lev
    !    H2_lev(i)%density   = EXP(H2_lev(i)%density )
    !    H2_lev(i)%Col_dens  = EXP(H2_lev(i)%Col_dens)
    ! ENDDO

    DensityN = DensityN_2 * rhoN / rhoN_2
    compr_n  = compr_n_2  * rhoN / rhoN_2
    DO i = 1, Nspec
       speci(i)%density = speci_2(i)%density * rhoN / rhoN_2
    ENDDO
    DO i = 1, NH2_lev
       H2_lev(i)%density = H2_lev_2(i)%density * rhoN / rhoN_2
    ENDDO
    DO i = 1, NCO_lev
       CO_lev(i)%density = CO_lev_2(i)%density * rhoN / rhoN_2
    ENDDO



    dist_step  = distance - distance_1

    v_variab(iv_Vn        ) = Vn
    v_variab(iv_Vi        ) = Vi
    v_variab(iv_RhoN      ) = RhoN
    v_variab(iv_RhoI      ) = RhoI
    v_variab(iv_RhoA      ) = RhoA
    v_variab(iv_RhoNEG    ) = RhoNEG
    v_variab(iv_Tn        ) = Tn
    v_variab(iv_Ti        ) = Ti
    v_variab(iv_Te        ) = Te
    v_variab(iv_DensityN  ) = DensityN
    v_variab(iv_DensityI  ) = DensityI
    v_variab(iv_DensityA  ) = DensityA
    v_variab(iv_DensityNeg) = DensityNeg
    v_variab(iv_gv        ) = grad_V
    v_variab(iv_nh        ) = coldens_h
    v_variab(iv_nh2       ) = coldens_h2
    v_variab(iv_nco       ) = coldens_co
    v_variab(iv_vCH       ) = v_CH
    v_variab(iv_vS        ) = v_S
    v_variab(iv_vSH       ) = v_SH
    v_variab(iv_compr_n   ) = compr_n
    v_variab(iv_compr_i   ) = compr_i

    v_variab(bv_speci:ev_speci) = speci(1:Nspec)%density

    v_variab(bv_H2_lev:ev_H2_lev) = H2_lev(1:NH2_lev_var)%density

    WHERE (v_variab > 0)
       v_lvariab = LOG(v_variab)
    ELSEWHERE
       v_lvariab = minus_infinity
    END WHERE

    ! -------------------------------------------
    ! Deallocate all tables
    ! -------------------------------------------
    DEALLOCATE (H2_lev_1)
    DEALLOCATE (H2_lev_2)
    DEALLOCATE (CO_lev_1)
    DEALLOCATE (CO_lev_2)
    DEALLOCATE(speci_1)
    DEALLOCATE(speci_2)

  END SUBROUTINE INTERP_JUMP_CONDITIONS


  SUBROUTINE COMPUTE_ORTHO_TRAJ(imax,eps,coef1,coef2)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     compute the vector (DVn,DVi) over a segment of eps*distance
    !     and compute a vector orthogonal to this one, with a norm = Vi
    !     by convention, we impose the orthogonal vector to have a positive coef2
    ! subroutine/function needed :
    ! input variables :
    !     imax = the index of current distance
    !     eps  = determines the length of the vector
    ! output variables :
    !     coef1, coef2 = orthogonal vector coefficient
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_PROFIL_TABLES

    IMPLICIT none

    INTEGER,       INTENT(in)                :: imax
    REAL(KIND=DP), INTENT(in)                :: eps
    REAL(KIND=DP), INTENT(out)               :: coef1
    REAL(KIND=DP), INTENT(out)               :: coef2
    REAL(KIND=DP)                            :: coefa
    REAL(KIND=DP)                            :: coefb
    REAL(KIND=DP)                            :: dmin
    REAL(KIND=DP)                            :: dmax
    ! REAL(KIND=DP)                            :: vnmin, vimin
    ! REAL(KIND=DP)                            :: vnmax, vimax
    INTEGER                                  :: imin

    ! -------------------------------------------
    ! Find the range of indices over
    ! wich the vector is computed
    ! -------------------------------------------
    dmax = traj_main(imax)%distance
    imin = imax
    dmin = traj_main(imin)%distance
    DO WHILE( ABS(dmax-dmin) / dmax < eps )
       imin = imin - 1
       dmin = traj_main(imin)%distance
    ENDDO

    ! -------------------------------------------
    ! Find the range of indices over
    ! wich the vector is computed
    ! -------------------------------------------
    ! vnmax = traj_main(imax)%Vn
    ! vimax = traj_main(imax)%Vi
    ! imin = imax
    ! vnmin = traj_main(imin)%Vn
    ! vimin = traj_main(imin)%Vi
    ! DO WHILE( ABS(vnmax-vnmin) / vnmin < eps .AND. ABS(vimax-vimin) / vimin < eps )
    !    imin = imin - 1
    !    vnmin = traj_main(imin)%Vn
    !    vimin = traj_main(imin)%Vi
    ! ENDDO

    coefa = traj_main(imax)%Vn - traj_main(imin)%Vn
    coefb = traj_main(imax)%Vi - traj_main(imin)%Vi

    IF      ( ABS(coefb) > 1e4_dp * ABS(coefa) ) THEN
       IF ( coefa / coefb < 0.0_dp ) THEN
          coef1 =   traj_main(imax)%Vi / SQRT( (coefa*coefa) / (coefb*coefb) + 1.0_dp )
          coef2 = - coef1 * coefa / coefb 
       ELSE
          coef1 = - traj_main(imax)%Vi / SQRT( (coefa*coefa) / (coefb*coefb) + 1.0_dp )
          coef2 =   coef1 * coefa / coefb 
       ENDIF
    ELSE IF ( ABS(coefa) > 1e4_dp * ABS(coefb) ) THEN
       coef2 =   traj_main(imax)%Vi / SQRT( (coefb*coefb) / (coefa*coefa) + 1.0_dp )
       coef1 = - coef2 * coefb / coefa
    ELSE
       coef2 =   traj_main(imax)%Vi / SQRT( (coefb*coefb) / (coefa*coefa) + 1.0_dp )
       coef1 = - coef2 * coefb / coefa
    ENDIF

  END SUBROUTINE COMPUTE_ORTHO_TRAJ


  SUBROUTINE ROTATE_ORTHO_TRAJ(coef1,coef2,angle)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     rotate the orthogonal vector by an input angle
    ! subroutine/function needed :
    ! input variables :
    !     coef1, coef2 = orthogonal vector coefficient
    !     angle = rotation angle
    ! output variables :
    !     coef1, coef2 = rotated orthogonal vector coefficient
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_PROFIL_TABLES

    IMPLICIT none

    REAL(KIND=DP), INTENT(inout)             :: coef1
    REAL(KIND=DP), INTENT(inout)             :: coef2
    REAL(KIND=DP), INTENT(in)                :: angle
    REAL(KIND=DP)                            :: a
    REAL(KIND=DP)                            :: b

    ! -------------------------------------------
    ! Find the range of indices over
    ! wich the vector is computed
    ! -------------------------------------------
    a = cos(angle) * coef1 - sin(angle) * coef2
    b = sin(angle) * coef1 + cos(angle) * coef2

    coef1 = a
    coef2 = b

  END SUBROUTINE ROTATE_ORTHO_TRAJ


  SUBROUTINE AVERAGE_TRAJEC(traj1, traj2, trajf, N, Ntraj2, ka, kb, adjst)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     average two trajectories from indices ka to kb (in first trajectory)
    ! subroutine/function needed :
    ! input variables :
    !     traj1 and traj2 = the two trajectories to average
    !     ka and kb = range of indices in the first trajectory
    !     N = size og the traj tables
    !     Ntraj2 = max iteration in the second trajectory
    !     adjst  = exception index - do not interpolate traj2 on this point
    ! output variables :
    ! results :
    !     trajf from ka to kb
    !---------------------------------------------------------------------------
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_CHEM_REACT
    USE MODULE_H2
    USE MODULE_CO
    USE MODULE_LINE_EXCIT
    USE MODULE_PROFIL_TABLES

    IMPLICIT NONE

    INTEGER,                              INTENT(in)    :: N
    INTEGER,                              INTENT(in)    :: Ntraj2
    INTEGER,                              INTENT(in)    :: ka
    INTEGER,                              INTENT(in)    :: kb
    INTEGER,                              INTENT(in)    :: adjst
    TYPE (TYPE_TRAJECTORY), DIMENSION(N), INTENT(in)    :: traj1
    TYPE (TYPE_TRAJECTORY), DIMENSION(N), INTENT(in)    :: traj2
    TYPE (TYPE_TRAJECTORY), DIMENSION(N), INTENT(inout) :: trajf
    REAL(KIND=DP)                                       :: d1
    REAL(KIND=DP)                                       :: d2
    REAL(KIND=DP)                                       :: dum
    INTEGER                                             :: k1
    INTEGER                                             :: k2
    INTEGER                                             :: i

    DO k1 = ka, kb
       d1 = traj1(k1)%distance
       CALL FIND_DISTANCE(d1, N, traj2,  Ntraj2, k2)
       ! k2 = MAX(k2,ka)
       k2 = MAX(k2, adjst+1)
       d2 = traj2(k2)%distance
       IF ( traj2(k2+1)%distance - d2 /= 0 .AND. k2 /= Ntraj2 ) THEN
          dum = (d1 - d2) / (traj2(k2+1)%distance - d2)
       ELSE
          dum = 0.0_dp
       ENDIF
       !---------------------------------------------------------------------------
       ! time, distance
       !---------------------------------------------------------------------------
       trajf(k1)%slab         = traj1(k1)%slab
       trajf(k1)%distance     = traj1(k1)%distance
       trajf(k1)%dist_step    = traj1(k1)%dist_step
       trajf(k1)%Av           = AVER(traj1(k1)%Av          ,traj2(k2)%Av            ,traj2(k2+1)%Av           ,dum)
       trajf(k1)%timeN        = AVER(traj1(k1)%timeN       ,traj2(k2)%timeN         ,traj2(k2+1)%timeN        ,dum)
       trajf(k1)%timeI        = AVER(traj1(k1)%timeI       ,traj2(k2)%timeI         ,traj2(k2+1)%timeI        ,dum)

       !---------------------------------------------------------------------------
       ! densities
       !---------------------------------------------------------------------------
       trajf(k1)%nH           = AVER(traj1(k1)%nH          ,traj2(k2)%nH            ,traj2(k2+1)%nH           ,dum)
       trajf(k1)%rhoN         = AVER(traj1(k1)%rhoN        ,traj2(k2)%rhoN          ,traj2(k2+1)%rhoN         ,dum)
       trajf(k1)%rhoI         = AVER(traj1(k1)%rhoI        ,traj2(k2)%rhoI          ,traj2(k2+1)%rhoI         ,dum)
       trajf(k1)%rhoA         = AVER(traj1(k1)%rhoA        ,traj2(k2)%rhoA          ,traj2(k2+1)%rhoA         ,dum)
       trajf(k1)%rhoNeg       = AVER(traj1(k1)%rhoNeg      ,traj2(k2)%rhoNeg        ,traj2(k2+1)%rhoNeg       ,dum)
       trajf(k1)%DensityN     = AVER(traj1(k1)%DensityN    ,traj2(k2)%DensityN      ,traj2(k2+1)%DensityN     ,dum)
       trajf(k1)%DensityI     = AVER(traj1(k1)%DensityI    ,traj2(k2)%DensityI      ,traj2(k2+1)%DensityI     ,dum)
       trajf(k1)%DensityA     = AVER(traj1(k1)%DensityA    ,traj2(k2)%DensityA      ,traj2(k2+1)%DensityA     ,dum)
       trajf(k1)%DensityNeg   = AVER(traj1(k1)%DensityNeg  ,traj2(k2)%DensityNeg    ,traj2(k2+1)%DensityNeg   ,dum)
       trajf(k1)%muN          = AVER(traj1(k1)%muN         ,traj2(k2)%muN           ,traj2(k2+1)%muN          ,dum)
       trajf(k1)%muI          = AVER(traj1(k1)%muI         ,traj2(k2)%muI           ,traj2(k2+1)%muI          ,dum)
       trajf(k1)%muA          = AVER(traj1(k1)%muA         ,traj2(k2)%muA           ,traj2(k2+1)%muA          ,dum)
       trajf(k1)%muNeg        = AVER(traj1(k1)%muNeg       ,traj2(k2)%muNeg         ,traj2(k2+1)%muNeg        ,dum)
       trajf(k1)%iondeg       = AVER(traj1(k1)%iondeg      ,traj2(k2)%iondeg        ,traj2(k2+1)%iondeg       ,dum)
       trajf(k1)%ionfrac      = AVER(traj1(k1)%ionfrac     ,traj2(k2)%ionfrac       ,traj2(k2+1)%ionfrac      ,dum)

       !---------------------------------------------------------------------------
       ! velocity
       !---------------------------------------------------------------------------
       trajf(k1)%Vsound       = AVER(traj1(k1)%Vsound      ,traj2(k2)%Vsound        ,traj2(k2+1)%Vsound       ,dum)
       trajf(k1)%Vmagnet      = AVER(traj1(k1)%Vmagnet     ,traj2(k2)%Vmagnet       ,traj2(k2+1)%Vmagnet      ,dum)
       trajf(k1)%Valfven      = AVER(traj1(k1)%Valfven     ,traj2(k2)%Valfven       ,traj2(k2+1)%Valfven      ,dum)
       trajf(k1)%Vn           = AVER(traj1(k1)%Vn          ,traj2(k2)%Vn            ,traj2(k2+1)%Vn           ,dum)
       trajf(k1)%Vi           = AVER(traj1(k1)%Vi          ,traj2(k2)%Vi            ,traj2(k2+1)%Vi           ,dum)
       trajf(k1)%dVn          = AVER(traj1(k1)%dVn         ,traj2(k2)%dVn           ,traj2(k2+1)%dVn          ,dum)
       trajf(k1)%dVi          = AVER(traj1(k1)%dVi         ,traj2(k2)%dVi           ,traj2(k2+1)%dVi          ,dum)
       trajf(k1)%Vgrad        = AVER(traj1(k1)%Vgrad       ,traj2(k2)%Vgrad         ,traj2(k2+1)%Vgrad        ,dum)
       trajf(k1)%grad_V       = AVER(traj1(k1)%grad_V      ,traj2(k2)%grad_V        ,traj2(k2+1)%grad_V       ,dum)
       trajf(k1)%v_CH         = AVER(traj1(k1)%v_CH        ,traj2(k2)%v_CH          ,traj2(k2+1)%v_CH         ,dum)
       trajf(k1)%v_S          = AVER(traj1(k1)%v_S         ,traj2(k2)%v_S           ,traj2(k2+1)%v_S          ,dum)
       trajf(k1)%v_SH         = AVER(traj1(k1)%v_SH        ,traj2(k2)%v_SH          ,traj2(k2+1)%v_SH         ,dum)

       !---------------------------------------------------------------------------
       ! temperature
       !---------------------------------------------------------------------------
       trajf(k1)%Tn           = AVER(traj1(k1)%Tn           ,traj2(k2)%Tn           ,traj2(k2+1)%Tn           ,dum)
       trajf(k1)%Ti           = AVER(traj1(k1)%Ti           ,traj2(k2)%Ti           ,traj2(k2+1)%Ti           ,dum)
       trajf(k1)%Te           = AVER(traj1(k1)%Te           ,traj2(k2)%Te           ,traj2(k2+1)%Te           ,dum)
       trajf(k1)%T_gr         = AVER(traj1(k1)%T_gr         ,traj2(k2)%T_gr         ,traj2(k2+1)%T_gr         ,dum)
       trajf(k1)%Teff_gr      = AVER(traj1(k1)%Teff_gr      ,traj2(k2)%Teff_gr      ,traj2(k2+1)%Teff_gr      ,dum)
       trajf(k1)%oph2         = AVER(traj1(k1)%oph2         ,traj2(k2)%oph2         ,traj2(k2+1)%oph2         ,dum)

       !---------------------------------------------------------------------------
       ! grains
       !---------------------------------------------------------------------------
       trajf(k1)%nlay_gr      = AVER(traj1(k1)%nlay_gr      ,traj2(k2)%nlay_gr      ,traj2(k2+1)%nlay_gr      ,dum)
       trajf(k1)%n_gr         = AVER(traj1(k1)%n_gr         ,traj2(k2)%n_gr         ,traj2(k2+1)%n_gr         ,dum)
       trajf(k1)%mu_gr        = AVER(traj1(k1)%mu_gr        ,traj2(k2)%mu_gr        ,traj2(k2+1)%mu_gr        ,dum)
       trajf(k1)%much_gr      = AVER(traj1(k1)%much_gr      ,traj2(k2)%much_gr      ,traj2(k2+1)%much_gr      ,dum)
       trajf(k1)%r_gr         = AVER(traj1(k1)%r_gr         ,traj2(k2)%r_gr         ,traj2(k2+1)%r_gr         ,dum)
       trajf(k1)%m_gr         = AVER(traj1(k1)%m_gr         ,traj2(k2)%m_gr         ,traj2(k2+1)%m_gr         ,dum)
       trajf(k1)%m_grc        = AVER(traj1(k1)%m_grc        ,traj2(k2)%m_grc        ,traj2(k2+1)%m_grc        ,dum)
       trajf(k1)%m_grm        = AVER(traj1(k1)%m_grm        ,traj2(k2)%m_grm        ,traj2(k2+1)%m_grm        ,dum)
       trajf(k1)%d_ero        = AVER(traj1(k1)%d_ero        ,traj2(k2)%d_ero        ,traj2(k2+1)%d_ero        ,dum)
       trajf(k1)%d_ads        = AVER(traj1(k1)%d_ads        ,traj2(k2)%d_ads        ,traj2(k2+1)%d_ads        ,dum)
       trajf(k1)%rsq_grc      = AVER(traj1(k1)%rsq_grc      ,traj2(k2)%rsq_grc      ,traj2(k2+1)%rsq_grc      ,dum)
       trajf(k1)%rsq_grm      = AVER(traj1(k1)%rsq_grm      ,traj2(k2)%rsq_grm      ,traj2(k2+1)%rsq_grm      ,dum)
       trajf(k1)%compr_n      = AVER(traj1(k1)%compr_n      ,traj2(k2)%compr_n      ,traj2(k2+1)%compr_n      ,dum)
       trajf(k1)%compr_i      = AVER(traj1(k1)%compr_i      ,traj2(k2)%compr_i      ,traj2(k2+1)%compr_i      ,dum)

       !---------------------------------------------------------------------------
       ! energetics
       !---------------------------------------------------------------------------
       trajf(k1)%mom_flux_tot = AVER(traj1(k1)%mom_flux_tot ,traj2(k2)%mom_flux_tot ,traj2(k2+1)%mom_flux_tot ,dum)
       trajf(k1)%mom_flux_kin = AVER(traj1(k1)%mom_flux_kin ,traj2(k2)%mom_flux_kin ,traj2(k2+1)%mom_flux_kin ,dum)
       trajf(k1)%mom_flux_the = AVER(traj1(k1)%mom_flux_the ,traj2(k2)%mom_flux_the ,traj2(k2+1)%mom_flux_the ,dum)
       trajf(k1)%mom_flux_mag = AVER(traj1(k1)%mom_flux_mag ,traj2(k2)%mom_flux_mag ,traj2(k2+1)%mom_flux_mag ,dum)
       trajf(k1)%mom_flux_vis = AVER(traj1(k1)%mom_flux_vis ,traj2(k2)%mom_flux_vis ,traj2(k2+1)%mom_flux_vis ,dum)
       trajf(k1)%nrj_flux_tot = AVER(traj1(k1)%nrj_flux_tot ,traj2(k2)%nrj_flux_tot ,traj2(k2+1)%nrj_flux_tot ,dum)
       trajf(k1)%nrj_flux_kin = AVER(traj1(k1)%nrj_flux_kin ,traj2(k2)%nrj_flux_kin ,traj2(k2+1)%nrj_flux_kin ,dum)
       trajf(k1)%nrj_flux_the = AVER(traj1(k1)%nrj_flux_the ,traj2(k2)%nrj_flux_the ,traj2(k2+1)%nrj_flux_the ,dum)
       trajf(k1)%nrj_flux_mag = AVER(traj1(k1)%nrj_flux_mag ,traj2(k2)%nrj_flux_mag ,traj2(k2+1)%nrj_flux_mag ,dum)
       trajf(k1)%nrj_flux_vis = AVER(traj1(k1)%nrj_flux_vis ,traj2(k2)%nrj_flux_vis ,traj2(k2+1)%nrj_flux_vis ,dum)
       trajf(k1)%nrj_flux_int = AVER(traj1(k1)%nrj_flux_int ,traj2(k2)%nrj_flux_int ,traj2(k2+1)%nrj_flux_int ,dum)
       trajf(k1)%nrj_flux_src = AVER(traj1(k1)%nrj_flux_src ,traj2(k2)%nrj_flux_src ,traj2(k2+1)%nrj_flux_src ,dum)
       trajf(k1)%nrj_flux_pho = AVER(traj1(k1)%nrj_flux_pho ,traj2(k2)%nrj_flux_pho ,traj2(k2+1)%nrj_flux_pho ,dum)
       trajf(k1)%mag_flux_corr= AVER(traj1(k1)%mag_flux_corr,traj2(k2)%mag_flux_corr,traj2(k2+1)%mag_flux_corr,dum)

       !---------------------------------------------------------------------------
       ! heating / cooling
       !---------------------------------------------------------------------------
       trajf(k1)%line_rad_tot      = AVER(traj1(k1)%line_rad_tot     ,traj2(k2)%line_rad_tot     ,traj2(k2+1)%line_rad_tot     ,dum)
       trajf(k1)%line_rad_n        = AVER(traj1(k1)%line_rad_n       ,traj2(k2)%line_rad_n       ,traj2(k2+1)%line_rad_n       ,dum)
       trajf(k1)%line_rad_i        = AVER(traj1(k1)%line_rad_i       ,traj2(k2)%line_rad_i       ,traj2(k2+1)%line_rad_i       ,dum)
       trajf(k1)%line_rad_e        = AVER(traj1(k1)%line_rad_e       ,traj2(k2)%line_rad_e       ,traj2(k2+1)%line_rad_e       ,dum)
       trajf(k1)%line_rad_n_molec  = AVER(traj1(k1)%line_rad_n_molec ,traj2(k2)%line_rad_n_molec ,traj2(k2+1)%line_rad_n_molec ,dum)
       trajf(k1)%line_rad_n_H2     = AVER(traj1(k1)%line_rad_n_H2    ,traj2(k2)%line_rad_n_H2    ,traj2(k2+1)%line_rad_n_H2    ,dum)
       trajf(k1)%line_rad_n_13CO   = AVER(traj1(k1)%line_rad_n_13CO  ,traj2(k2)%line_rad_n_13CO  ,traj2(k2+1)%line_rad_n_13CO  ,dum)
       trajf(k1)%line_rad_n_OH     = AVER(traj1(k1)%line_rad_n_OH    ,traj2(k2)%line_rad_n_OH    ,traj2(k2+1)%line_rad_n_OH    ,dum)
       trajf(k1)%line_rad_n_NH3    = AVER(traj1(k1)%line_rad_n_NH3   ,traj2(k2)%line_rad_n_NH3   ,traj2(k2+1)%line_rad_n_NH3   ,dum)
       trajf(k1)%line_rad_n_rCO    = AVER(traj1(k1)%line_rad_n_rCO   ,traj2(k2)%line_rad_n_rCO   ,traj2(k2+1)%line_rad_n_rCO   ,dum)
       trajf(k1)%line_rad_n_vCO    = AVER(traj1(k1)%line_rad_n_vCO   ,traj2(k2)%line_rad_n_vCO   ,traj2(k2+1)%line_rad_n_vCO   ,dum)
       trajf(k1)%line_rad_n_CO     = AVER(traj1(k1)%line_rad_n_CO    ,traj2(k2)%line_rad_n_CO    ,traj2(k2+1)%line_rad_n_CO    ,dum)
       trajf(k1)%line_rad_n_roH2O  = AVER(traj1(k1)%line_rad_n_roH2O ,traj2(k2)%line_rad_n_roH2O ,traj2(k2+1)%line_rad_n_roH2O ,dum)
       trajf(k1)%line_rad_n_rpH2O  = AVER(traj1(k1)%line_rad_n_rpH2O ,traj2(k2)%line_rad_n_rpH2O ,traj2(k2+1)%line_rad_n_rpH2O ,dum)
       trajf(k1)%line_rad_n_vH2O   = AVER(traj1(k1)%line_rad_n_vH2O  ,traj2(k2)%line_rad_n_vH2O  ,traj2(k2+1)%line_rad_n_vH2O  ,dum)
       trajf(k1)%line_rad_n_H2O    = AVER(traj1(k1)%line_rad_n_H2O   ,traj2(k2)%line_rad_n_H2O   ,traj2(k2+1)%line_rad_n_H2O   ,dum)
       trajf(k1)%line_rad_i_molec  = AVER(traj1(k1)%line_rad_i_molec ,traj2(k2)%line_rad_i_molec ,traj2(k2+1)%line_rad_i_molec ,dum)
       trajf(k1)%line_rad_i_H2     = AVER(traj1(k1)%line_rad_i_H2    ,traj2(k2)%line_rad_i_H2    ,traj2(k2+1)%line_rad_i_H2    ,dum)
       trajf(k1)%line_rad_e_molec  = AVER(traj1(k1)%line_rad_e_molec ,traj2(k2)%line_rad_e_molec ,traj2(k2+1)%line_rad_e_molec ,dum)
       trajf(k1)%line_rad_e_H2     = AVER(traj1(k1)%line_rad_e_H2    ,traj2(k2)%line_rad_e_H2    ,traj2(k2+1)%line_rad_e_H2    ,dum)
       trajf(k1)%line_rad_e_CO     = AVER(traj1(k1)%line_rad_e_CO    ,traj2(k2)%line_rad_e_CO    ,traj2(k2+1)%line_rad_e_CO    ,dum)
       trajf(k1)%line_rad_n_atoms  = AVER(traj1(k1)%line_rad_n_atoms ,traj2(k2)%line_rad_n_atoms ,traj2(k2+1)%line_rad_n_atoms ,dum)
       trajf(k1)%line_rad_n_Hat    = AVER(traj1(k1)%line_rad_n_Hat   ,traj2(k2)%line_rad_n_Hat   ,traj2(k2+1)%line_rad_n_Hat   ,dum)
       trajf(k1)%line_rad_n_Oat    = AVER(traj1(k1)%line_rad_n_Oat   ,traj2(k2)%line_rad_n_Oat   ,traj2(k2+1)%line_rad_n_Oat   ,dum)
       trajf(k1)%line_rad_n_Op     = AVER(traj1(k1)%line_rad_n_Op    ,traj2(k2)%line_rad_n_Op    ,traj2(k2+1)%line_rad_n_Op    ,dum)
       trajf(k1)%line_rad_n_Cat    = AVER(traj1(k1)%line_rad_n_Cat   ,traj2(k2)%line_rad_n_Cat   ,traj2(k2+1)%line_rad_n_Cat   ,dum)
       trajf(k1)%line_rad_n_Cp     = AVER(traj1(k1)%line_rad_n_Cp    ,traj2(k2)%line_rad_n_Cp    ,traj2(k2+1)%line_rad_n_Cp    ,dum)
       trajf(k1)%line_rad_n_Nat    = AVER(traj1(k1)%line_rad_n_Nat   ,traj2(k2)%line_rad_n_Nat   ,traj2(k2+1)%line_rad_n_Nat   ,dum)
       trajf(k1)%line_rad_n_Np     = AVER(traj1(k1)%line_rad_n_Np    ,traj2(k2)%line_rad_n_Np    ,traj2(k2+1)%line_rad_n_Np    ,dum)
       trajf(k1)%line_rad_n_Sat    = AVER(traj1(k1)%line_rad_n_Sat   ,traj2(k2)%line_rad_n_Sat   ,traj2(k2+1)%line_rad_n_Sat   ,dum)
       trajf(k1)%line_rad_n_Sp     = AVER(traj1(k1)%line_rad_n_Sp    ,traj2(k2)%line_rad_n_Sp    ,traj2(k2+1)%line_rad_n_Sp    ,dum)
       trajf(k1)%line_rad_n_Siat   = AVER(traj1(k1)%line_rad_n_Siat  ,traj2(k2)%line_rad_n_Siat  ,traj2(k2+1)%line_rad_n_Siat  ,dum)
       trajf(k1)%line_rad_n_Sip    = AVER(traj1(k1)%line_rad_n_Sip   ,traj2(k2)%line_rad_n_Sip   ,traj2(k2+1)%line_rad_n_Sip   ,dum)
       trajf(k1)%line_rad_i_atoms  = AVER(traj1(k1)%line_rad_i_atoms ,traj2(k2)%line_rad_i_atoms ,traj2(k2+1)%line_rad_i_atoms ,dum)
       trajf(k1)%line_rad_i_Hat    = AVER(traj1(k1)%line_rad_i_Hat   ,traj2(k2)%line_rad_i_Hat   ,traj2(k2+1)%line_rad_i_Hat   ,dum)
       trajf(k1)%line_rad_i_Oat    = AVER(traj1(k1)%line_rad_i_Oat   ,traj2(k2)%line_rad_i_Oat   ,traj2(k2+1)%line_rad_i_Oat   ,dum)
       trajf(k1)%line_rad_i_Op     = AVER(traj1(k1)%line_rad_i_Op    ,traj2(k2)%line_rad_i_Op    ,traj2(k2+1)%line_rad_i_Op    ,dum)
       trajf(k1)%line_rad_i_Cat    = AVER(traj1(k1)%line_rad_i_Cat   ,traj2(k2)%line_rad_i_Cat   ,traj2(k2+1)%line_rad_i_Cat   ,dum)
       trajf(k1)%line_rad_i_Cp     = AVER(traj1(k1)%line_rad_i_Cp    ,traj2(k2)%line_rad_i_Cp    ,traj2(k2+1)%line_rad_i_Cp    ,dum)
       trajf(k1)%line_rad_i_Nat    = AVER(traj1(k1)%line_rad_i_Nat   ,traj2(k2)%line_rad_i_Nat   ,traj2(k2+1)%line_rad_i_Nat   ,dum)
       trajf(k1)%line_rad_i_Np     = AVER(traj1(k1)%line_rad_i_Np    ,traj2(k2)%line_rad_i_Np    ,traj2(k2+1)%line_rad_i_Np    ,dum)
       trajf(k1)%line_rad_i_Sat    = AVER(traj1(k1)%line_rad_i_Sat   ,traj2(k2)%line_rad_i_Sat   ,traj2(k2+1)%line_rad_i_Sat   ,dum)
       trajf(k1)%line_rad_i_Sp     = AVER(traj1(k1)%line_rad_i_Sp    ,traj2(k2)%line_rad_i_Sp    ,traj2(k2+1)%line_rad_i_Sp    ,dum)
       trajf(k1)%line_rad_i_Siat   = AVER(traj1(k1)%line_rad_i_Siat  ,traj2(k2)%line_rad_i_Siat  ,traj2(k2+1)%line_rad_i_Siat  ,dum)
       trajf(k1)%line_rad_i_Sip    = AVER(traj1(k1)%line_rad_i_Sip   ,traj2(k2)%line_rad_i_Sip   ,traj2(k2+1)%line_rad_i_Sip   ,dum)
       trajf(k1)%line_rad_e_atoms  = AVER(traj1(k1)%line_rad_e_atoms ,traj2(k2)%line_rad_e_atoms ,traj2(k2+1)%line_rad_e_atoms ,dum)
       trajf(k1)%line_rad_e_Hat    = AVER(traj1(k1)%line_rad_e_Hat   ,traj2(k2)%line_rad_e_Hat   ,traj2(k2+1)%line_rad_e_Hat   ,dum)
       trajf(k1)%line_rad_e_Oat    = AVER(traj1(k1)%line_rad_e_Oat   ,traj2(k2)%line_rad_e_Oat   ,traj2(k2+1)%line_rad_e_Oat   ,dum)
       trajf(k1)%line_rad_e_Op     = AVER(traj1(k1)%line_rad_e_Op    ,traj2(k2)%line_rad_e_Op    ,traj2(k2+1)%line_rad_e_Op    ,dum)
       trajf(k1)%line_rad_e_Cat    = AVER(traj1(k1)%line_rad_e_Cat   ,traj2(k2)%line_rad_e_Cat   ,traj2(k2+1)%line_rad_e_Cat   ,dum)
       trajf(k1)%line_rad_e_Cp     = AVER(traj1(k1)%line_rad_e_Cp    ,traj2(k2)%line_rad_e_Cp    ,traj2(k2+1)%line_rad_e_Cp    ,dum)
       trajf(k1)%line_rad_e_Nat    = AVER(traj1(k1)%line_rad_e_Nat   ,traj2(k2)%line_rad_e_Nat   ,traj2(k2+1)%line_rad_e_Nat   ,dum)
       trajf(k1)%line_rad_e_Np     = AVER(traj1(k1)%line_rad_e_Np    ,traj2(k2)%line_rad_e_Np    ,traj2(k2+1)%line_rad_e_Np    ,dum)
       trajf(k1)%line_rad_e_Sat    = AVER(traj1(k1)%line_rad_e_Sat   ,traj2(k2)%line_rad_e_Sat   ,traj2(k2+1)%line_rad_e_Sat   ,dum)
       trajf(k1)%line_rad_e_Sp     = AVER(traj1(k1)%line_rad_e_Sp    ,traj2(k2)%line_rad_e_Sp    ,traj2(k2+1)%line_rad_e_Sp    ,dum)
       trajf(k1)%line_rad_e_Siat   = AVER(traj1(k1)%line_rad_e_Siat  ,traj2(k2)%line_rad_e_Siat  ,traj2(k2+1)%line_rad_e_Siat  ,dum)
       trajf(k1)%line_rad_e_Sip    = AVER(traj1(k1)%line_rad_e_Sip   ,traj2(k2)%line_rad_e_Sip   ,traj2(k2+1)%line_rad_e_Sip   ,dum)
       trajf(k1)%chem_DE_tot       = AVER(traj1(k1)%chem_DE_tot      ,traj2(k2)%chem_DE_tot      ,traj2(k2+1)%chem_DE_tot      ,dum)
       trajf(k1)%chem_DE_n         = AVER(traj1(k1)%chem_DE_n        ,traj2(k2)%chem_DE_n        ,traj2(k2+1)%chem_DE_n        ,dum)
       trajf(k1)%chem_DE_i         = AVER(traj1(k1)%chem_DE_i        ,traj2(k2)%chem_DE_i        ,traj2(k2+1)%chem_DE_i        ,dum)
       trajf(k1)%chem_DE_e         = AVER(traj1(k1)%chem_DE_e        ,traj2(k2)%chem_DE_e        ,traj2(k2+1)%chem_DE_e        ,dum)
       trajf(k1)%chem_DE_n_diss_H2 = AVER(traj1(k1)%chem_DE_n_diss_H2,traj2(k2)%chem_DE_n_diss_H2,traj2(k2+1)%chem_DE_n_diss_H2,dum)
       trajf(k1)%chem_DE_i_diss_H2 = AVER(traj1(k1)%chem_DE_i_diss_H2,traj2(k2)%chem_DE_i_diss_H2,traj2(k2+1)%chem_DE_i_diss_H2,dum)
       trajf(k1)%chem_DE_e_diss_H2 = AVER(traj1(k1)%chem_DE_e_diss_H2,traj2(k2)%chem_DE_e_diss_H2,traj2(k2+1)%chem_DE_e_diss_H2,dum)
       trajf(k1)%chem_DE_n_phgas   = AVER(traj1(k1)%chem_DE_n_phgas  ,traj2(k2)%chem_DE_n_phgas  ,traj2(k2+1)%chem_DE_n_phgas  ,dum)
       trajf(k1)%chem_DE_i_phgas   = AVER(traj1(k1)%chem_DE_i_phgas  ,traj2(k2)%chem_DE_i_phgas  ,traj2(k2+1)%chem_DE_i_phgas  ,dum)
       trajf(k1)%chem_DE_e_phgas   = AVER(traj1(k1)%chem_DE_e_phgas  ,traj2(k2)%chem_DE_e_phgas  ,traj2(k2+1)%chem_DE_e_phgas  ,dum)
       trajf(k1)%chem_DE_n_phgrn   = AVER(traj1(k1)%chem_DE_n_phgrn  ,traj2(k2)%chem_DE_n_phgrn  ,traj2(k2+1)%chem_DE_n_phgrn  ,dum)
       trajf(k1)%chem_DE_i_phgrn   = AVER(traj1(k1)%chem_DE_i_phgrn  ,traj2(k2)%chem_DE_i_phgrn  ,traj2(k2+1)%chem_DE_i_phgrn  ,dum)
       trajf(k1)%chem_DE_e_phgrn   = AVER(traj1(k1)%chem_DE_e_phgrn  ,traj2(k2)%chem_DE_e_phgrn  ,traj2(k2+1)%chem_DE_e_phgrn  ,dum)
       trajf(k1)%chem_DE_n_cosmic  = AVER(traj1(k1)%chem_DE_n_cosmic ,traj2(k2)%chem_DE_n_cosmic ,traj2(k2+1)%chem_DE_n_cosmic ,dum)
       trajf(k1)%chem_DE_i_cosmic  = AVER(traj1(k1)%chem_DE_i_cosmic ,traj2(k2)%chem_DE_i_cosmic ,traj2(k2+1)%chem_DE_i_cosmic ,dum)
       trajf(k1)%chem_DE_e_cosmic  = AVER(traj1(k1)%chem_DE_e_cosmic ,traj2(k2)%chem_DE_e_cosmic ,traj2(k2+1)%chem_DE_e_cosmic ,dum)
       trajf(k1)%chem_DE_n_other   = AVER(traj1(k1)%chem_DE_n_other  ,traj2(k2)%chem_DE_n_other  ,traj2(k2+1)%chem_DE_n_other  ,dum)
       trajf(k1)%chem_DE_i_other   = AVER(traj1(k1)%chem_DE_i_other  ,traj2(k2)%chem_DE_i_other  ,traj2(k2+1)%chem_DE_i_other  ,dum)
       trajf(k1)%chem_DE_e_other   = AVER(traj1(k1)%chem_DE_e_other  ,traj2(k2)%chem_DE_e_other  ,traj2(k2+1)%chem_DE_e_other  ,dum)
       trajf(k1)%elast_scat_tot    = AVER(traj1(k1)%elast_scat_tot   ,traj2(k2)%elast_scat_tot   ,traj2(k2+1)%elast_scat_tot   ,dum)
       trajf(k1)%elast_scat_n      = AVER(traj1(k1)%elast_scat_n     ,traj2(k2)%elast_scat_n     ,traj2(k2+1)%elast_scat_n     ,dum)
       trajf(k1)%elast_scat_i      = AVER(traj1(k1)%elast_scat_i     ,traj2(k2)%elast_scat_i     ,traj2(k2+1)%elast_scat_i     ,dum)
       trajf(k1)%elast_scat_e      = AVER(traj1(k1)%elast_scat_e     ,traj2(k2)%elast_scat_e     ,traj2(k2+1)%elast_scat_e     ,dum)
       trajf(k1)%exch_eint_tot     = AVER(traj1(k1)%exch_eint_tot    ,traj2(k2)%exch_eint_tot    ,traj2(k2+1)%exch_eint_tot    ,dum)
       trajf(k1)%exch_eint_n       = AVER(traj1(k1)%exch_eint_n      ,traj2(k2)%exch_eint_n      ,traj2(k2+1)%exch_eint_n      ,dum)
       trajf(k1)%exch_eint_i       = AVER(traj1(k1)%exch_eint_i      ,traj2(k2)%exch_eint_i      ,traj2(k2+1)%exch_eint_i      ,dum)
       trajf(k1)%exch_eint_e       = AVER(traj1(k1)%exch_eint_e      ,traj2(k2)%exch_eint_e      ,traj2(k2+1)%exch_eint_e      ,dum)
       trajf(k1)%therm_grain_tot   = AVER(traj1(k1)%therm_grain_tot  ,traj2(k2)%therm_grain_tot  ,traj2(k2+1)%therm_grain_tot  ,dum)
       trajf(k1)%mech_trsf_tot     = AVER(traj1(k1)%mech_trsf_tot    ,traj2(k2)%mech_trsf_tot    ,traj2(k2+1)%mech_trsf_tot    ,dum)
       trajf(k1)%mech_trsf_n       = AVER(traj1(k1)%mech_trsf_n      ,traj2(k2)%mech_trsf_n      ,traj2(k2+1)%mech_trsf_n      ,dum)
       trajf(k1)%mech_trsf_i       = AVER(traj1(k1)%mech_trsf_i      ,traj2(k2)%mech_trsf_i      ,traj2(k2+1)%mech_trsf_i      ,dum)
       trajf(k1)%mech_trsf_e       = AVER(traj1(k1)%mech_trsf_e      ,traj2(k2)%mech_trsf_e      ,traj2(k2+1)%mech_trsf_e      ,dum)
       trajf(k1)%mech_trsf_n_chem  = AVER(traj1(k1)%mech_trsf_n_chem ,traj2(k2)%mech_trsf_n_chem ,traj2(k2+1)%mech_trsf_n_chem ,dum)
       trajf(k1)%mech_trsf_i_chem  = AVER(traj1(k1)%mech_trsf_i_chem ,traj2(k2)%mech_trsf_i_chem ,traj2(k2+1)%mech_trsf_i_chem ,dum)
       trajf(k1)%mech_trsf_e_chem  = AVER(traj1(k1)%mech_trsf_e_chem ,traj2(k2)%mech_trsf_e_chem ,traj2(k2+1)%mech_trsf_e_chem ,dum)
       trajf(k1)%mech_trsf_n_compr = AVER(traj1(k1)%mech_trsf_n_compr,traj2(k2)%mech_trsf_n_compr,traj2(k2+1)%mech_trsf_n_compr,dum)
       trajf(k1)%mech_trsf_i_compr = AVER(traj1(k1)%mech_trsf_i_compr,traj2(k2)%mech_trsf_i_compr,traj2(k2+1)%mech_trsf_i_compr,dum)
       trajf(k1)%mech_trsf_e_compr = AVER(traj1(k1)%mech_trsf_e_compr,traj2(k2)%mech_trsf_e_compr,traj2(k2+1)%mech_trsf_e_compr,dum)
       trajf(k1)%mech_trsf_n_visc  = AVER(traj1(k1)%mech_trsf_n_visc ,traj2(k2)%mech_trsf_n_visc ,traj2(k2+1)%mech_trsf_n_visc ,dum)
       trajf(k1)%mech_trsf_i_visc  = AVER(traj1(k1)%mech_trsf_i_visc ,traj2(k2)%mech_trsf_i_visc ,traj2(k2+1)%mech_trsf_i_visc ,dum)
       trajf(k1)%mech_trsf_e_visc  = AVER(traj1(k1)%mech_trsf_e_visc ,traj2(k2)%mech_trsf_e_visc ,traj2(k2+1)%mech_trsf_e_visc ,dum)

       !---------------------------------------------------------------------------
       ! radiation
       !---------------------------------------------------------------------------
       trajf(k1)%fluph        = AVER(traj1(k1)%fluph        ,traj2(k2)%fluph        ,traj2(k2+1)%fluph        ,dum)
       trajf(k1)%coldens_h    = AVER(traj1(k1)%coldens_h    ,traj2(k2)%coldens_h    ,traj2(k2+1)%coldens_h    ,dum)
       trajf(k1)%coldens_h2   = AVER(traj1(k1)%coldens_h2   ,traj2(k2)%coldens_h2   ,traj2(k2+1)%coldens_h2   ,dum)
       trajf(k1)%coldens_co   = AVER(traj1(k1)%coldens_co   ,traj2(k2)%coldens_co   ,traj2(k2+1)%coldens_co   ,dum)

       !---------------------------------------------------------------------------
       ! chemical profiles
       !---------------------------------------------------------------------------
       DO i = 1, Nspec
          trajf(k1)%abon(i)    = AVER(traj1(k1)%abon(i)    ,traj2(k2)%abon(i)    ,traj2(k2+1)%abon(i)    ,dum)
          trajf(k1)%coldens(i) = AVER(traj1(k1)%coldens(i) ,traj2(k2)%coldens(i) ,traj2(k2+1)%coldens(i) ,dum)
       ENDDO
       ! DO i = 1, Nreact
       ! trajf(k1)%all_rate_chem(i) = AVER(traj1(k1)%all_rate_chem(i) ,traj2(k2)%all_rate_chem(i) ,traj2(k2+1)%all_rate_chem(i) ,dum)
       ! ENDDO

       !---------------------------------------------------------------------------
       ! excitation profiles
       !---------------------------------------------------------------------------
       DO i = 1, NH2_lev
          trajf(k1)%abon_h2lev(i)   = AVER(traj1(k1)%abon_h2lev(i)   ,traj2(k2)%abon_h2lev(i)   ,traj2(k2+1)%abon_h2lev(i)   ,dum)
          trajf(k1)%cold_h2lev(i)   = AVER(traj1(k1)%cold_h2lev(i)   ,traj2(k2)%cold_h2lev(i)   ,traj2(k2+1)%cold_h2lev(i)   ,dum)
          trajf(k1)%cdleft_h2lev(i) = AVER(traj1(k1)%cdleft_h2lev(i) ,traj2(k2)%cdleft_h2lev(i) ,traj2(k2+1)%cdleft_h2lev(i) ,dum)
       ENDDO
       DO i = 1, NH2_lines_out
          trajf(k1)%emis_h2lin(i)   = AVER(traj1(k1)%emis_h2lin(i)   ,traj2(k2)%emis_h2lin(i)   ,traj2(k2+1)%emis_h2lin(i)   ,dum)
          trajf(k1)%inta_h2lin(i)   = AVER(traj1(k1)%inta_h2lin(i)   ,traj2(k2)%inta_h2lin(i)   ,traj2(k2+1)%inta_h2lin(i)   ,dum)
       ENDDO
       DO i = 1, NCO_lev
          trajf(k1)%abon_colev(i)   = AVER(traj1(k1)%abon_colev(i)   ,traj2(k2)%abon_colev(i)   ,traj2(k2+1)%abon_colev(i)   ,dum)
          trajf(k1)%cold_colev(i)   = AVER(traj1(k1)%cold_colev(i)   ,traj2(k2)%cold_colev(i)   ,traj2(k2+1)%cold_colev(i)   ,dum)
          trajf(k1)%cdleft_colev(i) = AVER(traj1(k1)%cdleft_colev(i) ,traj2(k2)%cdleft_colev(i) ,traj2(k2+1)%cdleft_colev(i) ,dum)
       ENDDO
       DO i = 1, ntrhat 
          trajf(k1)%emis_h(i)       = AVER(traj1(k1)%emis_h(i)       ,traj2(k2)%emis_h(i)       ,traj2(k2+1)%emis_h(i)       ,dum)
          trajf(k1)%inta_h(i)       = AVER(traj1(k1)%inta_h(i)       ,traj2(k2)%inta_h(i)       ,traj2(k2+1)%inta_h(i)       ,dum)
       ENDDO
       DO i = 1, ntrcat 
          trajf(k1)%emis_c(i)       = AVER(traj1(k1)%emis_c(i)       ,traj2(k2)%emis_c(i)       ,traj2(k2+1)%emis_c(i)       ,dum)
          trajf(k1)%inta_c(i)       = AVER(traj1(k1)%inta_c(i)       ,traj2(k2)%inta_c(i)       ,traj2(k2+1)%inta_c(i)       ,dum)
       ENDDO
       DO i = 1, ntrnat 
          trajf(k1)%emis_n(i)       = AVER(traj1(k1)%emis_n(i)       ,traj2(k2)%emis_n(i)       ,traj2(k2+1)%emis_n(i)       ,dum)
          trajf(k1)%inta_n(i)       = AVER(traj1(k1)%inta_n(i)       ,traj2(k2)%inta_n(i)       ,traj2(k2+1)%inta_n(i)       ,dum)
       ENDDO
       DO i = 1, ntroat 
          trajf(k1)%emis_o(i)       = AVER(traj1(k1)%emis_o(i)       ,traj2(k2)%emis_o(i)       ,traj2(k2+1)%emis_o(i)       ,dum)
          trajf(k1)%inta_o(i)       = AVER(traj1(k1)%inta_o(i)       ,traj2(k2)%inta_o(i)       ,traj2(k2+1)%inta_o(i)       ,dum)
       ENDDO
       DO i = 1, ntrsat 
          trajf(k1)%emis_s(i)       = AVER(traj1(k1)%emis_s(i)       ,traj2(k2)%emis_s(i)       ,traj2(k2+1)%emis_s(i)       ,dum)
          trajf(k1)%inta_s(i)       = AVER(traj1(k1)%inta_s(i)       ,traj2(k2)%inta_s(i)       ,traj2(k2+1)%inta_s(i)       ,dum)
       ENDDO
       DO i = 1, ntrsiat
          trajf(k1)%emis_si(i)      = AVER(traj1(k1)%emis_si(i)      ,traj2(k2)%emis_si(i)      ,traj2(k2+1)%emis_si(i)      ,dum)
          trajf(k1)%inta_si(i)      = AVER(traj1(k1)%inta_si(i)      ,traj2(k2)%inta_si(i)      ,traj2(k2+1)%inta_si(i)      ,dum)
       ENDDO
       DO i = 1, ntrcpl 
          trajf(k1)%emis_cp(i)      = AVER(traj1(k1)%emis_cp(i)      ,traj2(k2)%emis_cp(i)      ,traj2(k2+1)%emis_cp(i)      ,dum)
          trajf(k1)%inta_cp(i)      = AVER(traj1(k1)%inta_cp(i)      ,traj2(k2)%inta_cp(i)      ,traj2(k2+1)%inta_cp(i)      ,dum)
       ENDDO
       DO i = 1, ntrnpl 
          trajf(k1)%emis_np(i)      = AVER(traj1(k1)%emis_np(i)      ,traj2(k2)%emis_np(i)      ,traj2(k2+1)%emis_np(i)      ,dum)
          trajf(k1)%inta_np(i)      = AVER(traj1(k1)%inta_np(i)      ,traj2(k2)%inta_np(i)      ,traj2(k2+1)%inta_np(i)      ,dum)
       ENDDO
       DO i = 1, ntropl 
          trajf(k1)%emis_op(i)      = AVER(traj1(k1)%emis_op(i)      ,traj2(k2)%emis_op(i)      ,traj2(k2+1)%emis_op(i)      ,dum)
          trajf(k1)%inta_op(i)      = AVER(traj1(k1)%inta_op(i)      ,traj2(k2)%inta_op(i)      ,traj2(k2+1)%inta_op(i)      ,dum)
       ENDDO
       DO i = 1, ntrspl 
          trajf(k1)%emis_sp(i)      = AVER(traj1(k1)%emis_sp(i)      ,traj2(k2)%emis_sp(i)      ,traj2(k2+1)%emis_sp(i)      ,dum)
          trajf(k1)%inta_sp(i)      = AVER(traj1(k1)%inta_sp(i)      ,traj2(k2)%inta_sp(i)      ,traj2(k2+1)%inta_sp(i)      ,dum)
       ENDDO
       DO i = 1, ntrsipl
          trajf(k1)%emis_sip(i)     = AVER(traj1(k1)%emis_sip(i)     ,traj2(k2)%emis_sip(i)     ,traj2(k2+1)%emis_sip(i)     ,dum)
          trajf(k1)%inta_sip(i)     = AVER(traj1(k1)%inta_sip(i)     ,traj2(k2)%inta_sip(i)     ,traj2(k2+1)%inta_sip(i)     ,dum)
       ENDDO
       
    ENDDO

  END SUBROUTINE AVERAGE_TRAJEC


  SUBROUTINE TEST_CROSS_SONIC(imax, eps, success)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     test if we can cross the sonic point or not
    ! subroutine/function needed :
    ! input variables :
    !     imax = maximal index of input trajectory
    !     eps  = relative precision on variables
    ! output variables :
    !     success = boolean indicating if we can or not cross the sonic point
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_PROFIL_TABLES

    IMPLICIT none

    INTEGER,       INTENT(in)                :: imax
    REAL(KIND=DP), INTENT(in)                :: eps
    LOGICAL,       INTENT(out)               :: success
    REAL(KIND=DP)                            :: dmin
    REAL(KIND=DP)                            :: dmax
    INTEGER                                  :: imin
    INTEGER                                  :: i

    ! Variables for regression
    INTEGER                                  :: N
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: X
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: Y
    REAL(KIND=DP)                            :: a
    REAL(KIND=DP)                            :: b
    REAL(KIND=DP)                            :: rsq

    ! -------------------------------------------
    ! Find the range of indices over
    ! wich the fits are computed
    ! -------------------------------------------
    dmax = traj_main(imax)%distance
    imin = imax
    dmin = traj_main(imin)%distance
    DO WHILE( ABS(dmax-dmin) / dmax < eps )
       imin = imin - 1
       dmin = traj_main(imin)%distance
    ENDDO

    ! -------------------------------------------
    ! Allocate X and Y tables, fill X table
    ! -------------------------------------------
    N = imax - imin + 1
    ALLOCATE (X(N))
    ALLOCATE (Y(N))
    DO i = imin, imax
       X(i-imin+1) = traj_main(i)%distance
    ENDDO

    ! -------------------------------------------
    ! Fit neutral & ion velocities & temperatures
    ! Y = a X + B
    ! -------------------------------------------
    DO i = imin, imax
       Y(i-imin+1) = traj_main(i)%Vsound - traj_main(i)%Vn
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)

    ! -------------------------------------------
    ! If the slope is negative, we can cross
    ! If not, well ... we can't
    ! -------------------------------------------
    IF ( a < 0.0_dp ) THEN
       success = .true.
    ELSE
       success = .false.
    ENDIF

    ! -------------------------------------------
    ! Deallocate X and Y tables
    ! -------------------------------------------
    DEALLOCATE (X)
    DEALLOCATE (Y)
  END SUBROUTINE TEST_CROSS_SONIC


  SUBROUTINE CROSS_SONIC_POINT(imax, eps, dist_cross, step_cross)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     linearize main variables along the trajectory in order
    !     to cross a sonic point
    ! subroutine/function needed :
    ! input variables :
    !     imax = maximal index of input trajectory
    !     eps  = relative precision on variables
    ! output variables :
    !     dist_cross = distance at which variables are extrapolated
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS
    USE MODULE_PHYS_VAR
    USE MODULE_GAMMA
    USE MODULE_PROFIL_TABLES
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_H2
    USE MODULE_CO

    IMPLICIT none

    INTEGER,       INTENT(in)                :: imax
    REAL(KIND=DP), INTENT(in)                :: eps
    REAL(KIND=DP), INTENT(in)                :: step_cross
    REAL(KIND=DP), INTENT(out)               :: dist_cross
    REAL(KIND=DP)                            :: dmin
    REAL(KIND=DP)                            :: dmax
    INTEGER                                  :: imin
    INTEGER                                  :: i
    INTEGER                                  :: j

    ! Variables for regression
    INTEGER                                  :: N
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: X
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: Y
    REAL(KIND=DP)                            :: a
    REAL(KIND=DP)                            :: b
    REAL(KIND=DP)                            :: rsq

    ! regression coefficients
    REAL(KIND=DP)                            :: a_Vn, b_Vn, r2_Vn
    REAL(KIND=DP)                            :: a_Vi, b_Vi, r2_Vi
    REAL(KIND=DP)                            :: a_Tn, b_Tn, r2_Tn
    REAL(KIND=DP)                            :: a_Ti, b_Ti, r2_Ti

    ! dummy variables
    REAL(KIND=DP)                            :: cte
    REAL(KIND=DP)                            :: dum

    ! -------------------------------------------
    ! Find the range of indices over
    ! wich the fits are computed
    ! -------------------------------------------
    dmax = traj_main(imax)%distance
    imin = imax
    dmin = traj_main(imin)%distance
    DO WHILE( ABS(dmax-dmin) / dmax < eps )
       imin = imin - 1
       dmin = traj_main(imin)%distance
    ENDDO

    ! -------------------------------------------
    ! Allocate X and Y tables, fill X table
    ! -------------------------------------------
    N = imax - imin + 1
    ALLOCATE (X(N))
    ALLOCATE (Y(N))
    DO i = imin, imax
       X(i-imin+1) = LOG(traj_main(i)%distance)
    ENDDO

    ! -------------------------------------------
    ! Fit neutral & ion velocities & temperatures
    ! Y = a X + B
    ! -------------------------------------------
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Vn)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    a_Vn  = a
    b_Vn  = b
    r2_Vn = Rsq
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Vi)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    a_Vi  = a
    b_Vi  = b
    r2_Vi = Rsq
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Tn)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    a_Tn  = a
    b_Tn  = b
    r2_Tn = Rsq
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Ti)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    a_Ti  = a
    b_Ti  = b
    r2_Ti = Rsq

    ! -------------------------------------------
    ! Find distance where Vn crosses Cs
    ! -------------------------------------------
    cte = Gamma * kB / muN
    dum = EXP(b_Vn) / ( cte * EXP(b_Tn) )**0.5_dp
    dist_cross = dum ** ( 1.0_dp / (a_Tn / 2.0_dp - a_Vn) )
    IF( (dist_cross - dmax) < 0.0_dp .OR. (dist_cross - dmax) / dmax > 1.0_dp ) THEN
       WRITE(*,*) "There is a problem in the crossing of sonic point"
       WRITE(*,*) "-> code stops, please advise"
       STOP
    ENDIF

    ! -------------------------------------------
    ! Compute Vn, Vi, Tn, and Ti at 
    ! dist_cross + (dist_cross - dmax)
    ! Deduce all other variables
    ! We assume no source terms to keep things
    ! as simple as possible
    ! -------------------------------------------
    dist_cross = dist_cross + step_cross * (dist_cross - dmax)
    Vn         = EXP(b_Vn + a_Vn * LOG(dist_cross))
    Vi         = EXP(b_Vi + a_Vi * LOG(dist_cross))
    Tn         = EXP(b_Tn + a_Tn * LOG(dist_cross))
    Ti         = EXP(b_Ti + a_Ti * LOG(dist_cross))
    ! Ti         = Tn
    ! Te         = Ti
    ! rhoN       = traj_main(imax)%rhoN       * traj_main(imax)%Vn      / Vn
    ! rhoI       = traj_main(imax)%rhoI       * traj_main(imax)%Vi      / Vi
    ! rhoA       = traj_main(imax)%rhoA       * traj_main(imax)%Vi      / Vi
    ! rhoNeg     = traj_main(imax)%rhoNeg     * traj_main(imax)%Vi      / Vi
    ! DensityN   = traj_main(imax)%DensityN   * rhoN   / traj_main(imax)%rhoN
    ! DensityI   = traj_main(imax)%DensityI   * rhoI   / traj_main(imax)%rhoI
    ! DensityA   = traj_main(imax)%DensityA   * rhoA   / traj_main(imax)%rhoA
    ! DensityNeg = traj_main(imax)%DensityNeg * rhoNeg / traj_main(imax)%rhoNeg
    ! V_CH       = Vn
    ! V_S        = Vn
    ! V_SH       = Vn
    ! DO i = b_neu, e_neu
    !    speci(i)%density = speci(i)%density * rhoN / traj_main(imax)%rhoN
    ! ENDDO
    ! DO i = b_ion, e_ion
    !    speci(i)%density = speci(i)%density * rhoI / traj_main(imax)%rhoI
    ! ENDDO
    ! DO i = b_ani, e_ani
    !    speci(i)%density = speci(i)%density * rhoI / traj_main(imax)%rhoI
    ! ENDDO
    ! DO i = b_neg, e_neg
    !    speci(i)%density = speci(i)%density * rhoI / traj_main(imax)%rhoI
    ! ENDDO
    ! DO i = b_gra, e_gra
    !    speci(i)%density = speci(i)%density * rhoI / traj_main(imax)%rhoI
    ! ENDDO
    ! DO i = b_cor, e_cor
    !    speci(i)%density = speci(i)%density * rhoI / traj_main(imax)%rhoI
    ! ENDDO
    ! nH = speci(ind_H)%density + 2._DP * speci(ind_H2)%density + speci(ind_Hplus)%density
    ! DO i = 1, NH2_lev
    !    H2_lev(i)%density = H2_lev(i)%density * rhoN / traj_main(imax)%rhoN
    ! ENDDO
    ! DO i = 1, NCO_lev
    !    CO_lev(i)%density = CO_lev(i)%density * rhoN / traj_main(imax)%rhoN
    ! ENDDO
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Vn)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    Vn = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Vi)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    Vi = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%RhoN)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    RhoN = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%RhoI)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    RhoI = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%RhoA)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    RhoA = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%RhoNeg)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    RhoNeg = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Tn)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    Tn = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Ti)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    Ti = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%Te)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    Te = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%DensityN)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    DensityN = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%DensityI)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    DensityI = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%DensityA)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    DensityA = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%DensityNeg)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    DensityNeg = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%v_CH)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    v_CH = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%v_S)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    v_S = EXP(b + a * LOG(dist_cross))
    DO i = imin, imax
       Y(i-imin+1) = LOG(traj_main(i)%v_SH)
    ENDDO
    CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    v_SH = EXP(b + a * LOG(dist_cross))
    DO j = 1, Nspec
       DO i = imin, imax
          Y(i-imin+1) = LOG(traj_main(i)%abon(j))
       ENDDO
       CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
       speci(j)%density = EXP(b + a * LOG(dist_cross))
    ENDDO
    DO j = 1, NH2_lev
       DO i = imin, imax
          Y(i-imin+1) = LOG(traj_main(i)%abon_h2lev(j))
       ENDDO
       CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
       H2_lev(j)%density = EXP(b + a * LOG(dist_cross))
    ENDDO
    DO j = 1, NCO_lev
       DO i = imin, imax
          Y(i-imin+1) = LOG(traj_main(i)%abon_colev(j))
       ENDDO
       CALL LINEAR_REGRESSION(N, X, Y, a, b, rsq)
       CO_lev(j)%density = EXP(b + a * LOG(dist_cross))
    ENDDO

    ! -------------------------------------------
    ! reset dvode variables
    ! -------------------------------------------
    v_variab(iv_Vn        ) = Vn
    v_variab(iv_Vi        ) = Vi
    v_variab(iv_RhoN      ) = RhoN
    v_variab(iv_RhoI      ) = RhoI
    v_variab(iv_RhoA      ) = RhoA
    v_variab(iv_RhoNEG    ) = RhoNEG
    v_variab(iv_Tn        ) = Tn
    v_variab(iv_Ti        ) = Ti
    v_variab(iv_Te        ) = Te
    v_variab(iv_DensityN  ) = DensityN
    v_variab(iv_DensityI  ) = DensityI
    v_variab(iv_DensityA  ) = DensityA
    v_variab(iv_DensityNeg) = DensityNeg
    v_variab(iv_vCH       ) = v_CH
    v_variab(iv_vS        ) = v_S
    v_variab(iv_vSH       ) = v_SH
    v_variab(bv_speci:ev_speci) = speci(1:Nspec)%density
    v_variab(bv_H2_lev:ev_H2_lev) = H2_lev(1:NH2_lev_var)%density
    WHERE (v_variab > 0)
       v_lvariab = LOG(v_variab)
    ELSEWHERE
       v_lvariab = minus_infinity
    END WHERE

    ! -------------------------------------------
    ! Deallocate X and Y tables
    ! -------------------------------------------
    DEALLOCATE (X)
    DEALLOCATE (Y)
  END SUBROUTINE CROSS_SONIC_POINT


  SUBROUTINE LINEAR_REGRESSION(N, X, Y, a, b, rsq)
    !---------------------------------------------------------------------------
    ! called by :
    !     CROSS_SONIC_POINT
    ! purpose :
    !     perform a simple linear regression of tables X and Y
    !     Y = a X + b
    ! subroutine/function needed :
    ! input variables :
    !     N    = number of points
    !     X, Y = tables to fit
    ! output variables :
    !     a, b = fits coefficients
    !     rsq  = determination coefficient
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT none

    INTEGER,                       INTENT(in)  :: N
    REAL(KIND=DP), DIMENSION(1:N), INTENT(in)  :: X
    REAL(KIND=DP), DIMENSION(1:N), INTENT(in)  :: Y
    REAL(KIND=DP),                 INTENT(out) :: a
    REAL(KIND=DP),                 INTENT(out) :: b
    REAL(KIND=DP),                 INTENT(out) :: rsq

    REAL(KIND=DP)                              :: sx
    REAL(KIND=DP)                              :: sxx
    REAL(KIND=DP)                              :: sy
    REAL(KIND=DP)                              :: sxy
    REAL(KIND=DP)                              :: delta
    REAL(KIND=DP)                              :: s

    s   = DBLE(N)
    sx  = SUM( X(1:N) )
    sy  = SUM( Y(1:N) )
    sxx = SUM( X(1:N)**2.0_dp )
    sxy = SUM( X(1:N)*Y(1:N) )

    delta = s * sxx - sx**2.0_dp
    b = ( sxx * sy  - sx * sxy ) / delta
    a = ( s   * sxy - sx * sy  ) / delta

    rsq = 0.0_dp
  END SUBROUTINE LINEAR_REGRESSION


  REAL(KIND=DP) FUNCTION AVER(a,b1,b2,s)
    !---------------------------------------------------------------------------
    ! called by :
    !     AVERAGE_TRAJEC
    ! purpose :
    !     perform a simple average of a and (b1 + (b2 - b1) * s)
    ! subroutine/function needed :
    ! input variables :
    !     a, b1, b2, s
    ! output variables :
    ! results :
    !     averaged value
    !---------------------------------------------------------------------------
    IMPLICIT none
    REAL(KIND=DP), INTENT(in) :: a
    REAL(KIND=DP), INTENT(in) :: b1
    REAL(KIND=DP), INTENT(in) :: b2
    REAL(KIND=DP), INTENT(in) :: s
    AVER = (a + b1 + (b2 - b1) * s) / 2.0_dp
  END FUNCTION AVER


  SUBROUTINE ESTIMATE_DECOUPLING(imax, eps, aver)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     compute average value of the decoupling over a short 
    !     distance (d*eps) along the main trajectory
    ! subroutine/function needed :
    ! input variables :
    !     imax = maximal index of input trajectory
    !     eps  = relative precision on variables
    ! output variables :
    !     aver = decoupling average
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_PROFIL_TABLES

    IMPLICIT none

    INTEGER,       INTENT(in)                :: imax
    REAL(KIND=DP), INTENT(in)                :: eps
    REAL(KIND=DP), INTENT(out)               :: aver
    REAL(KIND=DP)                            :: dmin
    REAL(KIND=DP)                            :: dmax
    INTEGER                                  :: imin
    INTEGER                                  :: i

    ! Variables for computing average
    INTEGER                                  :: N
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: X
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: Y
    REAL(KIND=DP)                            :: intg

    ! -------------------------------------------
    ! return 0 if we are at the first point of 
    ! the trajectory
    ! -------------------------------------------
    IF ( imax == 0 ) THEN
       aver = 0.0_dp
       RETURN
    ENDIF

    ! -------------------------------------------
    ! Find the range of indices over
    ! wich the average is computed
    ! -------------------------------------------
    dmax = traj_main(imax)%distance
    imin = imax
    dmin = traj_main(imin)%distance
    DO WHILE( ABS(dmax-dmin) / dmax < eps .AND. imin > 1 )
       imin = imin - 1
       dmin = traj_main(imin)%distance
    ENDDO

    ! -------------------------------------------
    ! Allocate X and Y tables, fill X and Y table
    ! -------------------------------------------
    N = imax - imin + 1
    ALLOCATE (X(N))
    ALLOCATE (Y(N))
    DO i = imin, imax
       X(i-imin+1) = traj_main(i)%distance
       Y(i-imin+1) = traj_main(i)%Vn - traj_main(i)%Vi
    ENDDO

    ! -------------------------------------------
    ! Compute average of decoupling
    ! -------------------------------------------
    IF ( imin /= imax ) THEN
       intg = INTEG_TABLE(N, X, Y)
       aver = intg / (dmax - dmin)
    ELSE
       aver = Y(imin)
    ENDIF

    ! -------------------------------------------
    ! Deallocate X and Y tables
    ! -------------------------------------------
    DEALLOCATE (X)
    DEALLOCATE (Y)
  END SUBROUTINE ESTIMATE_DECOUPLING


  REAL(KIND=DP) FUNCTION INTEG_TABLE(N,X,Y)
    !---------------------------------------------------------------------------
    ! called by :
    !     ESTIMATE_DECOUPLING
    ! purpose :
    !     perform a simple integral of Y over variable X
    ! subroutine/function needed :
    ! input variables :
    !     N = number of points
    !     X and Y = tables
    ! output variables :
    ! results :
    !     integral value
    !---------------------------------------------------------------------------
    IMPLICIT none
    INTEGER,                       INTENT(in) :: N
    REAL(KIND=DP), DIMENSION(1:N), INTENT(in) :: X
    REAL(KIND=DP), DIMENSION(1:N), INTENT(in) :: Y
    INTEGER                                   :: i
    INTEG_TABLE = 0.0_dp
    INTEG_TABLE = INTEG_TABLE + Y(1) * ( X(2) - X(1)   ) / 2.0_dp
    INTEG_TABLE = INTEG_TABLE + Y(N) * ( X(N) - X(N-1) ) / 2.0_dp
    DO i = 2, N-1
       INTEG_TABLE = INTEG_TABLE + Y(i) * ( X(i+1) - X(i-1) ) / 2.0_dp
    ENDDO
  END FUNCTION INTEG_TABLE

END MODULE MODULE_SWITCH_CJ
