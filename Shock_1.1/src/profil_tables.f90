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

MODULE MODULE_PROFIL_TABLES

  !*****************************************************************************
  !** The module 'MODULE_PROFIL_TABLES' contains all                          **
  !** the tables declaration and associated function to save the shock profil **
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! ----------------------------------------------------------------------------
  ! data type for shock trajectory
  ! ----------------------------------------------------------------------------
  TYPE TYPE_TRAJECTORY
     ! -------------------------------------------------------------------------
     ! time, distance, density, velocity, temperature
     ! -------------------------------------------------------------------------
     INTEGER       :: slab
     REAL(KIND=DP) :: distance
     REAL(KIND=DP) :: dist_step
     REAL(KIND=DP) :: Av
     REAL(KIND=DP) :: timeN
     REAL(KIND=DP) :: timeI
     REAL(KIND=DP) :: nH
     REAL(KIND=DP) :: Vsound
     REAL(KIND=DP) :: Vmagnet
     REAL(KIND=DP) :: Valfven
     REAL(KIND=DP) :: rhoN
     REAL(KIND=DP) :: rhoI
     REAL(KIND=DP) :: rhoA
     REAL(KIND=DP) :: rhoNeg
     REAL(KIND=DP) :: DensityN
     REAL(KIND=DP) :: DensityI
     REAL(KIND=DP) :: DensityA
     REAL(KIND=DP) :: DensityNeg
     REAL(KIND=DP) :: muN
     REAL(KIND=DP) :: muI
     REAL(KIND=DP) :: muA
     REAL(KIND=DP) :: muNeg
     REAL(KIND=DP) :: Tn
     REAL(KIND=DP) :: Ti
     REAL(KIND=DP) :: Te
     REAL(KIND=DP) :: Vn
     REAL(KIND=DP) :: Vi
     REAL(KIND=DP) :: dVn
     REAL(KIND=DP) :: dVi
     REAL(KIND=DP) :: Vgrad
     REAL(KIND=DP) :: grad_V
     REAL(KIND=DP) :: v_CH
     REAL(KIND=DP) :: v_S
     REAL(KIND=DP) :: v_SH
     REAL(KIND=DP) :: iondeg
     REAL(KIND=DP) :: ionfrac
     REAL(KIND=DP) :: oph2
     REAL(KIND=DP) :: T_gr
     REAL(KIND=DP) :: Teff_gr
     REAL(KIND=DP) :: nlay_gr
     REAL(KIND=DP) :: n_gr
     REAL(KIND=DP) :: mu_gr
     REAL(KIND=DP) :: much_gr
     REAL(KIND=DP) :: r_gr
     REAL(KIND=DP) :: m_gr
     REAL(KIND=DP) :: m_grc
     REAL(KIND=DP) :: m_grm
     REAL(KIND=DP) :: d_ero
     REAL(KIND=DP) :: d_ads
     REAL(KIND=DP) :: rsq_grc
     REAL(KIND=DP) :: rsq_grm
     REAL(KIND=DP) :: compr_n
     REAL(KIND=DP) :: compr_i

     ! -------------------------------------------------------------------------
     ! heating and cooling
     ! -------------------------------------------------------------------------
     REAL(KIND=DP) :: line_rad_tot
     REAL(KIND=DP) :: line_rad_n
     REAL(KIND=DP) :: line_rad_i
     REAL(KIND=DP) :: line_rad_e
     REAL(KIND=DP) :: line_rad_n_molec
     REAL(KIND=DP) :: line_rad_n_H2
     REAL(KIND=DP) :: line_rad_n_13CO
     REAL(KIND=DP) :: line_rad_n_OH
     REAL(KIND=DP) :: line_rad_n_NH3
     REAL(KIND=DP) :: line_rad_n_rCO
     REAL(KIND=DP) :: line_rad_n_vCO
     REAL(KIND=DP) :: line_rad_n_CO
     REAL(KIND=DP) :: line_rad_n_roH2O
     REAL(KIND=DP) :: line_rad_n_rpH2O
     REAL(KIND=DP) :: line_rad_n_vH2O
     REAL(KIND=DP) :: line_rad_n_H2O
     REAL(KIND=DP) :: line_rad_i_molec
     REAL(KIND=DP) :: line_rad_i_H2
     REAL(KIND=DP) :: line_rad_e_molec
     REAL(KIND=DP) :: line_rad_e_H2
     REAL(KIND=DP) :: line_rad_e_CO
     REAL(KIND=DP) :: line_rad_n_atoms
     REAL(KIND=DP) :: line_rad_n_Hat
     REAL(KIND=DP) :: line_rad_n_Oat
     REAL(KIND=DP) :: line_rad_n_Op
     REAL(KIND=DP) :: line_rad_n_Cat
     REAL(KIND=DP) :: line_rad_n_Cp
     REAL(KIND=DP) :: line_rad_n_Nat
     REAL(KIND=DP) :: line_rad_n_Np
     REAL(KIND=DP) :: line_rad_n_Sat
     REAL(KIND=DP) :: line_rad_n_Sp
     REAL(KIND=DP) :: line_rad_n_Siat
     REAL(KIND=DP) :: line_rad_n_Sip
     REAL(KIND=DP) :: line_rad_i_atoms
     REAL(KIND=DP) :: line_rad_i_Hat
     REAL(KIND=DP) :: line_rad_i_Oat
     REAL(KIND=DP) :: line_rad_i_Op
     REAL(KIND=DP) :: line_rad_i_Cat
     REAL(KIND=DP) :: line_rad_i_Cp
     REAL(KIND=DP) :: line_rad_i_Nat
     REAL(KIND=DP) :: line_rad_i_Np
     REAL(KIND=DP) :: line_rad_i_Sat
     REAL(KIND=DP) :: line_rad_i_Sp
     REAL(KIND=DP) :: line_rad_i_Siat
     REAL(KIND=DP) :: line_rad_i_Sip
     REAL(KIND=DP) :: line_rad_e_atoms
     REAL(KIND=DP) :: line_rad_e_Hat
     REAL(KIND=DP) :: line_rad_e_Oat
     REAL(KIND=DP) :: line_rad_e_Op
     REAL(KIND=DP) :: line_rad_e_Cat
     REAL(KIND=DP) :: line_rad_e_Cp
     REAL(KIND=DP) :: line_rad_e_Nat
     REAL(KIND=DP) :: line_rad_e_Np
     REAL(KIND=DP) :: line_rad_e_Sat
     REAL(KIND=DP) :: line_rad_e_Sp
     REAL(KIND=DP) :: line_rad_e_Siat
     REAL(KIND=DP) :: line_rad_e_Sip
     REAL(KIND=DP) :: chem_DE_tot
     REAL(KIND=DP) :: chem_DE_n
     REAL(KIND=DP) :: chem_DE_i
     REAL(KIND=DP) :: chem_DE_e
     REAL(KIND=DP) :: chem_DE_n_diss_H2
     REAL(KIND=DP) :: chem_DE_i_diss_H2
     REAL(KIND=DP) :: chem_DE_e_diss_H2
     REAL(KIND=DP) :: chem_DE_n_phgas
     REAL(KIND=DP) :: chem_DE_i_phgas
     REAL(KIND=DP) :: chem_DE_e_phgas
     REAL(KIND=DP) :: chem_DE_n_phgrn
     REAL(KIND=DP) :: chem_DE_i_phgrn
     REAL(KIND=DP) :: chem_DE_e_phgrn
     REAL(KIND=DP) :: chem_DE_n_cosmic
     REAL(KIND=DP) :: chem_DE_i_cosmic
     REAL(KIND=DP) :: chem_DE_e_cosmic
     REAL(KIND=DP) :: chem_DE_n_other
     REAL(KIND=DP) :: chem_DE_i_other
     REAL(KIND=DP) :: chem_DE_e_other
     REAL(KIND=DP) :: elast_scat_tot
     REAL(KIND=DP) :: elast_scat_n
     REAL(KIND=DP) :: elast_scat_i
     REAL(KIND=DP) :: elast_scat_e
     REAL(KIND=DP) :: exch_eint_tot
     REAL(KIND=DP) :: exch_eint_n
     REAL(KIND=DP) :: exch_eint_i
     REAL(KIND=DP) :: exch_eint_e
     REAL(KIND=DP) :: therm_grain_tot
     REAL(KIND=DP) :: mech_trsf_tot
     REAL(KIND=DP) :: mech_trsf_n
     REAL(KIND=DP) :: mech_trsf_i
     REAL(KIND=DP) :: mech_trsf_e
     REAL(KIND=DP) :: mech_trsf_n_chem
     REAL(KIND=DP) :: mech_trsf_i_chem
     REAL(KIND=DP) :: mech_trsf_e_chem
     REAL(KIND=DP) :: mech_trsf_n_compr
     REAL(KIND=DP) :: mech_trsf_i_compr
     REAL(KIND=DP) :: mech_trsf_e_compr
     REAL(KIND=DP) :: mech_trsf_n_visc
     REAL(KIND=DP) :: mech_trsf_i_visc
     REAL(KIND=DP) :: mech_trsf_e_visc

     ! -------------------------------------------------------------------------
     ! energetics
     ! -------------------------------------------------------------------------
     REAL(KIND=DP) :: mom_flux_tot
     REAL(KIND=DP) :: mom_flux_kin
     REAL(KIND=DP) :: mom_flux_the
     REAL(KIND=DP) :: mom_flux_mag
     REAL(KIND=DP) :: mom_flux_vis
     REAL(KIND=DP) :: nrj_flux_tot
     REAL(KIND=DP) :: nrj_flux_kin
     REAL(KIND=DP) :: nrj_flux_the
     REAL(KIND=DP) :: nrj_flux_mag
     REAL(KIND=DP) :: nrj_flux_vis
     REAL(KIND=DP) :: nrj_flux_int
     REAL(KIND=DP) :: nrj_flux_src
     REAL(KIND=DP) :: nrj_flux_pho
     REAL(KIND=DP) :: mag_flux_corr

     ! -------------------------------------------------------------------------
     ! radiation
     ! -------------------------------------------------------------------------
     REAL(KIND=DP) :: fluph
     REAL(KIND=DP) :: coldens_h
     REAL(KIND=DP) :: coldens_h2
     REAL(KIND=DP) :: coldens_co

     ! -------------------------------------------------------------------------
     ! chemical profiles
     ! -------------------------------------------------------------------------
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: abon
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: coldens
     REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: all_rate_chem  ! reactions rates x abundances profiles (cm-3 s-1)

     ! -------------------------------------------------------------------------
     ! excitation profiles
     ! -------------------------------------------------------------------------
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: abon_h2lev
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: cold_h2lev
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: cdleft_h2lev
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_h2lin
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_h2lin
     ! CO
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: abon_colev
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: cold_colev
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: cdleft_colev
     ! atomic lines
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_h
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_h
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_c
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_c
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_n
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_n
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_o
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_o
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_s
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_s
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_si
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_si
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_cp
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_cp
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_np
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_np
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_op
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_op
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_sp
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_sp
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: emis_sip
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: inta_sip

  END TYPE TYPE_TRAJECTORY


  TYPE (TYPE_TRAJECTORY), ALLOCATABLE, DIMENSION(:) :: traj_main
  TYPE (TYPE_TRAJECTORY), ALLOCATABLE, DIMENSION(:) :: traj_high
  TYPE (TYPE_TRAJECTORY), ALLOCATABLE, DIMENSION(:) :: traj_down
  TYPE (TYPE_TRAJECTORY), ALLOCATABLE, DIMENSION(:) :: traj_curr

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE ALLOCATE_TABLES
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     allocates tables of physical variables to save the shock profil
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------

    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    ! USE MODULE_CHEM_REACT
    USE MODULE_LINE_EXCIT

    IMPLICIT NONE

    INTEGER :: i

    ALLOCATE (traj_main(1:Nstep_max) )
    ALLOCATE (traj_high(1:Nstep_max) )
    ALLOCATE (traj_down(1:Nstep_max) )
    ALLOCATE (traj_curr(1:Nstep_max) )

    DO i = 1, Nstep_max
       ALLOCATE (traj_main(i)%abon(1:Nspec)              )
       ALLOCATE (traj_main(i)%coldens(1:Nspec)           )
       ! ALLOCATE (traj_main(i)%all_rate_chem(1:Nreact)    ) ! to save memory
       ALLOCATE (traj_main(i)%abon_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_main(i)%cold_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_main(i)%cdleft_h2lev(1:NH2_lev)    )
       ALLOCATE (traj_main(i)%emis_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_main(i)%inta_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_main(i)%abon_colev(1:NCO_lev)      )
       ALLOCATE (traj_main(i)%cold_colev(1:NCO_lev)      )
       ALLOCATE (traj_main(i)%cdleft_colev(1:NCO_lev)    )
       ALLOCATE (traj_main(i)%emis_h  (1:ntrhat ))
       ALLOCATE (traj_main(i)%inta_h  (1:ntrhat ))
       ALLOCATE (traj_main(i)%emis_c  (1:ntrcat ))
       ALLOCATE (traj_main(i)%inta_c  (1:ntrcat ))
       ALLOCATE (traj_main(i)%emis_n  (1:ntrnat ))
       ALLOCATE (traj_main(i)%inta_n  (1:ntrnat ))
       ALLOCATE (traj_main(i)%emis_o  (1:ntroat ))
       ALLOCATE (traj_main(i)%inta_o  (1:ntroat ))
       ALLOCATE (traj_main(i)%emis_s  (1:ntrsat ))
       ALLOCATE (traj_main(i)%inta_s  (1:ntrsat ))
       ALLOCATE (traj_main(i)%emis_si (1:ntrsiat))
       ALLOCATE (traj_main(i)%inta_si (1:ntrsiat))
       ALLOCATE (traj_main(i)%emis_cp (1:ntrcpl ))
       ALLOCATE (traj_main(i)%inta_cp (1:ntrcpl ))
       ALLOCATE (traj_main(i)%emis_np (1:ntrnpl ))
       ALLOCATE (traj_main(i)%inta_np (1:ntrnpl ))
       ALLOCATE (traj_main(i)%emis_op (1:ntropl ))
       ALLOCATE (traj_main(i)%inta_op (1:ntropl ))
       ALLOCATE (traj_main(i)%emis_sp (1:ntrspl ))
       ALLOCATE (traj_main(i)%inta_sp (1:ntrspl ))
       ALLOCATE (traj_main(i)%emis_sip(1:ntrsipl))
       ALLOCATE (traj_main(i)%inta_sip(1:ntrsipl))

       ALLOCATE (traj_high(i)%abon(1:Nspec)              )
       ALLOCATE (traj_high(i)%coldens(1:Nspec)           )
       ! ALLOCATE (traj_high(i)%all_rate_chem(1:Nreact)    ) ! to save memory
       ALLOCATE (traj_high(i)%abon_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_high(i)%cold_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_high(i)%cdleft_h2lev(1:NH2_lev)    )
       ALLOCATE (traj_high(i)%emis_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_high(i)%inta_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_high(i)%abon_colev(1:NCO_lev)      )
       ALLOCATE (traj_high(i)%cold_colev(1:NCO_lev)      )
       ALLOCATE (traj_high(i)%cdleft_colev(1:NCO_lev)    )
       ALLOCATE (traj_high(i)%emis_h  (1:ntrhat ))
       ALLOCATE (traj_high(i)%inta_h  (1:ntrhat ))
       ALLOCATE (traj_high(i)%emis_c  (1:ntrcat ))
       ALLOCATE (traj_high(i)%inta_c  (1:ntrcat ))
       ALLOCATE (traj_high(i)%emis_n  (1:ntrnat ))
       ALLOCATE (traj_high(i)%inta_n  (1:ntrnat ))
       ALLOCATE (traj_high(i)%emis_o  (1:ntroat ))
       ALLOCATE (traj_high(i)%inta_o  (1:ntroat ))
       ALLOCATE (traj_high(i)%emis_s  (1:ntrsat ))
       ALLOCATE (traj_high(i)%inta_s  (1:ntrsat ))
       ALLOCATE (traj_high(i)%emis_si (1:ntrsiat))
       ALLOCATE (traj_high(i)%inta_si (1:ntrsiat))
       ALLOCATE (traj_high(i)%emis_cp (1:ntrcpl ))
       ALLOCATE (traj_high(i)%inta_cp (1:ntrcpl ))
       ALLOCATE (traj_high(i)%emis_np (1:ntrnpl ))
       ALLOCATE (traj_high(i)%inta_np (1:ntrnpl ))
       ALLOCATE (traj_high(i)%emis_op (1:ntropl ))
       ALLOCATE (traj_high(i)%inta_op (1:ntropl ))
       ALLOCATE (traj_high(i)%emis_sp (1:ntrspl ))
       ALLOCATE (traj_high(i)%inta_sp (1:ntrspl ))
       ALLOCATE (traj_high(i)%emis_sip(1:ntrsipl))
       ALLOCATE (traj_high(i)%inta_sip(1:ntrsipl))

       ALLOCATE (traj_down(i)%abon(1:Nspec)              )
       ALLOCATE (traj_down(i)%coldens(1:Nspec)           )
       ! ALLOCATE (traj_down(i)%all_rate_chem(1:Nreact)    ) ! to save memory
       ALLOCATE (traj_down(i)%abon_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_down(i)%cold_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_down(i)%cdleft_h2lev(1:NH2_lev)    )
       ALLOCATE (traj_down(i)%emis_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_down(i)%inta_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_down(i)%abon_colev(1:NCO_lev)      )
       ALLOCATE (traj_down(i)%cold_colev(1:NCO_lev)      )
       ALLOCATE (traj_down(i)%cdleft_colev(1:NCO_lev)    )
       ALLOCATE (traj_down(i)%emis_h  (1:ntrhat ))
       ALLOCATE (traj_down(i)%inta_h  (1:ntrhat ))
       ALLOCATE (traj_down(i)%emis_c  (1:ntrcat ))
       ALLOCATE (traj_down(i)%inta_c  (1:ntrcat ))
       ALLOCATE (traj_down(i)%emis_n  (1:ntrnat ))
       ALLOCATE (traj_down(i)%inta_n  (1:ntrnat ))
       ALLOCATE (traj_down(i)%emis_o  (1:ntroat ))
       ALLOCATE (traj_down(i)%inta_o  (1:ntroat ))
       ALLOCATE (traj_down(i)%emis_s  (1:ntrsat ))
       ALLOCATE (traj_down(i)%inta_s  (1:ntrsat ))
       ALLOCATE (traj_down(i)%emis_si (1:ntrsiat))
       ALLOCATE (traj_down(i)%inta_si (1:ntrsiat))
       ALLOCATE (traj_down(i)%emis_cp (1:ntrcpl ))
       ALLOCATE (traj_down(i)%inta_cp (1:ntrcpl ))
       ALLOCATE (traj_down(i)%emis_np (1:ntrnpl ))
       ALLOCATE (traj_down(i)%inta_np (1:ntrnpl ))
       ALLOCATE (traj_down(i)%emis_op (1:ntropl ))
       ALLOCATE (traj_down(i)%inta_op (1:ntropl ))
       ALLOCATE (traj_down(i)%emis_sp (1:ntrspl ))
       ALLOCATE (traj_down(i)%inta_sp (1:ntrspl ))
       ALLOCATE (traj_down(i)%emis_sip(1:ntrsipl))
       ALLOCATE (traj_down(i)%inta_sip(1:ntrsipl))

       ALLOCATE (traj_curr(i)%abon(1:Nspec)              )
       ALLOCATE (traj_curr(i)%coldens(1:Nspec)           )
       ! ALLOCATE (traj_curr(i)%all_rate_chem(1:Nreact)    ) ! to save memory
       ALLOCATE (traj_curr(i)%abon_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_curr(i)%cold_h2lev(1:NH2_lev)      )
       ALLOCATE (traj_curr(i)%cdleft_h2lev(1:NH2_lev)    )
       ALLOCATE (traj_curr(i)%emis_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_curr(i)%inta_h2lin(1:NH2_lines_out))
       ALLOCATE (traj_curr(i)%abon_colev(1:NCO_lev)      )
       ALLOCATE (traj_curr(i)%cold_colev(1:NCO_lev)      )
       ALLOCATE (traj_curr(i)%cdleft_colev(1:NCO_lev)    )
       ALLOCATE (traj_curr(i)%emis_h  (1:ntrhat ))
       ALLOCATE (traj_curr(i)%inta_h  (1:ntrhat ))
       ALLOCATE (traj_curr(i)%emis_c  (1:ntrcat ))
       ALLOCATE (traj_curr(i)%inta_c  (1:ntrcat ))
       ALLOCATE (traj_curr(i)%emis_n  (1:ntrnat ))
       ALLOCATE (traj_curr(i)%inta_n  (1:ntrnat ))
       ALLOCATE (traj_curr(i)%emis_o  (1:ntroat ))
       ALLOCATE (traj_curr(i)%inta_o  (1:ntroat ))
       ALLOCATE (traj_curr(i)%emis_s  (1:ntrsat ))
       ALLOCATE (traj_curr(i)%inta_s  (1:ntrsat ))
       ALLOCATE (traj_curr(i)%emis_si (1:ntrsiat))
       ALLOCATE (traj_curr(i)%inta_si (1:ntrsiat))
       ALLOCATE (traj_curr(i)%emis_cp (1:ntrcpl ))
       ALLOCATE (traj_curr(i)%inta_cp (1:ntrcpl ))
       ALLOCATE (traj_curr(i)%emis_np (1:ntrnpl ))
       ALLOCATE (traj_curr(i)%inta_np (1:ntrnpl ))
       ALLOCATE (traj_curr(i)%emis_op (1:ntropl ))
       ALLOCATE (traj_curr(i)%inta_op (1:ntropl ))
       ALLOCATE (traj_curr(i)%emis_sp (1:ntrspl ))
       ALLOCATE (traj_curr(i)%inta_sp (1:ntrspl ))
       ALLOCATE (traj_curr(i)%emis_sip(1:ntrsipl))
       ALLOCATE (traj_curr(i)%inta_sip(1:ntrsipl))
    ENDDO

  END SUBROUTINE ALLOCATE_TABLES


  SUBROUTINE INITIALIZE_TABLES(N, trajec)
    !---------------------------------------------------------------------------
    ! called by :
    !     INITIALIZE
    ! purpose :
    !     initializes all components of trajec structure to 0
    ! subroutine/function needed :
    ! input variables :
    !     trajec -> the trajectory structure
    !     N      -> maximum number of points in the trajectory
    ! output variables :
    ! results :
    !     trajec -> the trajectory structure
    !---------------------------------------------------------------------------

    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    ! USE MODULE_CHEM_REACT
    USE MODULE_H2
    USE MODULE_LINE_EXCIT

    IMPLICIT NONE

    INTEGER,                                INTENT(in)    :: N
    TYPE (TYPE_TRAJECTORY), DIMENSION(1:N), INTENT(inout) :: trajec
    INTEGER                                               :: j

    DO j = 1, N
       ! -------------------------------------------------------------------------
       ! time, distance, density, velocity, temperature
       ! -------------------------------------------------------------------------
       trajec(j)%slab       = 0
       trajec(j)%distance   = 0.0_dp
       trajec(j)%dist_step  = 0.0_dp
       trajec(j)%Av         = 0.0_dp
       trajec(j)%timeN      = 0.0_dp
       trajec(j)%timeI      = 0.0_dp
       trajec(j)%nH         = 0.0_dp
       trajec(j)%Vsound     = 0.0_dp
       trajec(j)%Vmagnet    = 0.0_dp
       trajec(j)%Valfven    = 0.0_dp
       trajec(j)%rhoN       = 0.0_dp
       trajec(j)%rhoI       = 0.0_dp
       trajec(j)%rhoA       = 0.0_dp
       trajec(j)%rhoNeg     = 0.0_dp
       trajec(j)%DensityN   = 0.0_dp
       trajec(j)%DensityI   = 0.0_dp
       trajec(j)%DensityA   = 0.0_dp
       trajec(j)%DensityNeg = 0.0_dp
       trajec(j)%muN        = 0.0_dp
       trajec(j)%muI        = 0.0_dp
       trajec(j)%muA        = 0.0_dp
       trajec(j)%muNeg      = 0.0_dp
       trajec(j)%Tn         = 0.0_dp
       trajec(j)%Ti         = 0.0_dp
       trajec(j)%Te         = 0.0_dp
       trajec(j)%Vn         = 0.0_dp
       trajec(j)%Vi         = 0.0_dp
       trajec(j)%dVn        = 0.0_dp
       trajec(j)%dVi        = 0.0_dp
       trajec(j)%Vgrad      = 0.0_dp
       trajec(j)%grad_V     = 0.0_dp
       trajec(j)%v_CH       = 0.0_dp
       trajec(j)%v_S        = 0.0_dp
       trajec(j)%v_SH       = 0.0_dp
       trajec(j)%iondeg     = 0.0_dp
       trajec(j)%ionfrac    = 0.0_dp
       trajec(j)%oph2       = 0.0_dp
       trajec(j)%T_gr       = 0.0_dp
       trajec(j)%Teff_gr    = 0.0_dp
       trajec(j)%nlay_gr    = 0.0_dp
       trajec(j)%n_gr       = 0.0_dp
       trajec(j)%mu_gr      = 0.0_dp
       trajec(j)%much_gr    = 0.0_dp
       trajec(j)%r_gr       = 0.0_dp
       trajec(j)%m_gr       = 0.0_dp
       trajec(j)%m_grc      = 0.0_dp
       trajec(j)%m_grm      = 0.0_dp
       trajec(j)%d_ero      = 0.0_dp
       trajec(j)%d_ads      = 0.0_dp
       trajec(j)%rsq_grc    = 0.0_dp
       trajec(j)%rsq_grm    = 0.0_dp
       trajec(j)%compr_n    = 0.0_dp
       trajec(j)%compr_i    = 0.0_dp
      
       ! -------------------------------------------------------------------------
       ! heating and cooling
       ! -------------------------------------------------------------------------
       trajec(j)%line_rad_tot      = 0.0_dp
       trajec(j)%line_rad_n        = 0.0_dp
       trajec(j)%line_rad_i        = 0.0_dp
       trajec(j)%line_rad_e        = 0.0_dp
       trajec(j)%line_rad_n_molec  = 0.0_dp
       trajec(j)%line_rad_n_H2     = 0.0_dp
       trajec(j)%line_rad_n_13CO   = 0.0_dp
       trajec(j)%line_rad_n_OH     = 0.0_dp
       trajec(j)%line_rad_n_NH3    = 0.0_dp
       trajec(j)%line_rad_n_rCO    = 0.0_dp
       trajec(j)%line_rad_n_vCO    = 0.0_dp
       trajec(j)%line_rad_n_CO     = 0.0_dp
       trajec(j)%line_rad_n_roH2O  = 0.0_dp
       trajec(j)%line_rad_n_rpH2O  = 0.0_dp
       trajec(j)%line_rad_n_vH2O   = 0.0_dp
       trajec(j)%line_rad_n_H2O    = 0.0_dp
       trajec(j)%line_rad_i_molec  = 0.0_dp
       trajec(j)%line_rad_i_H2     = 0.0_dp
       trajec(j)%line_rad_e_molec  = 0.0_dp
       trajec(j)%line_rad_e_H2     = 0.0_dp
       trajec(j)%line_rad_e_CO     = 0.0_dp
       trajec(j)%line_rad_n_atoms  = 0.0_dp
       trajec(j)%line_rad_n_Hat    = 0.0_dp
       trajec(j)%line_rad_n_Oat    = 0.0_dp
       trajec(j)%line_rad_n_Op     = 0.0_dp
       trajec(j)%line_rad_n_Cat    = 0.0_dp
       trajec(j)%line_rad_n_Cp     = 0.0_dp
       trajec(j)%line_rad_n_Nat    = 0.0_dp
       trajec(j)%line_rad_n_Np     = 0.0_dp
       trajec(j)%line_rad_n_Sat    = 0.0_dp
       trajec(j)%line_rad_n_Sp     = 0.0_dp
       trajec(j)%line_rad_n_Siat   = 0.0_dp
       trajec(j)%line_rad_n_Sip    = 0.0_dp
       trajec(j)%line_rad_i_atoms  = 0.0_dp
       trajec(j)%line_rad_i_Hat    = 0.0_dp
       trajec(j)%line_rad_i_Oat    = 0.0_dp
       trajec(j)%line_rad_i_Op     = 0.0_dp
       trajec(j)%line_rad_i_Cat    = 0.0_dp
       trajec(j)%line_rad_i_Cp     = 0.0_dp
       trajec(j)%line_rad_i_Nat    = 0.0_dp
       trajec(j)%line_rad_i_Np     = 0.0_dp
       trajec(j)%line_rad_i_Sat    = 0.0_dp
       trajec(j)%line_rad_i_Sp     = 0.0_dp
       trajec(j)%line_rad_i_Siat   = 0.0_dp
       trajec(j)%line_rad_i_Sip    = 0.0_dp
       trajec(j)%line_rad_e_atoms  = 0.0_dp
       trajec(j)%line_rad_e_Hat    = 0.0_dp
       trajec(j)%line_rad_e_Oat    = 0.0_dp
       trajec(j)%line_rad_e_Op     = 0.0_dp
       trajec(j)%line_rad_e_Cat    = 0.0_dp
       trajec(j)%line_rad_e_Cp     = 0.0_dp
       trajec(j)%line_rad_e_Nat    = 0.0_dp
       trajec(j)%line_rad_e_Np     = 0.0_dp
       trajec(j)%line_rad_e_Sat    = 0.0_dp
       trajec(j)%line_rad_e_Sp     = 0.0_dp
       trajec(j)%line_rad_e_Siat   = 0.0_dp
       trajec(j)%line_rad_e_Sip    = 0.0_dp
       trajec(j)%chem_DE_tot       = 0.0_dp
       trajec(j)%chem_DE_n         = 0.0_dp
       trajec(j)%chem_DE_i         = 0.0_dp
       trajec(j)%chem_DE_e         = 0.0_dp
       trajec(j)%chem_DE_n_diss_H2 = 0.0_dp
       trajec(j)%chem_DE_i_diss_H2 = 0.0_dp
       trajec(j)%chem_DE_e_diss_H2 = 0.0_dp
       trajec(j)%chem_DE_n_phgas   = 0.0_dp
       trajec(j)%chem_DE_i_phgas   = 0.0_dp
       trajec(j)%chem_DE_e_phgas   = 0.0_dp
       trajec(j)%chem_DE_n_phgrn   = 0.0_dp
       trajec(j)%chem_DE_i_phgrn   = 0.0_dp
       trajec(j)%chem_DE_e_phgrn   = 0.0_dp
       trajec(j)%chem_DE_n_cosmic  = 0.0_dp
       trajec(j)%chem_DE_i_cosmic  = 0.0_dp
       trajec(j)%chem_DE_e_cosmic  = 0.0_dp
       trajec(j)%chem_DE_n_other   = 0.0_dp
       trajec(j)%chem_DE_i_other   = 0.0_dp
       trajec(j)%chem_DE_e_other   = 0.0_dp
       trajec(j)%elast_scat_tot    = 0.0_dp
       trajec(j)%elast_scat_n      = 0.0_dp
       trajec(j)%elast_scat_i      = 0.0_dp
       trajec(j)%elast_scat_e      = 0.0_dp
       trajec(j)%exch_eint_tot     = 0.0_dp
       trajec(j)%exch_eint_n       = 0.0_dp
       trajec(j)%exch_eint_i       = 0.0_dp
       trajec(j)%exch_eint_e       = 0.0_dp
       trajec(j)%therm_grain_tot   = 0.0_dp
       trajec(j)%mech_trsf_tot     = 0.0_dp
       trajec(j)%mech_trsf_n       = 0.0_dp
       trajec(j)%mech_trsf_i       = 0.0_dp
       trajec(j)%mech_trsf_e       = 0.0_dp
       trajec(j)%mech_trsf_n_chem  = 0.0_dp
       trajec(j)%mech_trsf_i_chem  = 0.0_dp
       trajec(j)%mech_trsf_e_chem  = 0.0_dp
       trajec(j)%mech_trsf_n_compr = 0.0_dp
       trajec(j)%mech_trsf_i_compr = 0.0_dp
       trajec(j)%mech_trsf_e_compr = 0.0_dp
       trajec(j)%mech_trsf_n_visc  = 0.0_dp
       trajec(j)%mech_trsf_i_visc  = 0.0_dp
       trajec(j)%mech_trsf_e_visc  = 0.0_dp

       ! -------------------------------------------------------------------------
       ! energetics
       ! -------------------------------------------------------------------------
       trajec(j)%mom_flux_tot = 0.0_dp
       trajec(j)%mom_flux_kin = 0.0_dp
       trajec(j)%mom_flux_the = 0.0_dp
       trajec(j)%mom_flux_mag = 0.0_dp
       trajec(j)%mom_flux_vis = 0.0_dp
       trajec(j)%nrj_flux_tot = 0.0_dp
       trajec(j)%nrj_flux_kin = 0.0_dp
       trajec(j)%nrj_flux_the = 0.0_dp
       trajec(j)%nrj_flux_mag = 0.0_dp
       trajec(j)%nrj_flux_vis = 0.0_dp
       trajec(j)%nrj_flux_int = 0.0_dp
       trajec(j)%nrj_flux_src = 0.0_dp
       trajec(j)%nrj_flux_pho = 0.0_dp
       trajec(j)%mag_flux_corr= 0.0_dp

       ! -------------------------------------------------------------------------
       ! radiation
       ! -------------------------------------------------------------------------
       trajec(j)%fluph      = 0.0_dp
       trajec(j)%coldens_h  = 0.0_dp
       trajec(j)%coldens_h2 = 0.0_dp
       trajec(j)%coldens_co = 0.0_dp
      
       ! -------------------------------------------------------------------------
       ! chemical profiles
       ! -------------------------------------------------------------------------
       trajec(j)%abon(1:Nspec)           = 0.0_dp
       trajec(j)%coldens(1:Nspec)        = 0.0_dp
       ! To save memory, fill only at the end of computation
       ! trajec(j)%all_rate_chem(1:Nreact) = 0.0_dp
      
       ! -------------------------------------------------------------------------
       ! excitation profiles
       ! -------------------------------------------------------------------------
       trajec(j)%abon_h2lev(1:NH2_lev)       = 0.0_dp
       trajec(j)%cold_h2lev(1:NH2_lev)       = 0.0_dp
       trajec(j)%cdleft_h2lev(1:NH2_lev)     = 0.0_dp
       trajec(j)%emis_h2lin(1:NH2_lines_out) = 0.0_dp
       trajec(j)%inta_h2lin(1:NH2_lines_out) = 0.0_dp
       trajec(j)%abon_colev(1:NCO_lev)       = 0.0_dp
       trajec(j)%cold_colev(1:NCO_lev)       = 0.0_dp
       trajec(j)%cdleft_colev(1:NCO_lev)     = 0.0_dp
       trajec(j)%emis_h  (1:ntrhat )         = 0.0_dp
       trajec(j)%inta_h  (1:ntrhat )         = 0.0_dp
       trajec(j)%emis_c  (1:ntrcat )         = 0.0_dp
       trajec(j)%inta_c  (1:ntrcat )         = 0.0_dp
       trajec(j)%emis_n  (1:ntrnat )         = 0.0_dp
       trajec(j)%inta_n  (1:ntrnat )         = 0.0_dp
       trajec(j)%emis_o  (1:ntroat )         = 0.0_dp
       trajec(j)%inta_o  (1:ntroat )         = 0.0_dp
       trajec(j)%emis_s  (1:ntrsat )         = 0.0_dp
       trajec(j)%inta_s  (1:ntrsat )         = 0.0_dp
       trajec(j)%emis_si (1:ntrsiat)         = 0.0_dp
       trajec(j)%inta_si (1:ntrsiat)         = 0.0_dp
       trajec(j)%emis_cp (1:ntrcpl )         = 0.0_dp
       trajec(j)%inta_cp (1:ntrcpl )         = 0.0_dp
       trajec(j)%emis_np (1:ntrnpl )         = 0.0_dp
       trajec(j)%inta_np (1:ntrnpl )         = 0.0_dp
       trajec(j)%emis_op (1:ntropl )         = 0.0_dp
       trajec(j)%inta_op (1:ntropl )         = 0.0_dp
       trajec(j)%emis_sp (1:ntrspl )         = 0.0_dp
       trajec(j)%inta_sp (1:ntrspl )         = 0.0_dp
       trajec(j)%emis_sip(1:ntrsipl)         = 0.0_dp
       trajec(j)%inta_sip(1:ntrsipl)         = 0.0_dp
    ENDDO

  END SUBROUTINE INITIALIZE_TABLES


  SUBROUTINE DEALLOCATE_TABLE(N, trajec)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     deallocates tables of physical variables
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_LINE_EXCIT

    IMPLICIT NONE

    INTEGER,                                INTENT(in)    :: N
    TYPE (TYPE_TRAJECTORY), DIMENSION(1:N), INTENT(inout) :: trajec
    INTEGER                                               :: i

    DO i = 1, N
       IF ( ALLOCATED(trajec(i)%abon         ) ) DEALLOCATE (trajec(i)%abon         )
       IF ( ALLOCATED(trajec(i)%coldens      ) ) DEALLOCATE (trajec(i)%coldens      )
       IF ( ALLOCATED(trajec(i)%all_rate_chem) ) DEALLOCATE (trajec(i)%all_rate_chem)
       IF ( ALLOCATED(trajec(i)%abon_h2lev   ) ) DEALLOCATE (trajec(i)%abon_h2lev   )
       IF ( ALLOCATED(trajec(i)%cold_h2lev   ) ) DEALLOCATE (trajec(i)%cold_h2lev   )
       IF ( ALLOCATED(trajec(i)%cdleft_h2lev ) ) DEALLOCATE (trajec(i)%cdleft_h2lev )
       IF ( ALLOCATED(trajec(i)%emis_h2lin   ) ) DEALLOCATE (trajec(i)%emis_h2lin   )
       IF ( ALLOCATED(trajec(i)%inta_h2lin   ) ) DEALLOCATE (trajec(i)%inta_h2lin   )
       IF ( ALLOCATED(trajec(i)%abon_colev   ) ) DEALLOCATE (trajec(i)%abon_colev   )
       IF ( ALLOCATED(trajec(i)%cold_colev   ) ) DEALLOCATE (trajec(i)%cold_colev   )
       IF ( ALLOCATED(trajec(i)%cdleft_colev ) ) DEALLOCATE (trajec(i)%cdleft_colev )
       IF ( ALLOCATED(trajec(i)%emis_h       ) ) DEALLOCATE (trajec(i)%emis_h       )
       IF ( ALLOCATED(trajec(i)%inta_h       ) ) DEALLOCATE (trajec(i)%inta_h       )
       IF ( ALLOCATED(trajec(i)%emis_c       ) ) DEALLOCATE (trajec(i)%emis_c       )
       IF ( ALLOCATED(trajec(i)%inta_c       ) ) DEALLOCATE (trajec(i)%inta_c       )
       IF ( ALLOCATED(trajec(i)%emis_n       ) ) DEALLOCATE (trajec(i)%emis_n       )
       IF ( ALLOCATED(trajec(i)%inta_n       ) ) DEALLOCATE (trajec(i)%inta_n       )
       IF ( ALLOCATED(trajec(i)%emis_o       ) ) DEALLOCATE (trajec(i)%emis_o       )
       IF ( ALLOCATED(trajec(i)%inta_o       ) ) DEALLOCATE (trajec(i)%inta_o       )
       IF ( ALLOCATED(trajec(i)%emis_s       ) ) DEALLOCATE (trajec(i)%emis_s       )
       IF ( ALLOCATED(trajec(i)%inta_s       ) ) DEALLOCATE (trajec(i)%inta_s       )
       IF ( ALLOCATED(trajec(i)%emis_si      ) ) DEALLOCATE (trajec(i)%emis_si      )
       IF ( ALLOCATED(trajec(i)%inta_si      ) ) DEALLOCATE (trajec(i)%inta_si      )
       IF ( ALLOCATED(trajec(i)%emis_cp      ) ) DEALLOCATE (trajec(i)%emis_cp      )
       IF ( ALLOCATED(trajec(i)%inta_cp      ) ) DEALLOCATE (trajec(i)%inta_cp      )
       IF ( ALLOCATED(trajec(i)%emis_np      ) ) DEALLOCATE (trajec(i)%emis_np      )
       IF ( ALLOCATED(trajec(i)%inta_np      ) ) DEALLOCATE (trajec(i)%inta_np      )
       IF ( ALLOCATED(trajec(i)%emis_op      ) ) DEALLOCATE (trajec(i)%emis_op      )
       IF ( ALLOCATED(trajec(i)%inta_op      ) ) DEALLOCATE (trajec(i)%inta_op      )
       IF ( ALLOCATED(trajec(i)%emis_sp      ) ) DEALLOCATE (trajec(i)%emis_sp      )
       IF ( ALLOCATED(trajec(i)%inta_sp      ) ) DEALLOCATE (trajec(i)%inta_sp      )
       IF ( ALLOCATED(trajec(i)%emis_sip     ) ) DEALLOCATE (trajec(i)%emis_sip     )
       IF ( ALLOCATED(trajec(i)%inta_sip     ) ) DEALLOCATE (trajec(i)%inta_sip     )
    ENDDO

  END SUBROUTINE DEALLOCATE_TABLE


  SUBROUTINE FILL_PHYS_TABLES(trajec)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     fills tables of physical variables to save the shock profil (trajectory)
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    !     trajec = the point of trajectory to save
    ! results :
    !---------------------------------------------------------------------------

    USE MODULE_CONSTANTS
    USE MODULE_PHYS_VAR
    USE MODULE_GRAINS
    USE MODULE_CHEMICAL_SPECIES
    ! USE MODULE_CHEM_REACT
    USE MODULE_VAR_TH_BALANCE
    USE MODULE_LINE_EXCIT
    USE MODULE_H2
    USE MODULE_CO
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_ENERGETICS
    USE MODULE_RADIATION

    IMPLICIT NONE

    TYPE (TYPE_TRAJECTORY), INTENT(inout) :: trajec
    INTEGER                               :: i

    !---------------------------------------------------------------------------
    ! time, distance, density, velocity, temperature
    !---------------------------------------------------------------------------
    trajec%slab           = slab
    trajec%distance       = distance
    trajec%dist_step      = dist_step
    trajec%Av             = Av
    trajec%timeN          = timeN
    trajec%timeI          = timeI
    trajec%nH             = nH
    trajec%Vsound         = Vsound
    trajec%Vmagnet        = Vmagnet
    trajec%Valfven        = Valfven
    trajec%rhoN           = rhoN
    trajec%rhoI           = rhoI
    trajec%rhoA           = rhoA
    trajec%rhoNeg         = rhoNeg
    trajec%DensityN       = DensityN
    trajec%DensityI       = DensityI
    trajec%DensityA       = DensityA
    trajec%DensityNeg     = DensityNeg
    trajec%muN            = muN
    trajec%muI            = muI
    trajec%muA            = muA
    trajec%muNeg          = muNeg
    trajec%Tn             = Tn
    trajec%Ti             = Ti
    trajec%Te             = Te
    trajec%Vn             = Vn
    trajec%Vi             = Vi
    trajec%dVn            = dVn
    trajec%dVi            = dVi
    trajec%Vgrad          = Vgrad
    trajec%grad_V         = grad_V
    trajec%v_CH           = v_CH
    trajec%v_S            = v_S
    trajec%v_SH           = v_SH
    trajec%iondeg         = DensityI / (DensityN + DensityI + DensityNeg)
    trajec%ionfrac        = DensityI / nH
    trajec%oph2           = op_H2
    trajec%T_gr           = Tgrain
    trajec%Teff_gr        = Teff_grain
    trajec%nlay_gr        = Nlayers
    trajec%n_gr           = dens_grc
    trajec%mu_gr          = Mgrain
    trajec%much_gr        = Mgrain
    trajec%r_gr           = SQRT(rsq_grm)
    trajec%m_gr           = elements(ind_elem_G)%mass
    trajec%m_grc          = Mgrc
    trajec%m_grm          = Mgrm
    trajec%d_ero          = d_ero
    trajec%d_ads          = d_ads
    trajec%rsq_grc        = rsq_grc
    trajec%rsq_grm        = rsq_grm
    trajec%compr_n        = compr_n
    trajec%compr_i        = compr_i

    !---------------------------------------------------------------------------
    ! heating and cooling
    !---------------------------------------------------------------------------
    trajec%line_rad_tot      = line_rad_tot
    trajec%line_rad_n        = line_rad_n
    trajec%line_rad_i        = line_rad_i
    trajec%line_rad_e        = line_rad_e
    trajec%line_rad_n_molec  = line_rad_n_molec
    trajec%line_rad_n_H2     = line_rad_n_H2
    trajec%line_rad_n_13CO   = line_rad_n_13CO
    trajec%line_rad_n_OH     = line_rad_n_OH
    trajec%line_rad_n_NH3    = line_rad_n_NH3
    trajec%line_rad_n_rCO    = line_rad_n_rCO
    trajec%line_rad_n_vCO    = line_rad_n_vCO
    trajec%line_rad_n_CO     = line_rad_n_CO
    trajec%line_rad_n_roH2O  = line_rad_n_roH2O
    trajec%line_rad_n_rpH2O  = line_rad_n_rpH2O
    trajec%line_rad_n_vH2O   = line_rad_n_vH2O
    trajec%line_rad_n_H2O    = line_rad_n_H2O
    trajec%line_rad_i_molec  = line_rad_i_molec
    trajec%line_rad_i_H2     = line_rad_i_H2
    trajec%line_rad_e_molec  = line_rad_e_molec
    trajec%line_rad_e_H2     = line_rad_e_H2
    trajec%line_rad_e_CO     = line_rad_e_CO
    trajec%line_rad_n_atoms  = line_rad_n_atoms
    trajec%line_rad_n_Hat    = line_rad_n_Hat
    trajec%line_rad_n_Oat    = line_rad_n_Oat
    trajec%line_rad_n_Op     = line_rad_n_Op
    trajec%line_rad_n_Cat    = line_rad_n_Cat
    trajec%line_rad_n_Cp     = line_rad_n_Cp
    trajec%line_rad_n_Nat    = line_rad_n_Nat
    trajec%line_rad_n_Np     = line_rad_n_Np
    trajec%line_rad_n_Sat    = line_rad_n_Sat
    trajec%line_rad_n_Sp     = line_rad_n_Sp
    trajec%line_rad_n_Siat   = line_rad_n_Siat
    trajec%line_rad_n_Sip    = line_rad_n_Sip
    trajec%line_rad_i_atoms  = line_rad_i_atoms
    trajec%line_rad_i_Hat    = line_rad_i_Hat
    trajec%line_rad_i_Oat    = line_rad_i_Oat
    trajec%line_rad_i_Op     = line_rad_i_Op
    trajec%line_rad_i_Cat    = line_rad_i_Cat
    trajec%line_rad_i_Cp     = line_rad_i_Cp
    trajec%line_rad_i_Nat    = line_rad_i_Nat
    trajec%line_rad_i_Np     = line_rad_i_Np
    trajec%line_rad_i_Sat    = line_rad_i_Sat
    trajec%line_rad_i_Sp     = line_rad_i_Sp
    trajec%line_rad_i_Siat   = line_rad_i_Siat
    trajec%line_rad_i_Sip    = line_rad_i_Sip
    trajec%line_rad_e_atoms  = line_rad_e_atoms
    trajec%line_rad_e_Hat    = line_rad_e_Hat
    trajec%line_rad_e_Oat    = line_rad_e_Oat
    trajec%line_rad_e_Op     = line_rad_e_Op
    trajec%line_rad_e_Cat    = line_rad_e_Cat
    trajec%line_rad_e_Cp     = line_rad_e_Cp
    trajec%line_rad_e_Nat    = line_rad_e_Nat
    trajec%line_rad_e_Np     = line_rad_e_Np
    trajec%line_rad_e_Sat    = line_rad_e_Sat
    trajec%line_rad_e_Sp     = line_rad_e_Sp
    trajec%line_rad_e_Siat   = line_rad_e_Siat
    trajec%line_rad_e_Sip    = line_rad_e_Sip
    trajec%chem_DE_tot       = chem_DE_tot
    trajec%chem_DE_n         = chem_DE_n
    trajec%chem_DE_i         = chem_DE_i
    trajec%chem_DE_e         = chem_DE_e
    trajec%chem_DE_n_diss_H2 = chem_DE_n_diss_H2
    trajec%chem_DE_i_diss_H2 = chem_DE_i_diss_H2
    trajec%chem_DE_e_diss_H2 = chem_DE_e_diss_H2
    trajec%chem_DE_n_phgas   = chem_DE_n_phgas
    trajec%chem_DE_i_phgas   = chem_DE_i_phgas
    trajec%chem_DE_e_phgas   = chem_DE_e_phgas
    trajec%chem_DE_n_phgrn   = chem_DE_n_phgrn
    trajec%chem_DE_i_phgrn   = chem_DE_i_phgrn
    trajec%chem_DE_e_phgrn   = chem_DE_e_phgrn
    trajec%chem_DE_n_cosmic  = chem_DE_n_cosmic
    trajec%chem_DE_i_cosmic  = chem_DE_i_cosmic
    trajec%chem_DE_e_cosmic  = chem_DE_e_cosmic
    trajec%chem_DE_n_other   = chem_DE_n_other
    trajec%chem_DE_i_other   = chem_DE_i_other
    trajec%chem_DE_e_other   = chem_DE_e_other
    trajec%elast_scat_tot    = elast_scat_tot
    trajec%elast_scat_n      = elast_scat_n
    trajec%elast_scat_i      = elast_scat_i
    trajec%elast_scat_e      = elast_scat_e
    trajec%exch_eint_tot     = exch_eint_tot
    trajec%exch_eint_n       = exch_eint_n
    trajec%exch_eint_i       = exch_eint_i
    trajec%exch_eint_e       = exch_eint_e
    trajec%therm_grain_tot   = therm_grain_tot
    trajec%mech_trsf_tot     = mech_trsf_tot
    trajec%mech_trsf_n       = mech_trsf_n
    trajec%mech_trsf_i       = mech_trsf_i
    trajec%mech_trsf_e       = mech_trsf_e
    trajec%mech_trsf_n_chem  = mech_trsf_n_chem
    trajec%mech_trsf_i_chem  = mech_trsf_i_chem
    trajec%mech_trsf_e_chem  = mech_trsf_e_chem
    trajec%mech_trsf_n_compr = mech_trsf_n_compr
    trajec%mech_trsf_i_compr = mech_trsf_i_compr
    trajec%mech_trsf_e_compr = mech_trsf_e_compr
    trajec%mech_trsf_n_visc  = mech_trsf_n_visc
    trajec%mech_trsf_i_visc  = mech_trsf_i_visc
    trajec%mech_trsf_e_visc  = mech_trsf_e_visc

    !---------------------------------------------------------------------------
    ! energetics
    !---------------------------------------------------------------------------
    trajec%mom_flux_tot   = Momentum_flux
    trajec%mom_flux_kin   = Momentum_flux_kin
    trajec%mom_flux_the   = Momentum_flux_the
    trajec%mom_flux_mag   = Momentum_flux_mag
    trajec%mom_flux_vis   = Momentum_flux_vis
    trajec%nrj_flux_tot   = Energy_flux
    trajec%nrj_flux_kin   = Energy_flux_kin
    trajec%nrj_flux_the   = Energy_flux_the
    trajec%nrj_flux_mag   = Energy_flux_mag
    trajec%nrj_flux_vis   = Energy_flux_vis
    trajec%nrj_flux_int   = Energy_flux_int
    trajec%nrj_flux_src   = Energy_gain
    trajec%nrj_flux_pho   = Energy_flux_pho
    trajec%mag_flux_corr  = mag_flux_corr

    !---------------------------------------------------------------------------
    ! radiation
    !---------------------------------------------------------------------------
    trajec%coldens_h      = coldens_h
    trajec%coldens_h2     = coldens_h2
    trajec%coldens_co     = coldens_co
    trajec%fluph          = fluph

    !---------------------------------------------------------------------------
    ! chemical profiles
    !---------------------------------------------------------------------------
    DO i = 1, Nspec
       trajec%abon(i)     = speci(i)%density
       trajec%coldens(i)  = speci(i)%Col_dens
    ENDDO
    ! To save memory, fill only at the end of computation
    ! DO i = 1, Nreact
    !    trajec%all_rate_chem(i) = react(i)%rate
    ! ENDDO

    !---------------------------------------------------------------------------
    ! excitation profiles
    !---------------------------------------------------------------------------
    DO i = 1, NH2_lev
       trajec%abon_h2lev(i)   = H2_lev(i)%density
       trajec%cold_h2lev(i)   = H2_lev(i)%Col_dens
       trajec%cdleft_h2lev(i) = H2_lev(i)%cd_l
    ENDDO
    DO i = 1, NH2_lines_out
       trajec%emis_h2lin(i) = H2_lines(i)%emiss
       trajec%inta_h2lin(i) = H2_lines(i)%intensity
    ENDDO
    DO i = 1, NCO_lev
       trajec%abon_colev(i)   = CO_lev(i)%density
       trajec%cold_colev(i)   = CO_lev(i)%Col_dens
       trajec%cdleft_colev(i) = CO_lev(i)%cd_l
    ENDDO
    DO i = 1, ntrhat 
       trajec%emis_h(i) = emihat(i)
       trajec%inta_h(i) = inthat(i)
    ENDDO
    DO i = 1, ntrcat 
       trajec%emis_c(i) = emicat(i)
       trajec%inta_c(i) = intcat(i)
    ENDDO
    DO i = 1, ntrnat 
       trajec%emis_n(i) = eminat(i)
       trajec%inta_n(i) = intnat(i)
    ENDDO
    DO i = 1, ntroat 
       trajec%emis_o(i) = emioat(i)
       trajec%inta_o(i) = intoat(i)
    ENDDO
    DO i = 1, ntrsat 
       trajec%emis_s(i) = emisat(i)
       trajec%inta_s(i) = intsat(i)
    ENDDO
    DO i = 1, ntrsiat
       trajec%emis_si(i) = emisiat(i)
       trajec%inta_si(i) = intsiat(i)
    ENDDO
    DO i = 1, ntrcpl 
       trajec%emis_cp(i) = emicpl(i)
       trajec%inta_cp(i) = intcpl(i)
    ENDDO
    DO i = 1, ntrnpl 
       trajec%emis_np(i) = eminpl(i)
       trajec%inta_np(i) = intnpl(i)
    ENDDO
    DO i = 1, ntropl 
       trajec%emis_op(i) = emiopl(i)
       trajec%inta_op(i) = intopl(i)
    ENDDO
    DO i = 1, ntrspl 
       trajec%emis_sp(i) = emispl(i)
       trajec%inta_sp(i) = intspl(i)
    ENDDO
    DO i = 1, ntrsipl
       trajec%emis_sip(i) = emisipl(i)
       trajec%inta_sip(i) = intsipl(i)
    ENDDO

  END SUBROUTINE FILL_PHYS_TABLES


  SUBROUTINE FILL_CHEM_TABLE(N, trajec, Nitermax)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     reconstruct table of chemical rates along the trajectory
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    !     trajec(:)%all_rate_chem
    ! results :
    !---------------------------------------------------------------------------

    USE MODULE_CONSTANTS
    USE MODULE_PHYS_VAR
    USE MODULE_GRAINS
    USE MODULE_H2
    USE MODULE_CO
    USE MODULE_CHEM_REACT
    USE MODULE_DUST_TREATMENT
    USE MODULE_RADIATION
    USE MODULE_EVOLUTION
    USE MODULE_VAR_VODE,      ONLY : integ_type

    IMPLICIT NONE

    INTEGER,                                INTENT(in)    :: N
    INTEGER,                                INTENT(in)    :: Nitermax
    TYPE (TYPE_TRAJECTORY), DIMENSION(1:N), INTENT(inout) :: trajec
    INTEGER                                               :: i
    INTEGER                                               :: j

    DO j = 1, N
       ALLOCATE (trajec(j)%all_rate_chem(1:Nreact)    )
    ENDDO

    !----------------------------------
    ! Reinitialize radiation field
    ! and all opacities
    !----------------------------------
    CALL REINIT_RADIATION

    !----------------------------------
    ! recompute chemical rates
    !----------------------------------
    distance_old = 0.0_dp
    DO j = 1, Nitermax
       CALL REINIT_PHYS_VARIABLES(trajec(j))

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

       IF      (integ_type == 0) THEN
          CALL diffun(d_v_var,distance,v_lvariab,v_dvariab)
       ELSE IF (integ_type == 1) THEN
          CALL diffun(d_v_var,LOG(distance),v_lvariab,v_dvariab)
       ENDIF
       DO i = 1, Nreact
          trajec(j)%all_rate_chem(i) = react(i)%rate
       ENDDO
       distance_old = distance
    ENDDO

  END SUBROUTINE FILL_CHEM_TABLE


  SUBROUTINE REINIT_PHYS_VARIABLES(trajec)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     reinitialize all variables
    ! subroutine/function needed :
    ! input variables :
    !     trajec = the point of trajectory where to reinitialize
    ! output variables :
    ! results :
    !     all variables for reinitialization
    !---------------------------------------------------------------------------

    USE MODULE_CONSTANTS
    USE MODULE_PHYS_VAR
    USE MODULE_GRAINS
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_CHEM_REACT
    USE MODULE_DUST_TREATMENT
    USE MODULE_LINE_EXCIT
    USE MODULE_H2
    USE MODULE_CO
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_ENERGETICS
    USE MODULE_RADIATION
    USE MODULE_EVOLUTION, ONLY : v_dvariab, v_l_var, v_l_der 

    IMPLICIT NONE

    TYPE (TYPE_TRAJECTORY), INTENT(in) :: trajec
    INTEGER                            :: i

    !---------------------------------------------------------------------------
    ! time, distance
    !---------------------------------------------------------------------------
    slab       = trajec%slab
    distance   = trajec%distance
    dist_step  = trajec%dist_step
    Av         = trajec%Av
    timeN      = trajec%timeN
    timeI      = trajec%timeI

    !---------------------------------------------------------------------------
    ! densities
    !---------------------------------------------------------------------------
    nH         = trajec%nH
    rhoN       = trajec%rhoN
    rhoI       = trajec%rhoI
    rhoA       = trajec%rhoA
    rhoNeg     = trajec%rhoNeg
    DensityN   = trajec%DensityN
    DensityI   = trajec%DensityI
    DensityA   = trajec%DensityA
    DensityNeg = trajec%DensityNeg
    muN        = trajec%muN
    muI        = trajec%muI
    muA        = trajec%muA
    muNeg      = trajec%muNeg
    compr_n    = trajec%compr_n
    compr_i    = trajec%compr_i

    !---------------------------------------------------------------------------
    ! velocity
    !---------------------------------------------------------------------------
    Vsound     = trajec%Vsound
    Vmagnet    = trajec%Vmagnet
    Valfven    = trajec%Valfven
    Vn         = trajec%Vn
    Vi         = trajec%Vi
    dVn        = trajec%dVn
    dVi        = trajec%dVi
    Vgrad      = trajec%Vgrad
    grad_V     = trajec%grad_V
    v_CH       = trajec%v_CH
    v_S        = trajec%v_S
    v_SH       = trajec%v_SH

    !---------------------------------------------------------------------------
    ! temperature
    !---------------------------------------------------------------------------
    Tn         = trajec%Tn
    Ti         = trajec%Ti
    Te         = trajec%Te
    Tgrain     = trajec%T_gr
    Teff_grain = trajec%Teff_gr
    op_H2      = trajec%oph2

    !---------------------------------------------------------------------------
    ! grains
    !---------------------------------------------------------------------------
    Nlayers    = trajec%nlay_gr
    dens_grc   = trajec%n_gr
    Mgrain     = trajec%mu_gr
    Mgrain     = trajec%much_gr
    rsq_grm    = trajec%r_gr**2.0_dp
    Mgrc       = trajec%m_grc
    Mgrm       = trajec%m_grm
    d_ero      = trajec%d_ero
    d_ads      = trajec%d_ads
    rsq_grc    = trajec%rsq_grc
    rsq_grm    = trajec%rsq_grm

    !---------------------------------------------------------------------------
    ! energetics
    !---------------------------------------------------------------------------
    Momentum_flux     = trajec%mom_flux_tot
    Momentum_flux_kin = trajec%mom_flux_kin
    Momentum_flux_the = trajec%mom_flux_the
    Momentum_flux_mag = trajec%mom_flux_mag
    Momentum_flux_vis = trajec%mom_flux_vis
    Energy_flux       = trajec%nrj_flux_tot
    Energy_flux_kin   = trajec%nrj_flux_kin
    Energy_flux_the   = trajec%nrj_flux_the
    Energy_flux_mag   = trajec%nrj_flux_mag
    Energy_flux_vis   = trajec%nrj_flux_vis
    Energy_flux_int   = trajec%nrj_flux_int
    Energy_gain       = trajec%nrj_flux_src
    Energy_flux_pho   = trajec%nrj_flux_pho
    mag_flux_corr     = trajec%mag_flux_corr

    !---------------------------------------------------------------------------
    ! chemical profiles
    !---------------------------------------------------------------------------
    DO i = 1, Nspec
       speci(i)%density  = trajec%abon(i)
       speci(i)%Col_dens = trajec%coldens(i)
    ENDDO
    Col_Dens_nH = speci(ind_H)%Col_dens + 2._DP * speci(ind_H2)%Col_dens + speci(ind_Hplus)%Col_Dens
    ! To save memory, fill only at the end of computation
    ! DO i = 1, Nreact
    !    react(i)%rate = trajec%all_rate_chem(i)
    ! ENDDO

    !---------------------------------------------------------------------------
    ! H2 excitation profiles
    !---------------------------------------------------------------------------
    DO i = 1, NH2_lev
       H2_lev(i)%density  = trajec%abon_h2lev(i)
       H2_lev(i)%Col_dens = trajec%cold_h2lev(i)
       H2_lev(i)%cd_l     = trajec%cdleft_h2lev(i)
    ENDDO
    DO i = 1, NH2_lines_out
       H2_lines(i)%emiss     = trajec%emis_h2lin(i)
       H2_lines(i)%intensity = trajec%inta_h2lin(i)
    ENDDO

    !---------------------------------------------------------------------------
    ! CO excitation profiles
    !---------------------------------------------------------------------------
    DO i = 1, NCO_lev
       CO_lev(i)%density  = trajec%abon_colev(i)
       CO_lev(i)%Col_dens = trajec%cold_colev(i)
       CO_lev(i)%cd_l     = trajec%cdleft_colev(i)
    ENDDO

    !---------------------------------------------------------------------------
    ! atomic lines
    !---------------------------------------------------------------------------
    DO i = 1, ntrhat 
       emihat(i) = trajec%emis_h(i)
       inthat(i) = trajec%inta_h(i)
    ENDDO
    DO i = 1, ntrcat 
       emicat(i) = trajec%emis_c(i)
       intcat(i) = trajec%inta_c(i)
    ENDDO
    DO i = 1, ntrnat 
       eminat(i) = trajec%emis_n(i)
       intnat(i) = trajec%inta_n(i)
    ENDDO
    DO i = 1, ntroat 
       emioat(i) = trajec%emis_o(i)
       intoat(i) = trajec%inta_o(i)
    ENDDO
    DO i = 1, ntrsat 
       emisat(i) = trajec%emis_s(i)
       intsat(i) = trajec%inta_s(i)
    ENDDO
    DO i = 1, ntrsiat
       emisiat(i) = trajec%emis_si(i)
       intsiat(i) = trajec%inta_si(i)
    ENDDO
    DO i = 1, ntrcpl 
       emicpl(i) = trajec%emis_cp(i)
       intcpl(i) = trajec%inta_cp(i)
    ENDDO
    DO i = 1, ntrnpl 
       eminpl(i) = trajec%emis_np(i)
       intnpl(i) = trajec%inta_np(i)
    ENDDO
    DO i = 1, ntropl 
       emiopl(i) = trajec%emis_op(i)
       intopl(i) = trajec%inta_op(i)
    ENDDO
    DO i = 1, ntrspl 
       emispl(i) = trajec%emis_sp(i)
       intspl(i) = trajec%inta_sp(i)
    ENDDO
    DO i = 1, ntrsipl
       emisipl(i) = trajec%emis_sip(i)
       intsipl(i) = trajec%inta_sip(i)
    ENDDO

    !---------------------------------------------------------------------------
    ! radiation
    !---------------------------------------------------------------------------
    ! IF (F_AV /= 0) THEN
    !    WRITE(*,*) "Careful - you cannot switch to a CJ type shock and compute properly"
    !    WRITE(*,*) "          the radiative transfer : reinitialization of specific    "
    !    WRITE(*,*) "          intensities as function of angles not done yet !         "
    !    WRITE(*,*) "          -> need approximation to compute radiation based on fluph"
    !    WRITE(*,*) "          -> taking only the grains absorption into account        "
    !    WRITE(*,*) "          -> code stops                                            "
    !    ! The following is absolutely not correct - to be replaced"
    !    ! DO i = 1, nangle
    !    !    irf(i,1:nwlg) = irf_old(i,1:nwlg)
    !    ! ENDDO
    !    STOP
    ! ENDIF
    IF (F_COUP_RAD == 2) THEN
       CALL SET_RAD_GRID(slab)
       CALL FGKCOEF
       DO i = 1, Ncross
          phdest_rate(i)    = SECT_INT(phdest(i))
          phheating_rate(i) = HEATING_PHOTOCHEM_SECT_INT(phdest(i))
       ENDDO
    ENDIF
    fluph = trajec%fluph
    ! Careful, these correspond to the shielding column densities
    coldens_h  = trajec%coldens_h
    coldens_h2 = trajec%coldens_h2
    coldens_co = trajec%coldens_co

    !---------------------------------------------------------------------------
    ! dvode variables and derivatives
    !---------------------------------------------------------------------------
    v_variab(1:d_v_var)     = Zero
    v_lvariab(1:d_v_var)    = Zero
    v_dz_lvariab(1:d_v_var) = Zero
    v_dvariab(1:d_v_var)    = Zero
    v_l_var(1:d_v_var)      = Zero
    v_l_der(1:d_v_var)      = Zero

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

  END SUBROUTINE REINIT_PHYS_VARIABLES


  SUBROUTINE FIND_DISTANCE(z, N, trajec, Nmax, i)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     find the index of a given distance z in the table trajec(:)%distance
    ! subroutine/function needed :
    ! input variables :
    !     distance z
    ! output variables :
    !     index i
    ! results :
    !---------------------------------------------------------------------------

    IMPLICIT none

    REAL(KIND=DP),                          INTENT(in)  :: z
    INTEGER,                                INTENT(in)  :: N
    INTEGER,                                INTENT(in)  :: Nmax
    TYPE (TYPE_TRAJECTORY), DIMENSION(1:N), INTENT(in)  :: trajec
    INTEGER,                                INTENT(out) :: i

    i = 1
    DO WHILE (trajec(i)%distance < z .AND. i < Nmax)
       i = i + 1
    ENDDO
    ! IF (i < N) THEN
    !    i = i - 1
    ! ENDIF
    IF ( trajec(i)%distance >= z ) THEN
       i = i - 1
    ENDIF

  END SUBROUTINE FIND_DISTANCE


  SUBROUTINE COMPUTE_CHOC_SIZE(finstep)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     computes the size of the shock and the corresponding table index
    ! subroutine/function needed :
    ! input variables :
    !     final integration index (counter in mhd_vode.f90)
    ! output variables :
    ! results :
    !     size_shock, time_shock, isize_shock
    !---------------------------------------------------------------------------

    USE MODULE_PHYS_VAR, ONLY : size_shock, time_shock, isize_shock

    IMPLICIT NONE

    INTEGER,      INTENT(in) :: finstep
    REAL(KIND=DP)            :: DE_kin          ! variation of kinetic   energy flux between first and last point
    REAL(KIND=DP)            :: DE_the          ! variation of thermal   energy flux between first and last point
    REAL(KIND=DP)            :: DE_mag          ! variation of magnetic  energy flux between first and last point
    REAL(KIND=DP)            :: DE_vis          ! variation of viscous   energy flux between first and last point
    REAL(KIND=DP)            :: DE_int          ! variation of H2 intern energy flux between first and last point
    REAL(KIND=DP)            :: DE_sum          ! sum of all previous variation  (except for viscous & H2 int nrj, see below)
    REAL(KIND=DP)            :: flux_init_sum   ! sum of all initial energy flux (except for viscous & H2 int nrj, see below)

    REAL(KIND=DP)            :: DE_gain         ! variation of energy gain flux between first and current point
    REAL(KIND=DP), PARAMETER :: thres = 0.99_dp ! threshold of energy variation at which shock is considered as done
    REAL(KIND=DP), PARAMETER :: multsize = 2_dp ! shock size multiplication factor (to include postshock plateau)

    INTEGER                  :: i               ! dummy integer
    INTEGER                  :: ii              ! dummy integer

    DE_kin = traj_main(finstep)%nrj_flux_kin - traj_main(1)%nrj_flux_kin
    DE_the = traj_main(finstep)%nrj_flux_the - traj_main(1)%nrj_flux_the
    DE_mag = traj_main(finstep)%nrj_flux_mag - traj_main(1)%nrj_flux_mag
    DE_vis = traj_main(finstep)%nrj_flux_vis - traj_main(1)%nrj_flux_vis
    DE_int = traj_main(finstep)%nrj_flux_int - traj_main(1)%nrj_flux_int

    !---------------------------------------------------------------------------
    ! Compute sum of variation of all energy fluxes between first and last point
    !
    ! 1 - viscous energy is not included because viscosity is a heating source
    !     term rather than a real energy buffer carried by the gas. We therefore
    !     include viscous flux directly in Bn at the end of diffun (like 
    !     molec_cool, see evolution.f90)
    !
    ! 2 - H2 internal energy is not included because internal energy U is not
    !     included in the energy conservation equation. Only the source term 
    !     H2_energy is taken into account in Bn but dU/dz is not included
    !     anywhere while it should be (See Guillet thesis Eq. A.4 and Draine 
    !     1980).
    !     SHOULD BE MODIFIED - SEE BG'S REMARKS IN EVOLUTION.F90
    !     Proof that treatment is incorrect is that if we add the variation of
    !     internal energy flux, energy is not conserved anymore.
    !---------------------------------------------------------------------------
    DE_sum        = DE_kin + DE_the + DE_mag
    flux_init_sum = traj_main(1)%nrj_flux_kin + traj_main(1)%nrj_flux_the + traj_main(1)%nrj_flux_mag

    !---------------------------------------------------------------------------
    ! Compute shock size
    ! A shock (J or C) is considered as a structure whose main property is to
    ! convert energy : kinetic => magnetic, thermal, visc., internal, radiative
    ! A shock is thus done when the variation of all the energy buffer (kin,
    ! mag, the, int) has effectively been radiated away (the only pure energy
    ! loss of the system).
    ! Parameters : 
    !    - threshold = threshold of nrj variation to cut the shock
    !    - multsize  = shock size mult. factor (to include postshock plateau)
    ! Since internal energy is not correctly treated in the energy conservation
    ! equation, we only consider variation of all the other energy buffer
    !---------------------------------------------------------------------------
    ii = finstep
    DO i = 1, finstep
       DE_gain = traj_main(i)%nrj_flux_src - traj_main(1)%nrj_flux_src
       IF ( (-DE_gain) > thres * (-DE_sum) ) THEN
          ii = i
          EXIT
       ENDIF
    ENDDO
    size_shock = multsize * traj_main(ii)%distance
    isize_shock = finstep
    DO i = 1, finstep
       IF ( traj_main(i)%distance > size_shock ) THEN
          isize_shock = i
          EXIT
       ENDIF
    ENDDO
    size_shock = traj_main(isize_shock)%distance
    time_shock = traj_main(isize_shock)%timeN

  END SUBROUTINE COMPUTE_CHOC_SIZE


END MODULE MODULE_PROFIL_TABLES
