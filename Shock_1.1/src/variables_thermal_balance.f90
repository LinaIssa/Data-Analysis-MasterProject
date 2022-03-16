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

MODULE MODULE_VAR_TH_BALANCE
  !****************************************************************************
  !** The module 'MODULE_VAR_TH_BALANCE contains variables to save the 
  !** heating/cooling terms in the output ascii file
  !** heating/cooloing term (erg cm-3 s-1)
  !** Rules for naming the variables : "a_b_c_d"
  !         1) a_b = category of the processus:
  !                   - line_rad    = line radiative processes
  !                   - chem_DE     = energy defect of any reaction (classical or photo)
  !                   - elst_scat   = elastic scattering btw fluids due to drift or thermalization
  !                   - exch_eint   = exchange of internal energy
  !                   - therm_grain = thermalisation of the gas with grains
  !                   - mech_trsf   = mechanical processes
  !         2) _c = fluid: _i, _n, _e, _tot (sum over the three fluids)
  !         3) _d = specific name of the processus
  !****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !------------------------------------------------------------------
  ! A - line radiative processes
  !------------------------------------------------------------------

  REAL(KIND=DP) :: line_rad_tot      = 0.0_DP ! total radiative heating / cooling (sum of the three following terms)
  REAL(KIND=DP) :: line_rad_n        = 0.0_DP ! radiative heating / cooling of the neutral     fluid
  REAL(KIND=DP) :: line_rad_i        = 0.0_DP ! radiative heating / cooling of the pos charged fluid
  REAL(KIND=DP) :: line_rad_e        = 0.0_DP ! radiative heating / cooling of the neg charged fluid

  REAL(KIND=DP) :: line_rad_n_molec  = 0.0_DP ! radiative heating / cooling of the neutral     fluid by all molecules
  REAL(KIND=DP) :: line_rad_n_H2     = 0.0_DP ! radiative heating / cooling of the neutral     fluid by H2
  REAL(KIND=DP) :: line_rad_n_13CO   = 0.0_DP ! radiative heating / cooling of the neutral     fluid by 13CO
  REAL(KIND=DP) :: line_rad_n_OH     = 0.0_DP ! radiative heating / cooling of the neutral     fluid by OH
  REAL(KIND=DP) :: line_rad_n_NH3    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by NH3
  REAL(KIND=DP) :: line_rad_n_rCO    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by rCO
  REAL(KIND=DP) :: line_rad_n_vCO    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by vCO
  REAL(KIND=DP) :: line_rad_n_CO     = 0.0_DP ! radiative heating / cooling of the neutral     fluid by CO
  REAL(KIND=DP) :: line_rad_n_roH2O  = 0.0_DP ! radiative heating / cooling of the neutral     fluid by roH2O
  REAL(KIND=DP) :: line_rad_n_rpH2O  = 0.0_DP ! radiative heating / cooling of the neutral     fluid by rpH2O
  REAL(KIND=DP) :: line_rad_n_vH2O   = 0.0_DP ! radiative heating / cooling of the neutral     fluid by vH2O
  REAL(KIND=DP) :: line_rad_n_H2O    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by H2O

  REAL(KIND=DP) :: line_rad_i_molec  = 0.0_DP ! radiative heating / cooling of the pos charged fluid by all molecules
  REAL(KIND=DP) :: line_rad_i_H2     = 0.0_DP ! radiative heating / cooling of the pos charged fluid by H2

  REAL(KIND=DP) :: line_rad_e_molec  = 0.0_DP ! radiative heating / cooling of the neg charged fluid by all molecules
  REAL(KIND=DP) :: line_rad_e_H2     = 0.0_DP ! radiative heating / cooling of the neg charged fluid by H2
  REAL(KIND=DP) :: line_rad_e_CO     = 0.0_DP ! radiative heating / cooling of the neg charged fluid by CO

  REAL(KIND=DP) :: line_rad_n_atoms  = 0.0_DP ! radiative heating / cooling of the neutral     fluid by all atoms
  REAL(KIND=DP) :: line_rad_n_Hat    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by H
  REAL(KIND=DP) :: line_rad_n_Oat    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by O
  REAL(KIND=DP) :: line_rad_n_Op     = 0.0_DP ! radiative heating / cooling of the neutral     fluid by O+
  REAL(KIND=DP) :: line_rad_n_Cat    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by C
  REAL(KIND=DP) :: line_rad_n_Cp     = 0.0_DP ! radiative heating / cooling of the neutral     fluid by C+
  REAL(KIND=DP) :: line_rad_n_Nat    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by N
  REAL(KIND=DP) :: line_rad_n_Np     = 0.0_DP ! radiative heating / cooling of the neutral     fluid by N+
  REAL(KIND=DP) :: line_rad_n_Sat    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by S
  REAL(KIND=DP) :: line_rad_n_Sp     = 0.0_DP ! radiative heating / cooling of the neutral     fluid by S+
  REAL(KIND=DP) :: line_rad_n_Siat   = 0.0_DP ! radiative heating / cooling of the neutral     fluid by Si
  REAL(KIND=DP) :: line_rad_n_Sip    = 0.0_DP ! radiative heating / cooling of the neutral     fluid by Si+

  REAL(KIND=DP) :: line_rad_i_atoms  = 0.0_DP ! radiative heating / cooling of the pos charged fluid by all atoms
  REAL(KIND=DP) :: line_rad_i_Hat    = 0.0_DP ! radiative heating / cooling of the pos charged fluid by H
  REAL(KIND=DP) :: line_rad_i_Oat    = 0.0_DP ! radiative heating / cooling of the pos charged fluid by O
  REAL(KIND=DP) :: line_rad_i_Op     = 0.0_DP ! radiative heating / cooling of the pos charged fluid by O+
  REAL(KIND=DP) :: line_rad_i_Cat    = 0.0_DP ! radiative heating / cooling of the pos charged fluid by C
  REAL(KIND=DP) :: line_rad_i_Cp     = 0.0_DP ! radiative heating / cooling of the pos charged fluid by C+
  REAL(KIND=DP) :: line_rad_i_Nat    = 0.0_DP ! radiative heating / cooling of the pos charged fluid by N
  REAL(KIND=DP) :: line_rad_i_Np     = 0.0_DP ! radiative heating / cooling of the pos charged fluid by N+
  REAL(KIND=DP) :: line_rad_i_Sat    = 0.0_DP ! radiative heating / cooling of the pos charged fluid by S
  REAL(KIND=DP) :: line_rad_i_Sp     = 0.0_DP ! radiative heating / cooling of the pos charged fluid by S+
  REAL(KIND=DP) :: line_rad_i_Siat   = 0.0_DP ! radiative heating / cooling of the pos charged fluid by Si
  REAL(KIND=DP) :: line_rad_i_Sip    = 0.0_DP ! radiative heating / cooling of the pos charged fluid by Si+

  REAL(KIND=DP) :: line_rad_e_atoms  = 0.0_DP ! radiative heating / cooling of the neg charged fluid by all atoms
  REAL(KIND=DP) :: line_rad_e_Hat    = 0.0_DP ! radiative heating / cooling of the neg charged fluid by H
  REAL(KIND=DP) :: line_rad_e_Oat    = 0.0_DP ! radiative heating / cooling of the neg charged fluid by O
  REAL(KIND=DP) :: line_rad_e_Op     = 0.0_DP ! radiative heating / cooling of the neg charged fluid by O+
  REAL(KIND=DP) :: line_rad_e_Cat    = 0.0_DP ! radiative heating / cooling of the neg charged fluid by C
  REAL(KIND=DP) :: line_rad_e_Cp     = 0.0_DP ! radiative heating / cooling of the neg charged fluid by C+
  REAL(KIND=DP) :: line_rad_e_Nat    = 0.0_DP ! radiative heating / cooling of the neg charged fluid by N
  REAL(KIND=DP) :: line_rad_e_Np     = 0.0_DP ! radiative heating / cooling of the neg charged fluid by N+
  REAL(KIND=DP) :: line_rad_e_Sat    = 0.0_DP ! radiative heating / cooling of the neg charged fluid by S
  REAL(KIND=DP) :: line_rad_e_Sp     = 0.0_DP ! radiative heating / cooling of the neg charged fluid by S+
  REAL(KIND=DP) :: line_rad_e_Siat   = 0.0_DP ! radiative heating / cooling of the neg charged fluid by Si
  REAL(KIND=DP) :: line_rad_e_Sip    = 0.0_DP ! radiative heating / cooling of the neg charged fluid by Si+


  !------------------------------------------------------------------
  ! B - energy defect of any reaction (including photoreactions)
  !------------------------------------------------------------------

  REAL(KIND=DP) :: chem_DE_tot       = 0.0_DP ! total chemical heating / cooling (sum of the three following terms)
  REAL(KIND=DP) :: chem_DE_n         = 0.0_DP ! chemical heating / cooling of the neutral     fluid
  REAL(KIND=DP) :: chem_DE_i         = 0.0_DP ! chemical heating / cooling of the pos charged fluid
  REAL(KIND=DP) :: chem_DE_e         = 0.0_DP ! chemical heating / cooling of the neg charged fluid

  REAL(KIND=DP) :: chem_DE_n_diss_H2 = 0.0_DP ! chemical heating / cooling of the neutral     fluid by H2 dissociation
  REAL(KIND=DP) :: chem_DE_i_diss_H2 = 0.0_DP ! chemical heating / cooling of the pos charged fluid by H2 dissociation
  REAL(KIND=DP) :: chem_DE_e_diss_H2 = 0.0_DP ! chemical heating / cooling of the neg charged fluid by H2 dissociation

  REAL(KIND=DP) :: chem_DE_n_phgas   = 0.0_DP ! chemical heating / cooling of the neutral     fluid by photoreactions
  REAL(KIND=DP) :: chem_DE_i_phgas   = 0.0_DP ! chemical heating / cooling of the pos charged fluid by photoreactions
  REAL(KIND=DP) :: chem_DE_e_phgas   = 0.0_DP ! chemical heating / cooling of the neg charged fluid by photoreactions

  REAL(KIND=DP) :: chem_DE_n_phgrn   = 0.0_DP ! chemical heating / cooling of the neutral     fluid by photoelectric effect (grains)
  REAL(KIND=DP) :: chem_DE_i_phgrn   = 0.0_DP ! chemical heating / cooling of the pos charged fluid by photoelectric effect (grains)
  REAL(KIND=DP) :: chem_DE_e_phgrn   = 0.0_DP ! chemical heating / cooling of the neg charged fluid by photoelectric effect (grains)

  REAL(KIND=DP) :: chem_DE_n_cosmic  = 0.0_DP ! chemical heating / cooling of the neutral     fluid by cosmic ray ionization
  REAL(KIND=DP) :: chem_DE_i_cosmic  = 0.0_DP ! chemical heating / cooling of the pos charged fluid by cosmic ray ionization
  REAL(KIND=DP) :: chem_DE_e_cosmic  = 0.0_DP ! chemical heating / cooling of the neg charged fluid by cosmic ray ionization

  REAL(KIND=DP) :: chem_DE_n_other   = 0.0_DP ! chemical heating / cooling of the neutral     fluid by all other chemical process
  REAL(KIND=DP) :: chem_DE_i_other   = 0.0_DP ! chemical heating / cooling of the pos charged fluid by all other chemical process
  REAL(KIND=DP) :: chem_DE_e_other   = 0.0_DP ! chemical heating / cooling of the neg charged fluid by all other chemical process


  !------------------------------------------------------------------
  ! C - elastic scattering btw fluids due to drift or thermalization
  !------------------------------------------------------------------

  REAL(KIND=DP) :: elast_scat_tot    = 0.0_DP ! total elastic scat heating / cooling (sum of the three following terms)
  REAL(KIND=DP) :: elast_scat_n      = 0.0_DP ! elastic scat heating / cooling of the neutral     fluid
  REAL(KIND=DP) :: elast_scat_i      = 0.0_DP ! elastic scat heating / cooling of the pos charged fluid
  REAL(KIND=DP) :: elast_scat_e      = 0.0_DP ! elastic scat heating / cooling of the neg charged fluid


  !------------------------------------------------------------------
  ! D - exchange of internal energy between fluids
  !------------------------------------------------------------------

  REAL(KIND=DP) :: exch_eint_tot     = 0.0_DP ! total heat exchange (sum of the three following terms)
  REAL(KIND=DP) :: exch_eint_n       = 0.0_DP ! heat exchange heating / cooling of the neutral     fluid
  REAL(KIND=DP) :: exch_eint_i       = 0.0_DP ! heat exchange heating / cooling of the pos charged fluid
  REAL(KIND=DP) :: exch_eint_e       = 0.0_DP ! heat exchange heating / cooling of the neg charged fluid


  !------------------------------------------------------------------
  ! E - thermalization of the gas with grains
  !------------------------------------------------------------------

  REAL(KIND=DP) :: therm_grain_tot   = 0.0_DP ! total heating / cooling by thermalization with grains

  !------------------------------------------------------------------
  ! F - mechanical processes
  !------------------------------------------------------------------
  REAL(KIND=DP) :: mech_trsf_tot     = 0.0_DP ! total mechanical heating / cooling (sum of the three following terms)
  REAL(KIND=DP) :: mech_trsf_n       = 0.0_DP ! mechanical heating / cooling of the neutral     fluid
  REAL(KIND=DP) :: mech_trsf_i       = 0.0_DP ! mechanical heating / cooling of the pos charged fluid
  REAL(KIND=DP) :: mech_trsf_e       = 0.0_DP ! mechanical heating / cooling of the neg charged fluid

  REAL(KIND=DP) :: mech_trsf_n_chem  = 0.0_DP ! mechanical heating / cooling of the neutral     fluid by chemical modification of kin / int energies
  REAL(KIND=DP) :: mech_trsf_i_chem  = 0.0_DP ! mechanical heating / cooling of the pos charged fluid by chemical modification of kin / int energies
  REAL(KIND=DP) :: mech_trsf_e_chem  = 0.0_DP ! mechanical heating / cooling of the neg charged fluid by chemical modification of kin / int energies

  REAL(KIND=DP) :: mech_trsf_n_compr = 0.0_DP ! mechanical heating / cooling of the neutral     fluid by compression (net contribution)
  REAL(KIND=DP) :: mech_trsf_i_compr = 0.0_DP ! mechanical heating / cooling of the pos charged fluid by compression (net contribution)
  REAL(KIND=DP) :: mech_trsf_e_compr = 0.0_DP ! mechanical heating / cooling of the neg charged fluid by compression (net contribution)
    
  REAL(KIND=DP) :: mech_trsf_n_visc  = 0.0_DP ! mechanical heating / cooling of the neutral     fluid by viscous dissipation
  REAL(KIND=DP) :: mech_trsf_i_visc  = 0.0_DP ! mechanical heating / cooling of the pos charged fluid by viscous dissipation
  REAL(KIND=DP) :: mech_trsf_e_visc  = 0.0_DP ! mechanical heating / cooling of the neg charged fluid by viscous dissipation

  !------------------------------------------------------------------
  ! G - H2 cooling analysis
  !------------------------------------------------------------------
  REAL(KIND=DP) :: h2excit_n_coh = 0.0_DP
  REAL(KIND=DP) :: h2excit_n_add = 0.0_DP
  REAL(KIND=DP) :: h2excit_i_add = 0.0_DP
  REAL(KIND=DP) :: h2excit_e_add = 0.0_DP
  REAL(KIND=DP) :: h2excit_e_LW  = 0.0_DP

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_VAR_TH_BALANCE
