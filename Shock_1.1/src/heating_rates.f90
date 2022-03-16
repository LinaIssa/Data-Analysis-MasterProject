MODULE MODULE_HEATING
  !****************************************************************************
  !** The module 'MODULE_HEATING contains variables used to save the heating **
  !** rates for the output files                                             **
  !****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  !---------------------------------
  ! heating rates (of neutral fluid)
  !---------------------------------
  REAL(KIND=DP) :: heat_n_tot     = 0.0_DP ! total heating term (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_chem    = 0.0_DP ! heating via chemistry (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_diff_in = 0.0_DP ! heating via elastic ion-neutral diffusion (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_diff_an = 0.0_DP ! heating via elastic anion-neutral diffusion (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_diff_en = 0.0_DP ! heating via elastic elect-neutral diffusion (erg cm-3 s-1) 
  REAL(KIND=DP) :: heat_n_diff_gn = 0.0_DP ! heating via elastic grain-neutral diffusion (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_inel_en = 0.0_DP ! heating via inelastic elect-neutral diffusion (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_ther_in = 0.0_DP ! heating via thermal ion-neutral exchange (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_ther_gn = 0.0_DP ! heating via thermal grain-neutral exchange (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_ph      = 0.0_DP ! heating via photoelectric effect (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_cr      = 0.0_DP ! heating via cosmic ray ionizations (erg cm-3 s-1)
  REAL(KIND=DP) :: heat_n_phchem  = 0.0_DP ! total photochemical heating term (erg cm-3 s-1) Tabone 01/18
  REAL(KIND=DP) :: heat_n_compr   = 0.0_DP ! heating via compression
  REAL(KIND=DP) :: heat_n_h2int   = 0.0_DP ! heating via variation of internal energy of H2
  REAL(KIND=DP) :: heat_n_visc    = 0.0_DP ! heating via viscosity (J shock only)

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_HEATING
