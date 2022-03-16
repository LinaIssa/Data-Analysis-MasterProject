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

MODULE MODULE_OUTPUTS
  !*****************************************************************************
  !** The module 'MODULE_OUTPUTS' contains output subroutines (screen, files) **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


CONTAINS


  SUBROUTINE WRITE_INFO
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     write informations about shock parameters, H2 levels, grains,
    !     chemical species and reactions in the file file_info
    ! subroutine/function needed :
    !     INFO_SPECY
    !     INFO_REACTION
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_TECHCONFIG
    USE MODULE_CHEM_REACT
    USE MODULE_PHYS_VAR
    USE MODULE_H2
    USE MODULE_VAR_VODE, ONLY : duration_max, Eps_V, T0_V, H0_V, &
                                Tout_V, MF_V, Itask_V, Istate_V, Hdone, Hnext
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_GRAINS
    USE MODULE_ENERGETICS
    USE MODULE_TOOLS, ONLY    : GET_FILE_NUMBER

    IMPLICIT NONE

    INTEGER(KIND=LONG), DIMENSION(8)         :: date
    INTEGER(KIND=LONG)                       :: i
    CHARACTER(LEN=lenfilename)               :: name_file_info = 'info_mhd.out'
    INTEGER                                  :: file_info

    ! file opening
    file_info = GET_FILE_NUMBER()
    name_file_info = TRIM(out_dir) // TRIM(name_file_info)
    OPEN(file_info,file=name_file_info,status='REPLACE',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    ! date, time
    CALL DATE_AND_TIME(values=date)
    WRITE(file_info,'("date : ",I2,"/",I2,"/",I4)',ADVANCE='NO') date(3),date(2),date(1)
    WRITE(file_info,'(5X,"time : ",I2," H ",I2," min")')         date(5),date(6)

    ! --- input parameters ---
    WRITE(file_info,*)
    WRITE(file_info,'("input files")')
    WRITE(file_info,'(T11,"modele           :",A)')    modele  
    WRITE(file_info,'(T11,"specfile         :",A)')    specfile
    WRITE(file_info,'(T11,"chemfile         :",A)')    chemfile
    WRITE(file_info,'(T11,"h2exfile         :",A)')    h2exfile
    WRITE(file_info,'(T11,"gridfile         :",A)')    gridfile
    WRITE(file_info,*)
    WRITE(file_info,'("shock parameters")')
    WRITE(file_info,'(T11,"shock type       :",A)')     shock_type_init
    WRITE(file_info,'(T11,"Nfluids          :",I9)')    Nfluids
    WRITE(file_info,'(T11,"B (micro Gauss)  :",ES9.2)') Bfield*1.D6
    WRITE(file_info,'(T11,"Vs (km.s-1)      :",F9.1)')  Vs_km
    WRITE(file_info,'(T11,"Dv (km.s-1)      :",F9.1)')  DeltaVmin
    WRITE(file_info,'(T11,"nH (cm-3)        :",ES9.2)') nH
    WRITE(file_info,'(T11,"Tn (K)           :",F9.1)')  Tn
    WRITE(file_info,'(T11,"o/p H2           :",ES9.2)') op_H2
    WRITE(file_info,*)
    WRITE(file_info,'("environmental parameters")')
    WRITE(file_info,'(T11,"Zeta (s-1)       :",ES9.2)') Zeta
    WRITE(file_info,'(T11,"ISRF             :",I9)')    F_ISRF
    WRITE(file_info,'(T11,"RAD              :",F9.2)')  RAD
    WRITE(file_info,'(T11,"Av (mag)         :",F9.2)')  Av0
    WRITE(file_info,'(T11,"F_COUP_RAD       :",I9)')    F_COUP_RAD
    WRITE(file_info,'(T11,"F_AV             :",I9)')    F_AV      
    WRITE(file_info,'(T11,"F_invAv          :",I9)')    F_invAv   
    WRITE(file_info,'(T11,"Av conv factor   :",F9.2)')  inv_Av_fac
    WRITE(file_info,'(T11,"N_H2_0 (cm-2)    :",ES9.3)') N_H2_0
    WRITE(file_info,'(T11,"N_CO_0 (cm-2)    :",ES9.3)') N_CO_0
    WRITE(file_info,'(T11,"vturb            :",ES9.3)') vturb
    WRITE(file_info,*)
    WRITE(file_info,'("grains parameters")')
    WRITE(file_info,'(T11,"F_TGR            :",I9)')    F_TGR
    WRITE(file_info,'(T11,"Tgrains          :",ES9.2)') Tgrain 
    WRITE(file_info,'(T11,"amin_mrn         :",ES9.2)') amin_mrn
    WRITE(file_info,'(T11,"amax_mrn         :",ES9.2)') amax_mrn
    WRITE(file_info,'(T11,"alph_mrn         :",ES9.2)') alph_mrn
    WRITE(file_info,'(T11,"rho_grc          :",ES9.2)') rho_grc 
    WRITE(file_info,'(T11,"rho_grm          :",ES9.2)') rho_grm 
    WRITE(file_info,*)
    WRITE(file_info,'("excitation & cooling")')
    WRITE(file_info,'(T11,"ieqth            :",I9)')    ieqth
    WRITE(file_info,'(T11,"Cool_KN          :",I9)')    Cool_KN
    WRITE(file_info,'(T11,"NH2_lev          :",I9)')    NH2_lev
    WRITE(file_info,'(T11,"NH2_lines        :",I9)')    NH2_lines_out
    WRITE(file_info,'(T11,"H_H2_flag        :",A)')     H_H2_flag
    WRITE(file_info,'(T11,"iforH2           :",I9)')    iforH2
    WRITE(file_info,'(T11,"ikinH2           :",I9)')    ikinH2
    WRITE(file_info,'(T11,"pumpH2           :",I9)')    pumpH2
    WRITE(file_info,'(T11,"NCO_lev          :",I9)')    NCO_lev
    WRITE(file_info,*)
    WRITE(file_info,'("integration parameters")')
    WRITE(file_info,'(T11,"Nstep_max        :",I9)')    Nstep_max
    WRITE(file_info,'(T11,"timeJ            :",ES9.2)') timeJ
    WRITE(file_info,'(T11,"duration_max     :",ES9.2)') duration_max
    WRITE(file_info,'(T11,"XLL              :",ES9.2)') XLL
    WRITE(file_info,*)
    WRITE(file_info,'("DVODE parameters ")')
    WRITE(file_info,'(T11,"Eps_V            :",ES9.2)') Eps_V
    WRITE(file_info,'(T11,"T0_V             :",ES9.2)') T0_V
    WRITE(file_info,'(T11,"H0_V             :",ES9.2)') H0_V
    WRITE(file_info,'(T11,"Tout_V           :",ES9.2)') Tout_V
    WRITE(file_info,'(T11,"MF_V             :",I9)')    MF_V
    WRITE(file_info,'(T11,"Itask_V          :",I9)')    Itask_V
    WRITE(file_info,'(T11,"Istate_V         :",I9)')    Istate_V
    WRITE(file_info,'(T11,"Hdone            :",ES9.2)') Hdone
    WRITE(file_info,'(T11,"Hnext            :",ES9.2)') Hnext
    WRITE(file_info,*)
    WRITE(file_info,'("outputs parameters")')
    WRITE(file_info,'(T11,"F_W_HDF5_STD     :",I9)')    F_W_HDF5_STD
    WRITE(file_info,'(T11,"F_W_HDF5_CHE     :",I9)')    F_W_HDF5_CHE
    WRITE(file_info,'(T11,"F_W_ASCII        :",I9)')    F_W_ASCII
    WRITE(file_info,'(T11,"Npthdf5          :",I9)')    Npthdf5
    WRITE(file_info,'(T11,"Nstep_w          :",I9)')    Nstep_w
    WRITE(file_info,'(T11,"speci_out        :",A)')     speci_out
    WRITE(file_info,'(T11,"H2_out           :",A)')     H2_out
    WRITE(file_info,'(T11,"line_out         :",A)')     line_out
    WRITE(file_info,'(T11,"flag_analysis    :",A)')     flag_analysis

    ! --- H2 levels ---
    WRITE(file_info,*)
    WRITE(file_info,'("H2 molecule ")')
    WRITE(file_info,'(T11,"number of levels : ",I7)') NH2_lev
    WRITE(file_info,'(T11,"E(v=",I2,",J=",I2,") (K) : ",F7.1,"  (last level)")') &
         H2_lev(NH2_lev)%V, &
         H2_lev(NH2_lev)%J, &
         H2_lev(NH2_lev)%energy
    WRITE(file_info,'(T11,"Vmax             : ", I7)') Vmax_H2
    WRITE(file_info,'(T11,"Jmax             : ", I7)') Jmax_H2
    WRITE(file_info,'(T11,"H-H2 collisions  : ", A4)') H_H2_flag
    WRITE(file_info,'(T11,"iforH2           : ", I7)') iforH2
    WRITE(file_info,'(T11,"ikinH2           : ", I7)') ikinH2
    WRITE(file_info,'(T11,"H2_int_E         : ", ES9.2," K")') H2_int_E

    ! --- elemental abundances ---
    WRITE(file_info,*)
    WRITE(file_info,'("elemental abundances (gas + mantles + PAH)")')
    DO i=1,Nelements
       WRITE(file_info,'(T11,A," : ",ES8.2, "   (ref : ",ES8.2,")")') &
            elements(i)%name, DBLE(elements(i)%ab_init), DBLE(elements(i)%ab_ref)
    END DO

    CALL info_elements(file_info)

    ! --- energetics ---
    WRITE(file_info,*)
    WRITE(file_info,'("energetics")')
    WRITE(file_info,'(T11,"mass flux (g/s/cm2)     :",ES9.2)') Mass_flux_init
    WRITE(file_info,'(T11,"momentum flux (erg/cm3) :",ES9.2)') Momentum_flux_init
    WRITE(file_info,'(T11,"energy flux (erg/s/cm2) :",ES9.2)') Energy_flux_init

    ! --- chemical species ---
    WRITE(file_info,*)
    WRITE(file_info,'(80("-"))')
    WRITE(file_info,'("-- ",I3," chemical species (+",I1," added)",T78," --")') Nspec, Nspec_plus-Nspec

    WRITE(file_info,'("-- incl. ",I2," neutrals, ",I2," on grains, &
         &", I2, " ions >0, ", I2, " ions <0",T78," --")') &
         Nneutrals, Nongrains, Nions, Nneg
    WRITE(file_info,'("-- name, composition, enthalpy(kCal/mol, -99.999=unknown)&
         &, Density (cm-3)",T78," --")')
    WRITE(file_info,'(80("-"))')
    DO i=1, Nspec
       CALL WRITE_SPECY(file_info,speci(i))
    ENDDO
    DO i=Nspec+1,Nspec_plus
       WRITE(file_info,'(I3,2X,A7)')speci(i)%index, speci(i)%name
    END DO

    ! --- chemical reactions ---
    WRITE(file_info,*)
    WRITE(file_info,'(80("-"))')
    WRITE(file_info,'("--- ",I4," chemical reactions",T77," ---")')Nreact
    WRITE(file_info,'("--- incl. ",I3," PHOTO, ",I3," CR_IO, ", &
         &I3," CR_DE, ",I3," H2_FO, ",I3," THREE, ",I3," SPUTT",T77," ---")') &
         &Nphoto, Ncr_io, Ncr_de, Nh2_fo, Nthree, Nsputt
    WRITE(file_info,'("---       ",I3," EROSI, ",I3," ADSOR, ",I3,&
         &" DISSO, ",I3," OTHER, and ",I3," REVER",T77," ---")') &
         &Nerosi, Nadsor, Ndisso, Nother, Nrever
    WRITE(file_info,'("--- R1 + R2 = P1 + P2 + P3 + P4, gamma (cm3.s-1), &
         &alpha, beta (K), DE (eV)",T77," ---")')
    WRITE(file_info,'(80("-"))')
    DO i=1, Nreact
       CALL WRITE_REACTION(file_info,react(i))
    END DO

    ! file closing
    CLOSE(file_info)

  END SUBROUTINE WRITE_INFO


  SUBROUTINE WRITE_INFO_END
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     add an ending sentence in info_mhd.out if the shock model run properly
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_TECHCONFIG
    USE MODULE_TOOLS, ONLY    : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR, ONLY : size_shock, time_shock, shock_type_fin

    IMPLICIT NONE

    CHARACTER(LEN=lenfilename)               :: name_file_info = 'info_mhd.out'
    INTEGER                                  :: file_info

    file_info = GET_FILE_NUMBER()
    name_file_info = TRIM(out_dir) // TRIM(name_file_info)
    OPEN(file_info,file=name_file_info,status='old',action='WRITE',position='append')
    WRITE(file_info,'(A25,ES18.9E3)') 'shock size (cm) ........ ', size_shock
    WRITE(file_info,'(A25,ES18.9E3)') 'shock timescale (yr) ... ', time_shock
    WRITE(file_info,'(A25,A2)')       'final shock type ....... ', shock_type_fin
    WRITE(file_info,'(A17)') 'END_OF_MODEL...ok'
    CLOSE(file_info)
  END SUBROUTINE WRITE_INFO_END


  SUBROUTINE WRITE_SCREEN
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     Write informations (counter, distance, time, temperatures) on screen
    !     during integration.
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CONSTANTS, ONLY : screen
    USE MODULE_PHYS_VAR
    USE MODULE_VAR_VODE

    IMPLICIT NONE

    WRITE(screen,'("step= ",I7, 2X)', ADVANCE='NO') counter
    WRITE(screen,'("z=", ES10.3, 2X)', ADVANCE='NO') distance
    WRITE(screen,'("Hdone=", ES9.3, 2X)', ADVANCE='NO') Hdone
    WRITE(screen,'("TimeN=", ES9.3, 2X)', ADVANCE='NO') timeN
    WRITE(screen,'("Av=", ES9.3, 2X)', ADVANCE='NO') Av
    ! WRITE(screen,'("NH =", ES9.3, 2X)', ADVANCE='NO') coldens_h
    ! WRITE(screen,'("NH2=", ES9.3, 2X)', ADVANCE='NO') coldens_h2
    ! WRITE(screen,'("NCO=", ES9.3, 2X)', ADVANCE='NO') coldens_co
    ! WRITE(screen,'("N_H=", ES9.3, 2X)', ADVANCE='NO') Col_Dens_nH
    WRITE(screen,'("Tn=", ES9.3, 2X)', ADVANCE='NO') Tn
    ! WRITE(screen,'("Ti=", ES9.3, 2X)', ADVANCE='NO') Ti
    WRITE(screen,'("dV=", ES10.3, 2X)', ADVANCE='YES') DeltaV

  END SUBROUTINE WRITE_SCREEN


  SUBROUTINE WRITE_OUTPUT(state)
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     write several variables in corresponding output files. 
    !     - if called with "first", opens the files and write header
    !     - if called with "blank", writes two blank lines in files
    !     - if called with "final", closes the files.
    !     - otherwise, writes results.
    ! subroutine/function needed :
    !     GET_FILE_NUMBER
    ! input variables :
    !     state (character) -> what should we do
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS,            ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR
    USE MODULE_GRAINS
    USE MODULE_H2
    USE MODULE_CHEMICAL_SPECIES, ONLY : Nspec, speci, ind_H2, ind_H, ind_Hplus &
                                      , ind_PAH0, ind_PAHp, ind_PAHm, ind_elem_H &
                                      , Nelements, elements
    USE MODULE_VAR_TH_BALANCE
    USE MODULE_MOLECULAR_COOLING
    USE MODULE_ENERGETICS
    USE MODULE_LINE_EXCIT
    USE MODULE_RADIATION
    IMPLICIT NONE

    CHARACTER(len=5), INTENT(in) :: state
    INTEGER(KIND=LONG)           :: i
    INTEGER(KIND=LONG)           :: ii
    CHARACTER(len=15)            :: str,str2
    CHARACTER(LEN=*), PARAMETER  :: format_header='(A18,1X)'
    CHARACTER(LEN=*), PARAMETER  :: format_lng='(ES18.9E3,1X)'
    CHARACTER(LEN=lenfilename)   :: name_file_phys        = 'mhd_phys.out'
    INTEGER(KIND=LONG), SAVE     :: file_phys
    CHARACTER(LEN=lenfilename)   :: name_file_speci       = 'mhd_speci.out'
    INTEGER(KIND=LONG), SAVE     :: file_speci
    CHARACTER(LEN=lenfilename)   :: name_file_coldens     = 'mhd_coldens.out'
    INTEGER(KIND=LONG), SAVE     :: file_coldens
    CHARACTER(LEN=lenfilename)   :: name_file_H2_lev      = 'H2_lev.out'
    INTEGER(KIND=LONG), SAVE     :: file_H2_lev
    CHARACTER(LEN=lenfilename)   :: name_file_H2_line     = 'H2_line.out'
    INTEGER(KIND=LONG), SAVE     :: file_H2_line
    CHARACTER(LEN=lenfilename)   :: name_file_therm_bal   = 'thermal_balance.out'
    INTEGER(KIND=LONG), SAVE     :: file_therm_bal
    CHARACTER(LEN=lenfilename)   :: name_file_energetics  = 'energetics.out'
    INTEGER(KIND=LONG), SAVE     :: file_energetics
    CHARACTER(LEN=lenfilename)   :: name_file_intensity   = 'intensity.out'
    INTEGER(KIND=LONG), SAVE     :: file_intensity
    CHARACTER(LEN=lenfilename)   :: name_file_populations = 'populations.out'
    INTEGER(KIND=LONG), SAVE     :: file_populations
    CHARACTER(LEN=lenfilename)   :: name_file_fe_pops     = 'fe_pops.out'
    INTEGER(KIND=LONG), SAVE     :: file_fe_pops
    CHARACTER(LEN=lenfilename)   :: name_file_fe_lines    = 'fe_lines.out'
    INTEGER(KIND=LONG), SAVE     :: file_fe_lines
    REAL(KIND=DP)                :: dum

    SELECT CASE (state)

    !----------------------------------------------------------------------------
    ! if its the first time, then open the different files and write the headers
    !----------------------------------------------------------------------------
    CASE ( "first" )
       !--- MHD variables + grains ---
       !------------------------------
       file_phys = GET_FILE_NUMBER()
       name_file_phys = TRIM(out_dir) // TRIM(name_file_phys)
       OPEN(file_phys,file = name_file_phys,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       ! mhd variables
       WRITE(file_phys,format_header,ADVANCE='NO')'distance'
       WRITE(file_phys,format_header,ADVANCE='NO')'timeN'
       WRITE(file_phys,format_header,ADVANCE='NO')'timeI'
       WRITE(file_phys,format_header,ADVANCE='NO')'nH'
       WRITE(file_phys,format_header,ADVANCE='NO')'NH'
       WRITE(file_phys,format_header,ADVANCE='NO')'NH2'
       WRITE(file_phys,format_header,ADVANCE='NO')'NCO'
       WRITE(file_phys,format_header,ADVANCE='NO')'Av'
       WRITE(file_phys,format_header,ADVANCE='NO')'op_H2'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vn'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vi'
       WRITE(file_phys,format_header,ADVANCE='NO')'VCH'     !CH!
       WRITE(file_phys,format_header,ADVANCE='NO')'VS '     !S !
       WRITE(file_phys,format_header,ADVANCE='NO')'VSH'     !SH!
       WRITE(file_phys,format_header,ADVANCE='NO')'Tn'
       WRITE(file_phys,format_header,ADVANCE='NO')'Ti'
       WRITE(file_phys,format_header,ADVANCE='NO')'Te'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vsound'
       WRITE(file_phys,format_header,ADVANCE='NO')'Vmagnet'
       WRITE(file_phys,format_header,ADVANCE='NO')'VAlfven'
       WRITE(file_phys,format_header,ADVANCE='NO')'rhoN'
       WRITE(file_phys,format_header,ADVANCE='NO')'rhoI'
       WRITE(file_phys,format_header,ADVANCE='NO')'rhoNEG'
       WRITE(file_phys,format_header,ADVANCE='NO')'DensityN'
       WRITE(file_phys,format_header,ADVANCE='NO')'DensityI'
       WRITE(file_phys,format_header,ADVANCE='NO')'DensityNEG'
       WRITE(file_phys,format_header,ADVANCE='NO')'muN'
       WRITE(file_phys,format_header,ADVANCE='NO')'muI'
       WRITE(file_phys,format_header,ADVANCE='NO')'muA'
       WRITE(file_phys,format_header,ADVANCE='NO')'muNEG'
       WRITE(file_phys,format_header,ADVANCE='NO')'-dVn/dz'
       WRITE(file_phys,format_header,ADVANCE='NO')'-dVi/dz'
       ! grains
       WRITE(file_phys,format_header,ADVANCE='NO')'Tgrain'
       WRITE(file_phys,format_header,ADVANCE='NO')'Teff_gr'
       WRITE(file_phys,format_header,ADVANCE='NO')'Mgrc'
       WRITE(file_phys,format_header,ADVANCE='NO')'Mgrm'
       WRITE(file_phys,format_header,ADVANCE='NO')'d-ero'
       WRITE(file_phys,format_header,ADVANCE='NO')'d-ads'
       WRITE(file_phys,format_header,ADVANCE='NO')'rsq-core'
       WRITE(file_phys,format_header,ADVANCE='NO')'rsq-mant'
       WRITE(file_phys,format_header,ADVANCE='NO')'ab-cor'
       WRITE(file_phys,format_header,ADVANCE='NO')'ab-ads'
       WRITE(file_phys,format_header,ADVANCE='NO')'Nlayers'
       WRITE(file_phys,format_header,ADVANCE='NO')'densgrc'
       WRITE(file_phys,format_header,ADVANCE='NO')'compr_n'
       WRITE(file_phys,format_header,ADVANCE='NO')'compr_i'
       ! Others
       WRITE(file_phys,format_header,ADVANCE='NO')'phH2'
       WRITE(file_phys,format_header,ADVANCE='NO')'phCO'
       WRITE(file_phys,format_header,ADVANCE='NO')'-DissH2_t'
       WRITE(file_phys,format_header,ADVANCE='NO')'F_gr_H2'
       WRITE(file_phys,format_header,ADVANCE='NO')'E_gr_H2_K'
       WRITE(file_phys,format_header,ADVANCE='NO')'x(H2)'
       WRITE(file_phys,format_header,ADVANCE='NO')'x(H)'
       WRITE(file_phys,format_header,ADVANCE='NO')'x(H+)'
       WRITE(file_phys,format_header,ADVANCE='NO')'sum(neut)'
       WRITE(file_phys,format_header,ADVANCE='NO')'sum(ions)'
       WRITE(file_phys,format_header,ADVANCE='NO')'sum(nega)'
       WRITE(file_phys,format_header,ADVANCE='NO')'grad_V'               !  Used only for J shocks
       WRITE(file_phys,format_header,ADVANCE='NO')'fluph'
       WRITE(file_phys,*)

       !--- chemical species ---
       !-------------------------------------------------
       file_speci = GET_FILE_NUMBER()
       name_file_speci = TRIM(out_dir) // TRIM(name_file_speci)
       OPEN(file_speci,file=name_file_speci,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       ! write column number
       WRITE(file_speci,'(A1)',ADVANCE='NO') '#'
       !WRITE(file_speci,'(2X,I16)',ADVANCE='NO') 1
       WRITE(file_speci,'(2X,I23)',ADVANCE='NO') 1
       DO i=1,Nspec+12+Nelements
          WRITE(file_speci,'(3X,I16)',ADVANCE='NO') i+1
       ENDDO
       WRITE(file_speci,*)
       
       ! densities (cm-3) or  fractional abundances
       WRITE(file_speci,'(A1)',ADVANCE='NO') ' '
       WRITE(file_speci,format_header,ADVANCE='NO')'distance'
       WRITE(file_speci,format_header,ADVANCE='NO')'timeN'
       WRITE(file_speci,format_header,ADVANCE='NO')'Tn'
       WRITE(file_speci,format_header,ADVANCE='NO')'nH'
       WRITE(file_speci,format_header,ADVANCE='NO')'Av'
       WRITE(file_speci,format_header,ADVANCE='NO')'NH'
       WRITE(file_speci,format_header,ADVANCE='NO')'NH2'
       WRITE(file_speci,format_header,ADVANCE='NO')'invAv'
       ! Ions
       SELECT CASE(TRIM(speci_out))
       CASE('AD','CD')
          str = 'n(I)'
       CASE('FD')
          str = 'x(I)'
       CASE DEFAULT
          STOP "*** Error in output specification ***"
       END SELECT
       WRITE(file_speci,format_header,ADVANCE='NO')TRIM(str)
       ! Other species
       DO i=1,Nspec
          SELECT CASE(TRIM(speci_out))
          CASE('AD','CD')
             str = 'n(' // TRIM(speci(i)%name) // ')'
          CASE('FD')
             str = 'x(' // TRIM(speci(i)%name) // ')'
          CASE DEFAULT
             STOP "*** Error in output specification ***"
          END SELECT
          WRITE(file_speci,format_header,ADVANCE='NO')TRIM(str)
       END DO
       ! Conservation laws
       WRITE(file_speci,format_header,ADVANCE='NO')'compr_n'
       WRITE(file_speci,format_header,ADVANCE='NO')'compr_i'
       DO i=1,Nelements
          WRITE(file_speci,format_header,ADVANCE='NO')"elab_"//TRIM(elements(i)%name)
       ENDDO
       WRITE(file_speci,format_header,ADVANCE='NO')"ab_"//"PAH"
       WRITE(file_speci,format_header,ADVANCE='NO')'sum(i)-sum(neg)'
       WRITE(file_speci,*)

       !--- column-densities ---
       !-------------------------------------------------
       file_coldens = GET_FILE_NUMBER()
       name_file_coldens = TRIM(out_dir) // TRIM(name_file_coldens)
       OPEN(file_coldens,file=name_file_coldens,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       !  column densities (cm-2)
       WRITE(file_coldens,format_header,ADVANCE='NO')'distance'
       WRITE(file_coldens,format_header,ADVANCE='NO')'timeN'
       WRITE(file_coldens,format_header,ADVANCE='NO')'Tn'
       WRITE(file_coldens,format_header,ADVANCE='NO')'nH'
       WRITE(file_coldens,format_header,ADVANCE='NO')'NH'
       WRITE(file_coldens,format_header,ADVANCE='NO')'NH2'
       WRITE(file_coldens,format_header,ADVANCE='NO')'NCO'
       WRITE(file_coldens,format_header,ADVANCE='NO')'Av'
       DO i=1,Nspec
          str = 'N(' // TRIM(speci(i)%name) // ')'
          WRITE(file_coldens,format_header,ADVANCE='NO')TRIM(str)
       END DO
       WRITE(file_coldens,*)


       !--- H2 levels and H2 lines ---
       !------------------------------
       ! populations
       file_H2_lev = GET_FILE_NUMBER()
       name_file_H2_lev = TRIM(out_dir) // TRIM(name_file_H2_lev)
       OPEN(file_H2_lev,file=name_file_H2_lev,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_H2_lev,format_header,ADVANCE='NO')'distance'
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'timeN'
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'Av'
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'Tn'
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'H2tot'
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'H2elec'
       ! densities (cm-3)
       DO i=1,NH2_lev
          WRITE(str,'(I2)') H2_lev(i)%V
          WRITE(str2,'(I2)') H2_lev(i)%J
          WRITE(file_H2_lev,format_header,ADVANCE='NO') '(' // TRIM(ADJUSTL(str)) // &
               ',' // TRIM(ADJUSTL(str2)) // ')'
       END DO
       WRITE(file_H2_lev,format_header,ADVANCE='NO')'sum_lev'
       WRITE(file_H2_lev,*)

       ! intensities (erg/s/cm2/sr)
       file_H2_line = GET_FILE_NUMBER()
       name_file_H2_line = TRIM(out_dir) // TRIM(name_file_H2_line)
       OPEN(file_H2_line,file=name_file_H2_line,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_H2_line,format_header,ADVANCE='NO')'distance'
       WRITE(file_H2_line,format_header,ADVANCE='NO')'timeN'
       WRITE(file_H2_line,format_header,ADVANCE='NO')'Tn'
       DO i=1,NH2_lines_out
          WRITE(file_H2_line,format_header,ADVANCE='NO')TRIM(H2_lines(i)%name)
       END DO
       WRITE(file_H2_line,format_header,ADVANCE='NO')'total_H2_em'
       ! WRITE(file_H2_line,format_header,ADVANCE='NO')'h2excit_n_coh'
       ! WRITE(file_H2_line,format_header,ADVANCE='NO')'h2excit_n_add'
       ! WRITE(file_H2_line,format_header,ADVANCE='NO')'h2excit_i_add'
       ! WRITE(file_H2_line,format_header,ADVANCE='NO')'h2excit_e_add'
       ! WRITE(file_H2_line,format_header,ADVANCE='NO')'h2excit_e_LW'
       WRITE(file_H2_line,*)

       !---------BT/BG 05-19: consistent outputs--------
       !------------- Themal balance -------------
       !--- heat/cooling rates (erg/s/cm3) ---
       !--------------------------------------
       file_therm_bal = GET_FILE_NUMBER()
       name_file_therm_bal = TRIM(out_dir) // TRIM(name_file_therm_bal)
       OPEN(file_therm_bal,file=name_file_therm_bal,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'distance'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'timeN'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'Tn'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'Ti'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'Te'

       !line radiative processes
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_tot'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_molec'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_H2'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_13CO'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_OH'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_NH3'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_rCO'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_vCO'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_CO'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_roH2O'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_rpH2O'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_vH2O'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_H2O'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_molec'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_H2'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_molec'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_H2'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_CO'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_atoms'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Hat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Oat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Op'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Cat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Cp'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Nat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Np'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Sat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Sp'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Siat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_n_Sip'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_atoms'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Hat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Oat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Op'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Cat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Cp'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Nat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Np'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Sat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Sp'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Siat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_i_Sip'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_atoms'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Hat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Oat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Op'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Cat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Cp'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Nat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Np'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Sat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Sp'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Siat'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'line_rad_e_Sip'

       !energy defect of any reaction (including photoreactions)
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_tot'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_n'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_i'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_e'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_n_diss_H2'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_i_diss_H2'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_e_diss_H2'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_n_phgas'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_i_phgas'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_e_phgas'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_n_phgrn'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_i_phgrn'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_e_phgrn'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_n_cosmic'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_i_cosmic'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_e_cosmic'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_n_other'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_i_other'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'chem_DE_e_other'

       !elastic scattering btw fluids due to drift or thermalization
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'elast_scat_tot'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'elast_scat_n'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'elast_scat_i'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'elast_scat_e'

       !exchange of internal energy between fluids
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'exch_eint_tot'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'exch_eint_n'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'exch_eint_i'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'exch_eint_e'

       !thermalization of the gas with grains
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'therm_grain_tot'

       !mechanical processes
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_tot'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_n'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_i'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_e'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_n_chem'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_i_chem'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_e_chem'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_n_compr'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_i_compr'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_e_compr'

       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_n_visc'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_i_visc'
       WRITE(file_therm_bal,format_header,ADVANCE='NO')'mech_trsf_e_visc'
       WRITE(file_therm_bal,*)


       !--- energetics ---
       !------------------
       file_energetics = GET_FILE_NUMBER()
       name_file_energetics = TRIM(out_dir) // TRIM(name_file_energetics)
       OPEN(file_energetics,file=name_file_energetics,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE(file_energetics,format_header,ADVANCE='NO')'distance'
       WRITE(file_energetics,format_header,ADVANCE='NO')'timeN'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Tn'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Mass'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Mass_cons'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_kin'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_the'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_mag'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment_vis'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Moment'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Mom_cons'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_kin'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_the'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_int'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_mag'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy_vis'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energy'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energ_gain'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energ_cons'
       WRITE(file_energetics,format_header,ADVANCE='NO')'Energ_pho'
       WRITE(file_energetics,format_header,ADVANCE='NO')'magn_corr'
       WRITE(file_energetics,*)

       !--- intensities---
       !------------------
       file_intensity = GET_FILE_NUMBER()
       name_file_intensity = TRIM(out_dir) // TRIM(name_file_intensity)
       OPEN(file_intensity,file=name_file_intensity,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')
 
       WRITE(file_intensity,format_header,ADVANCE='NO')'distance'
       WRITE(file_intensity,format_header,ADVANCE='NO')'timeN'
       WRITE(file_intensity,format_header,ADVANCE='NO')'Tn'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(158m)'                   ! 'C+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2324.7A)'                ! 'C+(3-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2323.5A)'                ! 'C+(4-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2328.1A)'                ! 'C+(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2326.9A)'                ! 'C+(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C+(2325.4A)'                ! 'C+(5-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'Si+(34.8m)'                 ! 'Si+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'H(1215.7A)'                 ! 'H(2-1)+(4-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'H(1025.7A)'                 ! 'H(5-1)+(6-1)+(7-1)+(8-1)+(9-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(609.8m)'                  ! 'C(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(370.4m)'                  ! 'C(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(9850A)'                   ! 'C(4-3)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'C(9824A)'                   ! 'C(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'Si(129.7m)'                 ! 'Si(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'Si(68.5m)'                  ! 'Si(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(63.2m)'                   ! 'O(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(145.3m)'                  ! 'O(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(6300A)'                   ! 'O(4-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'O(6363A)'                   ! 'O(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'S+(6731A)'                  ! 'S+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'S+(6716A)'                  ! 'S+(3-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(205.3m)'                 ! 'N+(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(121.8m)'                 ! 'N+(3-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(6527A)'                  ! 'N+(4-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(6548A)'                  ! 'N+(4-2)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N+(6583A)'                  ! 'N+(4-3)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N(5200A)'                   ! 'N(2-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'N(5197A)'                   ! 'N(3-1)'
       WRITE(file_intensity,format_header,ADVANCE='NO') 'S(25.2m)'                   ! 'S(2-1)'       
       WRITE(file_intensity,*)

       !--- populations---
       !------------------
       file_populations = GET_FILE_NUMBER()
       name_file_populations = TRIM(out_dir) // TRIM(name_file_populations)
       OPEN(file_populations,file=name_file_populations,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')
 
       WRITE(file_populations,format_header,ADVANCE='NO')'distance'
       WRITE(file_populations,format_header,ADVANCE='NO')'timeN'
       WRITE(file_populations,format_header,ADVANCE='NO')'Tn'
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(3P-J=0)'             ! C, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(3P-J=1)'             ! C, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(3P-J=2)'             ! C, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(1D-J=2)'             ! C, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'C(1S-J=0)'             ! C, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(4S-J=3/2)'           ! N, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2D-J=5/2)'           ! N, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2D-J=3/2)'           ! N, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2P-J=1/2)'           ! N, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'N(2P-J=3/2)'           ! N, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(3P-J=2)'             ! O, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(3P-J=1)'             ! O, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(3P-J=0)'             ! O, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(1D-J=2)'             ! O, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'O(1S-J=0)'             ! O, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(3P-J=2)'             ! S, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(3P-J=1)'             ! S, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(3P-J=0)'             ! S, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(1D-J=2)'             ! S, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'S(1S-J=0)'             ! S, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(3P-J=0)'            ! Si, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(3P-J=1)'            ! Si, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(3P-J=2)'            ! Si, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(1D-J=2)'            ! Si, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si(1S-J=0)'            ! Si, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(2P-J=1/2)'          ! C+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(2P-J=3/2)'          ! C+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(4P-J=1/2)'          ! C+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(4P-J=3/2)'          ! C+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'C+(4P-J=5/2)'          ! C+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(3P-J=0)'            ! N+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(3P-J=1)'            ! N+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(3P-J=2)'            ! N+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(1D-J=2)'            ! N+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'N+(1S-J=0)'            ! N+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(4S-J=3/2)'          ! O+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2D-J=5/2)'          ! O+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2D-J=3/2)'          ! O+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2P-J=3/2)'          ! O+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'O+(2P-J=1/2)'          ! O+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(4S-J=3/2)'          ! S+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2D-J=3/2)'          ! S+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2D-J=5/2)'          ! S+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2P-J=1/2)'          ! S+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'S+(2P-J=3/2)'          ! S+, lev = 5
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(2P-J=1/2)'         ! Si+, lev = 1
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(2P-J=3/2)'         ! Si+, lev = 2
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(4P-J=1/2)'         ! Si+, lev = 3
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(4P-J=3/2)'         ! Si+, lev = 4
       WRITE(file_populations,format_header,ADVANCE='NO') 'Si+(4P-J=5/2)'         ! Si+, lev = 5
       WRITE(file_populations,*)

       !--- Fe+ level populations ---
       !-----------------------------
       ! Relative to ground state
       file_fe_pops = GET_FILE_NUMBER()
       name_file_fe_pops = TRIM(out_dir) // TRIM(name_file_fe_pops)
       OPEN(file_fe_pops,file=name_file_fe_pops,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

       WRITE (file_fe_pops, format_header, ADVANCE='NO')'Distance'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'timeN'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'Tn'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a6D-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4F-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4D-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4P-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4P-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4P-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4P-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4P-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4P-J1'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J13'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J11'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4H-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4F-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J11'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J9'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'a4G-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J7'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J5'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J3'
       WRITE (file_fe_pops, format_header, ADVANCE='NO')'b4D-J1'
       WRITE (file_fe_pops,*)

       !--- [FeII] lines ---
       !--------------------
       ! Intensities (erg//s/cm2/sr)
       file_fe_lines = GET_FILE_NUMBER()
       name_file_fe_lines = TRIM(out_dir) // TRIM(name_file_fe_lines)
       OPEN(file_fe_lines,file=name_file_fe_lines,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'Distance'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'timeN'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'Tn'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.248'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.275'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.271'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.279'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.295'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.298'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.321'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.328'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.534'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.600'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.644'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.664'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.677'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.711'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.745'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.798'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.800'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'1.810'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'17.936'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'25.988'
       WRITE (file_fe_lines, format_header, ADVANCE='NO')'35.777'
       WRITE (file_fe_lines,*)

    !---------------------------------------------------------
    ! write blanks in output files (useful for CJ type)
    !---------------------------------------------------------
    CASE ( "blank" )

       WRITE(file_phys,*)
       WRITE(file_phys,*)
       WRITE(file_speci,*)
       WRITE(file_speci,*)
       WRITE(file_coldens,*)
       WRITE(file_coldens,*)
       WRITE(file_H2_line,*)
       WRITE(file_H2_line,*)
       WRITE(file_therm_bal,*)
       WRITE(file_therm_bal,*)
       WRITE(file_energetics,*)
       WRITE(file_energetics,*)
       WRITE(file_intensity,*)
       WRITE(file_intensity,*)
       WRITE(file_populations,*)
       WRITE(file_populations,*)
       WRITE(file_fe_pops,*)
       WRITE(file_fe_pops,*)
       WRITE(file_fe_lines,*)
       WRITE(file_fe_lines,*)

    !---------------------------------------------------------
    ! if it's the last calculation step, then close the files
    !---------------------------------------------------------
    CASE ( "final" )

       CLOSE(file_phys)
       CLOSE(file_speci)
       CLOSE(file_coldens)
       CLOSE(file_H2_lev)
       CLOSE(file_H2_line)
       CLOSE(file_therm_bal)
       CLOSE(file_energetics)
       CLOSE(file_intensity)
       CLOSE(file_populations)
       CLOSE(file_fe_pops)
       CLOSE(file_fe_lines)

    CASE DEFAULT

       !-------------------------------------------
       ! MHD variables + grains
       !-------------------------------------------
       ! mhd variables
       WRITE(file_phys,format_lng,ADVANCE='NO')distance
       WRITE(file_phys,format_lng,ADVANCE='NO')timeN
       WRITE(file_phys,format_lng,ADVANCE='NO')timeI
       WRITE(file_phys,format_lng,ADVANCE='NO')nH
       WRITE(file_phys,format_lng,ADVANCE='NO')coldens_h
       WRITE(file_phys,format_lng,ADVANCE='NO')coldens_h2
       WRITE(file_phys,format_lng,ADVANCE='NO')coldens_co
       WRITE(file_phys,format_lng,ADVANCE='NO')Av
       WRITE(file_phys,format_lng,ADVANCE='NO')op_H2
       WRITE(file_phys,format_lng,ADVANCE='NO')Vn
       WRITE(file_phys,format_lng,ADVANCE='NO')Vi
       WRITE(file_phys,format_lng,ADVANCE='NO')V_CH      !CH!
       WRITE(file_phys,format_lng,ADVANCE='NO')V_S       !S !
       WRITE(file_phys,format_lng,ADVANCE='NO')V_SH      !SH!
       WRITE(file_phys,format_lng,ADVANCE='NO')Tn
       WRITE(file_phys,format_lng,ADVANCE='NO')Ti
       WRITE(file_phys,format_lng,ADVANCE='NO')Te
       WRITE(file_phys,format_lng,ADVANCE='NO')Vsound
       WRITE(file_phys,format_lng,ADVANCE='NO')Vmagnet
       WRITE(file_phys,format_lng,ADVANCE='NO')Valfven
       WRITE(file_phys,format_lng,ADVANCE='NO')rhoN
       WRITE(file_phys,format_lng,ADVANCE='NO')rhoI
       WRITE(file_phys,format_lng,ADVANCE='NO')rhoNEG
       WRITE(file_phys,format_lng,ADVANCE='NO')DensityN
       WRITE(file_phys,format_lng,ADVANCE='NO')DensityI
       WRITE(file_phys,format_lng,ADVANCE='NO')DensityNEG
       WRITE(file_phys,format_lng,ADVANCE='NO')muN
       WRITE(file_phys,format_lng,ADVANCE='NO')muI
       WRITE(file_phys,format_lng,ADVANCE='NO')muA
       WRITE(file_phys,format_lng,ADVANCE='NO')muNEG
       WRITE(file_phys,format_lng,ADVANCE='NO')-dVn
       WRITE(file_phys,format_lng,ADVANCE='NO')-dVi
       ! grains
       WRITE(file_phys,format_lng,ADVANCE='NO')Tgrain
       WRITE(file_phys,format_lng,ADVANCE='NO')Teff_grain
       WRITE(file_phys,format_lng,ADVANCE='NO')Mgrc
       WRITE(file_phys,format_lng,ADVANCE='NO')Mgrm
       WRITE(file_phys,format_lng,ADVANCE='NO')d_ero
       WRITE(file_phys,format_lng,ADVANCE='NO')d_ads
       WRITE(file_phys,format_lng,ADVANCE='NO')rsq_grc
       WRITE(file_phys,format_lng,ADVANCE='NO')rsq_grm
       WRITE(file_phys,format_lng,ADVANCE='NO')ab_cor
       WRITE(file_phys,format_lng,ADVANCE='NO')ab_ads
       WRITE(file_phys,format_lng,ADVANCE='NO')Nlayers
       WRITE(file_phys,format_lng,ADVANCE='NO')dens_grc
       WRITE(file_phys,format_lng,ADVANCE='NO')compr_n
       WRITE(file_phys,format_lng,ADVANCE='NO')compr_i

       ! Others
       WRITE(file_phys,format_lng,ADVANCE='NO')photodiss_H2
       WRITE(file_phys,format_lng,ADVANCE='NO')photodiss_CO
       WRITE(file_phys,format_lng,ADVANCE='NO')-Sel_tot_H2
       WRITE(file_phys,format_lng,ADVANCE='NO')For_gr_H2
       WRITE(file_phys,format_lng,ADVANCE='NO')H2_int_E
       WRITE(file_phys,format_lng,ADVANCE='NO')speci(ind_H2)%Density/nH
       WRITE(file_phys,format_lng,ADVANCE='NO')speci(ind_H)%Density/nH
       WRITE(file_phys,format_lng,ADVANCE='NO')speci(ind_Hplus)%Density/nH
       WRITE(file_phys,format_lng,ADVANCE='NO')DensityN
       WRITE(file_phys,format_lng,ADVANCE='NO')DensityI
       WRITE(file_phys,format_lng,ADVANCE='NO')DensityNeg
       IF (shock_type == "J") THEN
         WRITE(file_phys,format_lng,ADVANCE='NO')grad_V
       ELSE
         WRITE(file_phys,format_lng,ADVANCE='NO')-dVn
       ENDIF
       WRITE(file_phys,format_lng,ADVANCE='NO')fluph

       ! Tests
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux1_dvn_save_tot
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux1_dvn_save_Sn
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux1_dvn_save_Vn
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux1_dvn_save_Bn
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux1_dvn_save_mol
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux2_dvn_save_tot
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux2_dvn_save_Cs
       ! WRITE(file_phys,format_lng,ADVANCE='NO') aux2_dvn_save_Vn
       WRITE(file_phys,*)

       !-----------------
       ! chemical species
       !-----------------
       WRITE(file_speci,format_lng,ADVANCE='NO')distance
       WRITE(file_speci,format_lng,ADVANCE='NO')timeN
       WRITE(file_speci,format_lng,ADVANCE='NO')Tn
       WRITE(file_speci,format_lng,ADVANCE='NO')nH
       WRITE(file_speci,format_lng,ADVANCE='NO')Av
       WRITE(file_speci,format_lng,ADVANCE='NO')coldens_h
       WRITE(file_speci,format_lng,ADVANCE='NO')coldens_h2
       WRITE(file_speci,format_lng,ADVANCE='NO')inv_Av_fac
       !Ions
       SELECT CASE(TRIM(speci_out))
       CASE('AD','CD')
          ! densities (cm-3)
          WRITE(file_speci,format_lng,ADVANCE='NO')DensityI
       CASE('FD')
          ! fractional abundances
          WRITE(file_speci,format_lng,ADVANCE='NO')DensityI/nH
       CASE DEFAULT
          STOP "*** Error in output specification ***"
       END SELECT
       ! Other species
       SELECT CASE(TRIM(speci_out))
       CASE('AD','CD')
          ! densities (cm-3)
          DO i=1,Nspec
             WRITE(file_speci,format_lng,ADVANCE='NO')speci(i)%Density
          END DO
       CASE('FD')
          ! fractional abundances
          DO i=1,Nspec
             WRITE(file_speci,format_lng,ADVANCE='NO')speci(i)%Density/nH
          END DO
       CASE DEFAULT
          STOP "*** Error in output specification ***"
       END SELECT

       WRITE(file_speci,format_lng,ADVANCE='NO')compr_n
       WRITE(file_speci,format_lng,ADVANCE='NO')compr_i
       DO i=1,Nelements
           dum = 0.0_dp
           DO ii = 1, Nspec
              IF ( speci(ii)%fluid == 1 ) THEN
                 dum = dum + DBLE(speci(ii)%formula(i)) * speci(ii)%density / compr_n
              ELSE
                 dum = dum + DBLE(speci(ii)%formula(i)) * speci(ii)%density / compr_i
              ENDIF
           ENDDO
           dum = dum / elements(ind_elem_H)%Dens_init
           WRITE(file_speci,format_lng,ADVANCE='NO')dum
       ENDDO
       dum = ( speci(ind_PAH0)%density / compr_n &
             + speci(ind_PAHp)%density / compr_i &
             + speci(ind_PAHm)%density / compr_i ) &
             & / elements(ind_elem_H)%Dens_init
       WRITE(file_speci,format_lng,ADVANCE='NO')dum
       WRITE(file_speci,format_lng,ADVANCE='NO')DensityI/nH-DensityNeg/nH
       WRITE(file_speci,*)

       !-----------------
       ! column-densities
       !-----------------
       WRITE(file_coldens,format_lng,ADVANCE='NO')distance
       WRITE(file_coldens,format_lng,ADVANCE='NO')timeN
       WRITE(file_coldens,format_lng,ADVANCE='NO')Tn
       WRITE(file_coldens,format_lng,ADVANCE='NO')nH
       WRITE(file_coldens,format_lng,ADVANCE='NO')coldens_h
       WRITE(file_coldens,format_lng,ADVANCE='NO')coldens_h2
       WRITE(file_coldens,format_lng,ADVANCE='NO')coldens_co
       WRITE(file_coldens,format_lng,ADVANCE='NO')Av
       ! column densities (cm-2)
       DO i=1,Nspec
          WRITE(file_coldens,format_lng,ADVANCE='NO')speci(i)%Col_dens
       END DO
       WRITE(file_coldens,*)

       !----------------------
       ! H2 levels + H2 lines
       !----------------------
       WRITE(file_H2_lev,format_lng,ADVANCE='NO') distance
       WRITE(file_H2_lev,format_lng,ADVANCE='NO') timeN
       WRITE(file_H2_lev,format_lng,ADVANCE='NO') Av
       WRITE(file_H2_lev,format_lng,ADVANCE='NO') Tn

       SELECT CASE(TRIM(H2_out))
       CASE('AD')
          ! densities (cm-3)
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')speci(ind_H2)%Density
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')h2elec
          DO i=1,NH2_lev
             WRITE(file_H2_lev,format_lng,ADVANCE='NO')H2_lev(i)%Density
          END DO
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')SUM(H2_lev(1:NH2_lev)%Density)
       CASE('CD')
          ! column densities (cm-2)
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')speci(ind_H2)%Col_dens
          DO i=1,NH2_lev
             WRITE(file_H2_lev,format_lng,ADVANCE='NO')H2_lev(i)%Col_dens
          END DO
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')SUM(H2_lev(1:NH2_lev)%Col_dens)
       CASE('ln(N/g)')
          ! excitation diagram
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')speci(ind_H2)%Col_dens
          DO i=1,NH2_lev
             WRITE(file_H2_lev,format_lng,ADVANCE='NO')log(H2_lev(i)%Col_dens/H2_lev(i)%weight)
          END DO
          WRITE(file_H2_lev,format_lng,ADVANCE='NO')log(SUM(H2_lev(1:NH2_lev)%Col_dens/H2_lev(1:NH2_lev)%weight))
       CASE DEFAULT
          STOP "*** Error in output specification ***"
       END SELECT
       WRITE(file_H2_lev,*)

       WRITE(file_H2_line,format_lng,ADVANCE='NO') distance
       WRITE(file_H2_line,format_lng,ADVANCE='NO') timeN
       WRITE(file_H2_line,format_lng,ADVANCE='NO') Tn
       SELECT CASE(TRIM(line_out))
       CASE('local')
          ! emissivities (erg/s/cm-3)
          dum = 0.0_dp
          DO i=1,NH2_lines_out
             WRITE(file_H2_line,format_lng,ADVANCE='NO')H2_lines(i)%emiss
             dum = dum + H2_lines(i)%emiss
          END DO
          WRITE(file_H2_line,format_lng,ADVANCE='NO')dum
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_n_coh
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_n_add
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_i_add
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_e_add
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_e_LW
       CASE('integrated')
          ! intensities (erg/s/cm2/sr)
          dum = 0.0_dp
          DO i=1,NH2_lines_out
             WRITE(file_H2_line,format_lng,ADVANCE='NO')H2_lines(i)%intensity
             dum = dum + H2_lines(i)%intensity
          END DO
          WRITE(file_H2_line,format_lng,ADVANCE='NO')dum
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_n_coh
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_n_add
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_i_add
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_e_add
          ! WRITE(file_H2_line,format_lng,ADVANCE='NO')h2excit_e_LW
       CASE DEFAULT
          STOP "*** Error in output specification ***"
       END SELECT
       WRITE(file_H2_line,*)

       !---------------------------------------------------
       ! heating and cooling rates (erg/s/cm3)
       !---------------------------------------------------
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')distance
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')timeN
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')Tn
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')Ti
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')Te

       !line radiative processes
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_tot
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e

       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_molec
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_H2
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_13CO
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_OH
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_NH3
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_rCO
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_vCO
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_CO
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_roH2O
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_rpH2O
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_vH2O
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_H2O

       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_molec
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_H2 

       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_molec
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_H2 
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_CO 

       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_atoms
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Hat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Oat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Op
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Cat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Cp
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Nat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Np
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Sat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Sp
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Siat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_n_Sip

       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_atoms
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Hat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Oat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Op
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Cat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Cp
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Nat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Np
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Sat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Sp
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Siat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_i_Sip

       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_atoms
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Hat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Oat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Op
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Cat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Cp
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Nat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Np
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Sat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Sp
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Siat
       WRITE(file_therm_bal, format_lng,ADVANCE='NO')line_rad_e_Sip

       !energy defect of any reaction (including photoreactions)
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_tot
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_n
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_i
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_e

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_n_diss_H2
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_i_diss_H2
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_e_diss_H2

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_n_phgas
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_i_phgas
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_e_phgas

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_n_phgrn
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_i_phgrn
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_e_phgrn

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_n_cosmic
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_i_cosmic
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_e_cosmic

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_n_other
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_i_other
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') chem_DE_e_other

       !elastic scattering btw fluids due to drift or thermalization
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') elast_scat_tot
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') elast_scat_n
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') elast_scat_i
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') elast_scat_e

       !exchange of internal energy between fluids
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') exch_eint_tot
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') exch_eint_n
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') exch_eint_i
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') exch_eint_e

       !thermalization of the gas with grains
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') therm_grain_tot

       !mechanical processes
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_tot
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_n
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_i
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_e

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_n_chem
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_i_chem
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_e_chem

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_n_compr
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_i_compr
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_e_compr

       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_n_visc
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_i_visc
       WRITE(file_therm_bal, format_lng,ADVANCE='NO') mech_trsf_e_visc
       WRITE(file_therm_bal,*)

       !------------
       ! energetics
       !------------
       WRITE(file_energetics,format_lng,ADVANCE='NO')distance
       WRITE(file_energetics,format_lng,ADVANCE='NO')timeN
       WRITE(file_energetics,format_lng,ADVANCE='NO')Tn
       WRITE(file_energetics,format_lng,ADVANCE='NO')Mass_flux
       WRITE(file_energetics,format_lng,ADVANCE='NO')Mass_cons
       WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_kin
       WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_the
       WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_mag
       WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux_vis
       WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_flux
       WRITE(file_energetics,format_lng,ADVANCE='NO')Momentum_cons
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_kin
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_the
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_int
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_mag
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_vis
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_gain
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_cons
       WRITE(file_energetics,format_lng,ADVANCE='NO')Energy_flux_pho
       WRITE(file_energetics,format_lng,ADVANCE='NO')mag_flux_corr
       WRITE(file_energetics,*)

       !------------
       ! Intensities
       !------------
       WRITE(file_intensity,format_lng,ADVANCE='NO')distance
       WRITE(file_intensity,format_lng,ADVANCE='NO')timeN
       WRITE(file_intensity,format_lng,ADVANCE='NO')Tn
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(1)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(5)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(6)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(2)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(3)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcpl(4)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intsipl(1)
       WRITE(file_intensity,format_lng,ADVANCE='NO')inthat(1)+inthat(3)
       WRITE(file_intensity,format_lng,ADVANCE='NO')inthat(4)+inthat(6)+inthat(9)+inthat(12)+inthat(14)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(1)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(2)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(4)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intcat(5)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intsiat(1)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intsiat(2)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(2)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(1)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(6)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intoat(5)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intspl(7)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intspl(8)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(1)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(2)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(6)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(5)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intnpl(4)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intnat(6)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intnat(7)
       WRITE(file_intensity,format_lng,ADVANCE='NO')intsat(2)
       WRITE(file_intensity,*)

       !------------
       ! Populations
       !------------
       WRITE(file_populations,format_lng,ADVANCE='NO')distance
       WRITE(file_populations,format_lng,ADVANCE='NO')timeN
       WRITE(file_populations,format_lng,ADVANCE='NO')Tn
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cat(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_nat(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_oat(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sat(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_siat(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_cpl(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_npl(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_opl(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_spl(5)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(1)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(2)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(3)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(4)
       WRITE(file_populations,format_lng,ADVANCE='NO')pop_sipl(5)
       WRITE(file_populations,*)

       !--- Fe+ level populations ---
       !-----------------------------
       WRITE(file_fe_pops,format_lng,ADVANCE='NO')distance
       WRITE(file_fe_pops,format_lng,ADVANCE='NO')timeN
       WRITE(file_fe_pops,format_lng,ADVANCE='NO')Tn
       DO i=1,nlvfepl
         WRITE(file_fe_pops,format_lng,ADVANCE='NO')pop_fepl(i)
       END DO
       WRITE (file_fe_pops,*)

       !--- [FeII] lines ---
       !--------------------
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')distance
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')timeN
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')Tn
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(35)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(19)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(43)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(36)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(27)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(37)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(20)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(28)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(29)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(38)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(22)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(44)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(30)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(39)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(45)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(40)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(31)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(23)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(10)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(1)
       WRITE(file_fe_lines,format_lng,ADVANCE='NO')intfepl(18)
       WRITE (file_fe_lines,*)

    END SELECT

  END SUBROUTINE WRITE_OUTPUT


  SUBROUTINE WRITE_RADIATION_FIELD(state)
     USE MODULE_TOOLS,            ONLY : GET_FILE_NUMBER
     USE MODULE_PHYS_VAR
     USE MODULE_RADIATION
     IMPLICIT NONE

     CHARACTER(len=5), INTENT(in) :: state
     INTEGER(KIND=LONG)           :: i
     CHARACTER(LEN=*), PARAMETER  :: format_header='(A18,1X)'
     CHARACTER(LEN=*), PARAMETER  :: format_lng='(ES18.9E3,1X)'
     CHARACTER(LEN=lenfilename)   :: name_file_irf = 'mhd_irf.out'
     INTEGER(KIND=LONG), SAVE     :: file_irf

     SELECT CASE (state)

     !----------------------------------------------------------------------------
     ! if its the first time, then open the different files and write the headers
     !----------------------------------------------------------------------------
     CASE ( "first" )
        file_irf = GET_FILE_NUMBER()
        name_file_irf = TRIM(out_dir) // TRIM(name_file_irf)
        OPEN(file_irf,file = name_file_irf,status='UNKNOWN',&
             access='SEQUENTIAL',form='FORMATTED',action='WRITE')
 
        ! Interstellar radiation field intensity
        WRITE(file_irf,format_header,ADVANCE='NO')'distance'
        WRITE(file_irf,format_header,ADVANCE='NO')'wavelgth'
        WRITE(file_irf,format_header,ADVANCE='NO')'I IRF'
        WRITE(file_irf,*)
 
     !---------------------------------------------------------
     ! if it's the last calculation step, then close the files
     !---------------------------------------------------------
     CASE ( "final" )
 
        CLOSE(file_irf)

     CASE DEFAULT
 
        DO i = 1, nwlg
           WRITE(file_irf,format_lng,ADVANCE='NO')distance
           WRITE(file_irf,format_lng,ADVANCE='NO')wlg(i)
           WRITE(file_irf,format_lng,ADVANCE='NO')irf(1,i)
           WRITE(file_irf,*)
        ENDDO
        WRITE(file_irf,*)

     END SELECT
      
  END SUBROUTINE WRITE_RADIATION_FIELD


  SUBROUTINE WRITE_EXCIT
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     write H2 excitation diagram
    ! subroutine/function needed :
    !     GET_FILE_NUMBER
    ! input variables :
    !     none
    !---------------------------------------------------------------------------
    USE MODULE_TECHCONFIG
    USE MODULE_TOOLS,    ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR, ONLY : NH2_lev
    USE MODULE_H2

    IMPLICIT NONE

    REAL(KIND=DP),     PARAMETER :: tiny = 1.0e-40_DP
    INTEGER(KIND=LONG)           :: i
    CHARACTER(LEN=35), PARAMETER :: format_head = '(4X,"V   J   Energy(K)   log(N/g)")'
    CHARACTER(LEN= 7), PARAMETER :: format_i    = '(i5,1X)'
    CHARACTER(LEN=13), PARAMETER :: format_r    = '(ES15.6E3,1X)'
    CHARACTER(LEN=lenfilename)   :: name_file
    INTEGER                      :: file_unit

    name_file   = 'excit.out'
    file_unit = GET_FILE_NUMBER()
    name_file = TRIM(out_dir) // TRIM(name_file)
    OPEN(file_unit,file=name_file,status='UNKNOWN',&
            access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    WRITE(file_unit,format_head)

    DO i = 1, NH2_lev
       WRITE(file_unit,format_i,ADVANCE='NO') H2_lev(i)%V
       WRITE(file_unit,format_i,ADVANCE='NO') H2_lev(i)%J
       WRITE(file_unit,format_r,ADVANCE='NO') H2_lev(i)%Energy
       WRITE(file_unit,format_r,ADVANCE='NO') log(tiny+H2_lev(i)%Col_dens/H2_lev(i)%weight)
       WRITE(file_unit,*)
    ENDDO

    CLOSE(file_unit)

  END SUBROUTINE WRITE_EXCIT


  SUBROUTINE WRITE_SPECIES
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     write final abundances (same format as input file)
    ! subroutine/function needed :
    !     GET_FILE_NUMBER
    ! input variables :
    !     none
    !---------------------------------------------------------------------------
    USE MODULE_TECHCONFIG
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    IMPLICIT NONE

    CHARACTER(len=lenfilename) :: name_file_in
    CHARACTER(len=lenfilename) :: name_file_out
    INTEGER(KIND=LONG)         :: file_in
    INTEGER(KIND=LONG)         :: file_out
    CHARACTER(len=80)          :: charact
    CHARACTER(len=34)          :: end_l
    INTEGER                    :: i
    REAL                       :: toto

    file_in = GET_FILE_NUMBER()
    name_file_in = TRIM(data_dir) // TRIM(specfile)
    OPEN(file_in,file=name_file_in,status='OLD',&
         access='SEQUENTIAL',form='FORMATTED',action='READ')

    file_out = GET_FILE_NUMBER()
    name_file_out = 'species.out'
    name_file_out = TRIM(out_dir) // TRIM(name_file_out)
    OPEN(file_out,file=name_file_out,status='REPLACE',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    ! Skip header in input file
    DO i = 1, nlines_comment-3
       READ(file_in,'(A80)') charact
    END DO

    ! Write header in output file
    WRITE(file_out,'(A55,F8.2,A18)')    "!---- list of chemical species --- Steady state at T = ",&
                                        Tn, " K ---------------"
    WRITE(file_out,'(A12,F8.2,A18)')    "!---- nH  = ", nH, " cm-3 ----------"
    WRITE(file_out,'(A12,F8.2,A18)')    "!---- Av  = ", Av, " mag -----------"
    WRITE(file_out,'(A12,F8.2,A18)')    "!---- RAD = ", RAD," ---------------"
    WRITE(file_out,'(A15,1pe10.3,A18)') "!---- N_H2_0 = ", N_H2_0," cm-2 ----"
    WRITE(file_out,'(A15,1pe10.3,A18)') "!---- N_CO_0 = ", N_CO_0," cm-2 ----"
    
    ! Copy last 3 lines of comments from input to output file.
    DO i = 1, 3
       READ(file_in,'(A80)') charact
       WRITE(file_out,'(A80)') charact
    END DO

    DO i = 1, Nspec-1
       READ(file_in,'(28X,E10.3,A34)') toto, end_l
       WRITE(file_out,'(I3,2X,A7,2X,14X,1P,D10.3,A34)') i, speci(i)%name, speci(i)%Density/nH, end_l
    END DO

    CLOSE(file_in)
    CLOSE(file_out)
  END SUBROUTINE WRITE_SPECIES


  SUBROUTINE WRITE_H2ABUND
    !---------------------------------------------------------------------------
    ! called by :
    !     MHD
    ! purpose :
    !     write final abundances of the H2 levels (same format as input file)
    ! subroutine/function needed :
    !     GET_FILE_NUMBER
    ! input variables :
    !     none
    !---------------------------------------------------------------------------
    USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
    USE MODULE_H2
    USE MODULE_PHYS_VAR
    USE MODULE_CHEMICAL_SPECIES
    IMPLICIT NONE

    CHARACTER(len=lenfilename) :: name_file_out
    INTEGER(KIND=LONG)         :: file_out
    INTEGER                    :: i
    CHARACTER (LEN=2)          :: str1, str2

    file_out = GET_FILE_NUMBER()
    name_file_out = 'h2levels.out'
    name_file_out = TRIM(out_dir) // TRIM(name_file_out)
    OPEN(file_out,file=name_file_out,status='REPLACE',&
         access='SEQUENTIAL',form='FORMATTED',action='WRITE')

    ! Write header in output file
    WRITE(file_out,'(A30)') "!============================="
    WRITE(file_out,'(A29)') "! H2* steady state abundances"
    WRITE(file_out,'(A30)') "!============================="
    WRITE(file_out,'(A15,F8.2,A2)')  "!--- T       = ", Tn,     " K"
    WRITE(file_out,'(A15,F8.2,A5)')  "!--- nH      = ", nH,     " cm-3"
    WRITE(file_out,'(A15,F8.2,A4)')  "!--- Av      = ", Av,     " mag"
    WRITE(file_out,'(A15,F8.2)')     "!--- RAD     = ", RAD
    WRITE(file_out,'(A15,ES8.2,A5)') "!--- N_H2_0  = ", N_H2_0, " cm-2"
    WRITE(file_out,'(A15,ES8.2,A5)') "!--- N_CO_0  = ", N_CO_0, " cm-2"
    WRITE(file_out,'(A15,I8)')       "!--- NH2_lev = ", NH2_lev
    WRITE(file_out,'(A30)') "!============================="
    
    ! Write the populations of H2 levels (normalized to H2)
    DO i = 1, NH2_lev
       WRITE(str1,'(I2)') H2_lev(i)%V
       WRITE(str2,'(I2)') H2_lev(i)%J
       WRITE(file_out,'(A9,3X,ES11.3E3)') 'H2(' // str1 // ',' // str2 // ')' , &
                                          H2_lev(i)%Density / speci(ind_H2)%Density
    END DO

    CLOSE(file_out)
  END SUBROUTINE WRITE_H2ABUND

END MODULE MODULE_OUTPUTS

