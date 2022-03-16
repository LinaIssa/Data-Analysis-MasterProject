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

MODULE OUTPUT_HDF5
  !*****************************************************************************
  !** The module 'OUTPUT_HDF5' contains subroutines for writting ASCII files  **
  !** which is an intermediate step to the production of hdf5 files           **
  !*****************************************************************************

  IMPLICIT NONE
  INCLUDE "precision.f90"

  PUBLIC :: WRITE_ASCII

  INTEGER, PARAMETER                            :: lencol = 1000000  ! Maximum number of characters
                                                                     ! in each line of output ASCII
  INTEGER                                       :: speidx            ! index of the species

  INTEGER                                       :: nbf               ! number of formation   reaction of a given species
  INTEGER                                       :: nbd               ! number of destruction reaction of a given species
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: form_rate         ! formation   rates of a given species at all positions
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: dest_rate         ! destruction rates of a given species at all positions
  INTEGER,  PUBLIC, ALLOCATABLE, DIMENSION(:)   :: idx_freac         ! indices of the formation   reactions
  INTEGER,  PUBLIC, ALLOCATABLE, DIMENSION(:)   :: idx_dreac         ! indices of the destruction reactions

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



  SUBROUTINE WRITE_ASCII(filedir, metatab, nhdf5, typehdf5)
    !---------------------------------------------------------------------------
    ! called by :
    !    MHD
    ! purpose :
    !    write most of the code output (+metadata, see above) in several
    !    ASCII files. Theses ASCII files will be used by shock_write_HDF5.py 
    !    to produce HDF5 outputs
    ! subroutine/function needed :
    ! input variables :
    !    filedir : name of the temporary directory containing ascii files
    !    metatab : table of metadata
    !    nhdf5 : dmiension of the metadata table
    !    typehdf5 : type of the HDF5 file ... standard ('s') or chemical ('c')
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------

    USE MODULE_TECHCONFIG
    USE MODULE_TOOLS, ONLY      : GET_FILE_NUMBER
    USE MODULE_CHEMICAL_SPECIES

    IMPLICIT none

    CHARACTER(len=lenfilename),                   INTENT(in) :: filedir   ! directory of the ASCII files
    INTEGER,                                      INTENT(in) :: nhdf5     ! dimension of the metatab table
    TYPE (METADATA),            DIMENSION(nhdf5), INTENT(in) :: metatab   ! table of metadata
    CHARACTER,                                    INTENT(in) :: typehdf5  ! 's' for standard or 'c' for chemical
    CHARACTER(len=lenfilename)                               :: filename  ! name of the ASCII file to create
    INTEGER                                                  :: ndatafile ! number of data to write per file
    INTEGER                                                  :: idx_star  ! dummy integer
    INTEGER                                                  :: i, j      ! dummy integers

    ! =====================================================================================
    ! Create a directory in ../out/modele/ to store temporary ASCII files
    ! =====================================================================================
    sh_cmd = 'sh -c "if [ ! -d '//TRIM(filedir)//" ]; then mkdir "
    sh_cmd = TRIM(sh_cmd)//" "//TRIM(filedir)//'; fi "'
    CALL SYSTEM(sh_cmd)

    ! =====================================================================================
    ! Find the number of data to write per ASCII file
    ! Find the number of ASCII files to create for each line of meta_xxxxxxxx.dat
    ! =====================================================================================
    i = 1
    DO WHILE(i <= nhdf5)
       ! ---------------------------------------------------
       ! Test if several lines in meta_standard.dat have the
       ! same output ASCII file name
       ! ---------------------------------------------------
       filename  = TRIM(metatab(i)%file)
       DO j = i, nhdf5
          IF ( TRIM(metatab(j)%file) /= TRIM(filename) ) EXIT
       ENDDO
       ndatafile = j - i
       filename  = TRIM(filedir)//TRIM(filename)
       ! ---------------------------------------------------
       ! Test if we are writing the chemical data
       ! If so, we need to open as many ASCII files as the
       ! number of species treated in the chemical network
       ! ---------------------------------------------------
       IF      ( ( INDEX(metatab(i)%IDname, 'dens_spec_* ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'frate_*     ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'drate_*     ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_frate_* ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_drate_* ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'compo_of_*  ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'mass_of_*   ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%file, '*', .false.)
          IF (   ( INDEX(metatab(i)%file, 'ChemDens*.dat ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%file, 'Compo*.dat    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%file, 'Mass*.dat     ', .false.) == 1 ) ) THEN
             DO speidx = 1, Nspec
                filename = metatab(i)%file(1:idx_star-1)//TRIM(speci(speidx)%SpecyID)//metatab(i)%file(idx_star+1:)
                filename = TRIM(filedir)//TRIM(filename)
                CALL WRITE_METADATA(filename, metatab, nhdf5, i, ndatafile, typehdf5)
             ENDDO
          ELSE
             DO speidx = 1, Nspec
                CALL FILL_RATES
                ! ------------------------------------------
                ! dont write the ascii file if the species
                ! has no formation or destruction pathways
                ! ------------------------------------------
                IF ( ( ( INDEX(metatab(i)%IDname, 'frate_*     ', .false.) == 1 ).AND.( nbf == 0 ) ).OR.&
                     ( ( INDEX(metatab(i)%IDname, 'drate_*     ', .false.) == 1 ).AND.( nbd == 0 ) ) ) THEN
                   DEALLOCATE (form_rate)
                   DEALLOCATE (dest_rate)
                   DEALLOCATE (idx_freac)
                   DEALLOCATE (idx_dreac)
                   CYCLE
                ENDIF
                filename = metatab(i)%file(1:idx_star-1)//TRIM(speci(speidx)%SpecyID)//metatab(i)%file(idx_star+1:)
                filename = TRIM(filedir)//TRIM(filename)
                CALL WRITE_METADATA(filename, metatab, nhdf5, i, ndatafile, typehdf5)
                DEALLOCATE (form_rate)
                DEALLOCATE (dest_rate)
                DEALLOCATE (idx_freac)
                DEALLOCATE (idx_dreac)
             ENDDO
          ENDIF
       ELSE
          CALL WRITE_METADATA(filename, metatab, nhdf5, i, ndatafile, typehdf5)
       ENDIF
       ! ---------------------------------------------------
       ! Skip the lines in the metadata file that we have
       ! already treated by defining ndatafile
       ! ---------------------------------------------------
       i = i + ndatafile
    ENDDO

  END SUBROUTINE WRITE_ASCII


  SUBROUTINE WRITE_METADATA(filename, metatab, nhdf5, ifil, ndatafile, typehdf5)
    !---------------------------------------------------------------------------
    ! called by :
    !    WRITE_ASCII
    ! purpose :
    !    Write a set of metadata + data in a given ascii file
    ! subroutine/function needed :
    ! input variables :
    !    filename : ascii file name
    !    metatab : table of metadata
    !    nhdf5 : dmiension of the metadata table
    !    ifil : index of the first output of the set
    !    ndatafile : number of outputs in the set to write in the filename file
    !    typehdf5 : type of the HDF5 file ... standard ('s') or chemical ('c')
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------

    USE MODULE_TECHCONFIG
    USE MODULE_CONSTANTS, ONLY : iwrtmp, amu
    USE MODULE_TOOLS, ONLY     : GET_FILE_NUMBER
    USE MODULE_CHEMICAL_SPECIES
    USE MODULE_CHEM_REACT
    USE MODULE_PHYS_VAR
    USE MODULE_PROFIL_TABLES
    USE MODULE_GRAINS
    USE MODULE_H2
    USE MODULE_LINE_EXCIT
    USE MODULE_VAR_VODE, ONLY : duration_max, Eps_V

    IMPLICIT none

    CHARACTER(len=lenfilename),         INTENT(in) :: filename        ! name of the ASCII file
    INTEGER,                            INTENT(in) :: nhdf5           ! dimension of the metatab table
    TYPE (METADATA), DIMENSION(nhdf5)              :: metatab         ! table of metadata
    INTEGER,                            INTENT(in) :: ifil            ! starting index in metatab_std
    INTEGER,                            INTENT(in) :: ndatafile       ! number of data to write in the ASCII file
    CHARACTER,                          INTENT(in) :: typehdf5        ! 's' for standard or 'c' for chemical
    INTEGER,         DIMENSION(ndatafile)          :: nsubcols        ! number of subcolumns per data
    INTEGER                                        :: nvals           ! number of values per column
    CHARACTER(len=17)                              :: strf            ! string dummy variable to write reals
    CHARACTER(len=10)                              :: stri            ! string dummy variable to write integers
    CHARACTER(len=100)                             :: strs            ! string dummy variable to strings
    CHARACTER(len=lenmetadesc)                     :: template        ! temporary string to store metadata
    CHARACTER(len=name_length)                     :: sp_lowcas       ! name of elemental species (lower case)
    CHARACTER(len=name_length)                     :: sp_uppcas       ! name of elemental species (upper case)
    CHARACTER(len=lenmetakey)                      :: key_temp        ! temporary key string
    INTEGER                                        :: idx_star        ! location of the star the quantities listed in metadata.dat
    INTEGER                                        :: i, j, k, l, j1  ! dummy indices
    INTEGER                                        :: Nstep_hdf5      ! number of steps between 2 outputs of local quantities

    fichier = TRIM(filename)
    iwrtmp = GET_FILE_NUMBER()
    OPEN(iwrtmp, FILE=fichier, status="unknown", action="write", RECL=lencol)

    !======================================================================================
    ! Find the number of steps between two outputs of local quantities
    !======================================================================================
    IF ( typehdf5 == 'c' ) THEN
       Nstep_hdf5 = Nstep_w
    ELSE
       Nstep_hdf5 = 1
    ENDIF

    !======================================================================================
    ! Test if the quantity to write in the ASCII file are local or global quantities
    !    -> deduce the number of values to write (npo + 1 or 1)
    !======================================================================================
    IF ( INDEX(metatab(ifil)%path, 'Local', .false.) /= 0 ) THEN
       IF (counter <= Npthdf5) THEN
          nvals = counter
       ELSE
          nvals = Npthdf5
       ENDIF
    ELSE
       nvals = 1
    ENDIF
    !======================================================================================
    ! Compute the number of subcolumns per data
    ! Note: if a data shouldn't be written in the output (e.g. input nH in the case of an
    !       isobaric model), then we set nsubcols = 0
    ! TODO: if F_DUST_P == 1, the grain parameter in pdr.in are not used and should not be
    !       consider as parameters in the output ascii files
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       nsubcols(i-ifil+1) = 0
       IF      ( ( INDEX(metatab(i)%IDname, 'n_*                  ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_prof_*            ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_*                 ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + Nspec
       ! ELSE IF ( ( INDEX(metatab(i)%IDname, 'coolrate_*           ', .false.) == 1 ) ) THEN
       !    nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + nspl
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'elem_*               ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + Nelements
       ! ELSE IF ( ( INDEX(metatab(i)%IDname, 'cd_lev_*             ', .false.) == 1 ).OR.&
       !           ( INDEX(metatab(i)%IDname, 'pop_*                ', .false.) == 1 ) ) THEN
       !    DO j1 = 1, nspl
       !       DO j2 = 1, spec_lin(j1)%use
       !          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + 1
       !       ENDDO
       !    ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'pop_h2_*             ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_h2_*              ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + NH2_lev
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h2_*             ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_h2_*            ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + NH2_lines_out
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h_*              ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_h_*             ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrhat
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_c_*              ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_c_*             ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrcat
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_n_*              ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_n_*             ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrnat
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_o_*              ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_o_*             ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntroat
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_s_*              ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_s_*             ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrsat
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_si_*             ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_si_*            ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrsiat
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_cp_*             ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_cp_*            ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrcpl
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_np_*             ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_np_*            ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrnpl
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_op_*             ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_op_*            ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntropl
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sp_*             ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_sp_*            ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrspl
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sip_*            ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_sip_*           ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + ntrsipl
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'frate_*              ', .false.) == 1 ) ) THEN
          DO j = 1, nbf
             nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + 1
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'drate_*              ', .false.) == 1 ) ) THEN
          DO j = 1, nbd
             nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + 1
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'compo_of_*           ', .false.) == 1 ) ) THEN
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + Nelements
       ELSE
          nsubcols(i-ifil+1) = nsubcols(i-ifil+1) + 1
       ENDIF
    ENDDO

    !======================================================================================
    ! 1 - Write the line containing all the paths in the output ASCII file
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       IF      ( ( INDEX(metatab(i)%IDname, 'dens_spec_*    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'frate_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'drate_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_frate_*    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_drate_*    ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%path, '*', .false.)
          template = metatab(i)%path(1:idx_star-1)//TRIM(speci(speidx)%Hname)//metatab(i)%path(idx_star+1:)
       ELSE
          template = metatab(i)%path
       ENDIF
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_path)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)
 
    !======================================================================================
    ! 2 - Write the line containing all the datasets in the output ASCII file
    !     Dataset names need to be build here only for lines intensities
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       IF      ( ( INDEX(metatab(i)%IDname, 'dens_spec_*    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'frate_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'drate_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_frate_*    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_drate_*    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'compo_of_*     ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'mass_of_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%datset, '*', .false.)
          template = metatab(i)%datset(1:idx_star-1)//TRIM(speci(speidx)%Hname)//metatab(i)%datset(idx_star+1:)
       ELSE
          template = metatab(i)%datset
       ENDIF
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_datset)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 3 - Write the line containing all the human readable names in the output ASCII file
    !     if a "*" is present in the pubID, generic names need to be created
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       IF      ( ( INDEX(metatab(i)%IDname, 'n_*            ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_prof_*      ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_*           ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, Nspec
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(speci(j)%Hname)//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'elem_*         ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, Nelements
             sp_uppcas = TRIM(ADJUSTL(elements(j)%name))
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(ADJUSTL(sp_uppcas))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ! ELSE IF ( ( INDEX(metatab(i)%IDname, 'cd_lev_*       ', .false.) == 1 ).OR.&
       !           ( INDEX(metatab(i)%IDname, 'pop_*          ', .false.) == 1 ) ) THEN
       !    idx_star = INDEX(metatab(i)%Hname, '*', .false.)
       !    DO j1 = 1, nspl
       !       DO j2 = 1, spec_lin(j1)%use
       !          template = metatab(i)%Hname(1:idx_star-1)//&
       !                     TRIM(spec_lin(j1)%spHname)//" "//&
       !                     TRIM(spec_lin(j1)%LeHname(j2))//&
       !                     metatab(i)%Hname(idx_star+1:)
       !          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
       !          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       !       ENDDO
       !    ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'pop_h2_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_h2_*        ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, NH2_lev
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(H2_lev(j)%Hname)//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h2_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_h2_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, NH2_lines_out
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(H2_lines(j)%Hname)//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_h_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrhat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrhat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_c_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_c_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrcat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrcat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_n_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_n_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrnat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrnat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_o_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_o_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntroat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtroat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_s_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_s_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrsat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrsat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_si_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_si_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrsiat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrsiat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_cp_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_cp_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrcpl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrcpl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_np_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_np_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrnpl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrnpl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_op_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_op_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntropl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtropl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sp_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_sp_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrspl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrspl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sip_*      ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_sip_*     ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrsipl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(namtrsipl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'dens_spec_*    ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          template = metatab(i)%Hname(1:idx_star-1)//&
                     TRIM(speci(speidx)%Hname)//&
                     metatab(i)%Hname(idx_star+1:)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'frate_*        ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, nbf
             ! WRITE(*,*)
             ! WRITE(*,*) speci(react(idx_freac(j))%R( 1))%Hname
             ! WRITE(*,*) speci(react(idx_freac(j))%R( 2))%Hname
             ! WRITE(*,*) speci(react(idx_freac(j))%P( 1))%Hname
             ! WRITE(*,*) speci(react(idx_freac(j))%P( 2))%Hname
             ! WRITE(*,*) speci(react(idx_freac(j))%P( 3))%Hname
             ! WRITE(*,*) speci(react(idx_freac(j))%P( 4))%Hname
             template = metatab(i)%Hname(1:idx_star-1)//TRIM(speci(speidx)%Hname)
             IF ( react(idx_freac(j))%R(1)  /= 0 ) template = TRIM(template)//" | "//&
                                                   TRIM(speci(react(idx_freac(j))%R( 1))%Hname)
             DO j1 = 2, 2
             IF ( react(idx_freac(j))%R(j1) /= 0 ) template = TRIM(template)//" + "//&
                                                   TRIM(speci(react(idx_freac(j))%R(j1))%Hname)
             ENDDO
             IF ( react(idx_freac(j))%P( 1) /= 0 ) template = TRIM(template)//" > "//&
                                                   TRIM(speci(react(idx_freac(j))%P( 1))%Hname)
             DO j1 = 2, 4
             IF ( react(idx_freac(j))%P(j1) /= 0 ) template = TRIM(template)//" + "//&
                                                   TRIM(speci(react(idx_freac(j))%P(j1))%Hname)
             ENDDO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'drate_*        ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, nbd
             template = metatab(i)%Hname(1:idx_star-1)//TRIM(speci(speidx)%Hname)
             IF ( react(idx_dreac(j))%R(1)  /= 0 ) template = TRIM(template)//" | "//&
                                                   TRIM(speci(react(idx_dreac(j))%R( 1))%Hname)
             DO j1 = 2, 2
             IF ( react(idx_dreac(j))%R(j1) /= 0 ) template = TRIM(template)//" + "//&
                                                   TRIM(speci(react(idx_dreac(j))%R(j1))%Hname)
             ENDDO
             IF ( react(idx_dreac(j))%P( 1) /= 0 ) template = TRIM(template)//" > "//&
                                                   TRIM(speci(react(idx_dreac(j))%P( 1))%Hname)
             DO j1 = 2, 4
             IF ( react(idx_dreac(j))%P(j1) /= 0 ) template = TRIM(template)//" + "//&
                                                   TRIM(speci(react(idx_dreac(j))%P(j1))%Hname)
             ENDDO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'tot_frate_*    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_drate_*    ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          template = metatab(i)%Hname(1:idx_star-1)//TRIM(speci(speidx)%Hname)
          DO j = 1, nsubcols(i-ifil+1)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'compo_of_*     ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, Nelements
             template = metatab(i)%Hname(1:idx_star-1)//TRIM(speci(speidx)%Hname)//' ('
             sp_uppcas = TRIM(ADJUSTL(elements(j)%name))
             template = TRIM(template)//TRIM(sp_uppcas)//')'
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'mass_of_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          template = metatab(i)%Hname(1:idx_star-1)//&
                     TRIM(speci(speidx)%Hname)//&
                     metatab(i)%Hname(idx_star+1:)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ELSE
          template = metatab(i)%Hname
          DO j = 1, nsubcols(i-ifil+1)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ENDIF
    ENDDO
    WRITE(iwrtmp,*)
 
    !======================================================================================
    ! 4 - Write the line containing all the unique ID names in the output ASCII file
    !     if a "*" is present in the pubID, generic names need to be created
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       IF      ( ( INDEX(metatab(i)%IDname, 'n_*            ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_prof_*      ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_*           ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          DO j = 1, Nspec
             template = metatab(i)%IDname(1:idx_star-1)//&
                        TRIM(speci(j)%SpecyID)//&
                        metatab(i)%IDname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'elem_*         ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          DO j = 1, Nelements
             IF ( elements(j)%name == 'H      ') sp_lowcas = 'h      '
             IF ( elements(j)%name == 'D      ') sp_lowcas = 'd      '
             IF ( elements(j)%name == 'He     ') sp_lowcas = 'he     '
             IF ( elements(j)%name == 'C      ') sp_lowcas = 'c      '
             IF ( elements(j)%name == 'N      ') sp_lowcas = 'n      '
             IF ( elements(j)%name == 'O      ') sp_lowcas = 'o      '
             IF ( elements(j)%name == 'Na     ') sp_lowcas = 'na     '
             IF ( elements(j)%name == 'Mg     ') sp_lowcas = 'mg     '
             IF ( elements(j)%name == 'Si     ') sp_lowcas = 'si     '
             IF ( elements(j)%name == 'S      ') sp_lowcas = 's      '
             IF ( elements(j)%name == 'Fe     ') sp_lowcas = 'fe     '
             IF ( elements(j)%name == 'G      ') sp_lowcas = 'g      '
             template = metatab(i)%IDname(1:idx_star-1)//&
                        TRIM(ADJUSTL(sp_lowcas))//&
                        metatab(i)%IDname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'pop_h2_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'cd_h2_*        ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          DO j = 1, NH2_lev
             template = metatab(i)%IDname(1:idx_star-1)//&
                        TRIM(H2_lev(j)%ID)//&
                        metatab(i)%IDname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h2_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_h2_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          DO j = 1, NH2_lines_out
             template = metatab(i)%IDname(1:idx_star-1)//&
                        TRIM(H2_lines(j)%ID)//&
                        metatab(i)%IDname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_h_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrhat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrhat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_c_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_c_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrcat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrcat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_n_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_n_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrnat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrnat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_o_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_o_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntroat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtroat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_s_*        ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_s_*       ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrsat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrsat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_si_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_si_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrsiat
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrsiat(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_cp_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_cp_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrcpl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrcpl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_np_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_np_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrnpl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrnpl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_op_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_op_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntropl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtropl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sp_*       ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_sp_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrspl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrspl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sip_*      ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'inta_sip_*     ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%Hname, '*', .false.)
          DO j = 1, ntrsipl
             template = metatab(i)%Hname(1:idx_star-1)//&
                        TRIM(idtrsipl(j))//&
                        metatab(i)%Hname(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_Hname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'dens_spec_*    ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          template = metatab(i)%IDname(1:idx_star-1)//&
                     TRIM(speci(speidx)%SpecyID)//&
                     metatab(i)%IDname(idx_star+1:)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'frate_*        ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          DO j = 1, nbf
             template = metatab(i)%IDname(1:idx_star-1)//TRIM(speci(speidx)%SpecyID)
             IF ( react(idx_freac(j))%R(1)  /= 0 ) template = TRIM(template)//"_"//&
                                                   TRIM(speci(react(idx_freac(j))%R( 1))%SpecyID)
             DO j1 = 2, 2
             IF ( react(idx_freac(j))%R(j1) /= 0 ) template = TRIM(template)//"_"//&
                                                   TRIM(speci(react(idx_freac(j))%R(j1))%SpecyID)
             ENDDO
             IF ( react(idx_freac(j))%P( 1) /= 0 ) template = TRIM(template)//"_gives_"//&
                                                   TRIM(speci(react(idx_freac(j))%P( 1))%SpecyID)
             DO j1 = 2, 4
             IF ( react(idx_freac(j))%P(j1) /= 0 ) template = TRIM(template)//"_"//&
                                                   TRIM(speci(react(idx_freac(j))%P(j1))%SpecyID)
             ENDDO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'drate_*        ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          DO j = 1, nbd
             template = metatab(i)%IDname(1:idx_star-1)//TRIM(speci(speidx)%SpecyID)
             IF ( react(idx_dreac(j))%R(1)  /= 0 ) template = TRIM(template)//"_"//&
                                                   TRIM(speci(react(idx_dreac(j))%R( 1))%SpecyID)
             DO j1 = 2, 2
             IF ( react(idx_dreac(j))%R(j1) /= 0 ) template = TRIM(template)//"_"//&
                                                   TRIM(speci(react(idx_dreac(j))%R(j1))%SpecyID)
             ENDDO
             IF ( react(idx_dreac(j))%P( 1) /= 0 ) template = TRIM(template)//"_gives_"//&
                                                   TRIM(speci(react(idx_dreac(j))%P( 1))%SpecyID)
             DO j1 = 2, 4
             IF ( react(idx_dreac(j))%P(j1) /= 0 ) template = TRIM(template)//"_"//&
                                                   TRIM(speci(react(idx_dreac(j))%P(j1))%SpecyID)
             ENDDO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'tot_frate_*    ', .false.) == 1 ).OR.&
                 ( INDEX(metatab(i)%IDname, 'tot_drate_*    ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          template = metatab(i)%IDname(1:idx_star-1)//TRIM(speci(speidx)%SpecyID)
          DO j = 1, nsubcols(i-ifil+1)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'compo_of_*     ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          DO j = 1, Nelements
             template = metatab(i)%IDname(1:idx_star-1)//TRIM(speci(speidx)%SpecyID)//'_'
             IF ( elements(j)%name == 'H      ') sp_lowcas = 'h      '
             IF ( elements(j)%name == 'D      ') sp_lowcas = 'd      '
             IF ( elements(j)%name == 'He     ') sp_lowcas = 'he     '
             IF ( elements(j)%name == 'C      ') sp_lowcas = 'c      '
             IF ( elements(j)%name == 'N      ') sp_lowcas = 'n      '
             IF ( elements(j)%name == 'O      ') sp_lowcas = 'o      '
             IF ( elements(j)%name == 'Na     ') sp_lowcas = 'na     '
             IF ( elements(j)%name == 'Mg     ') sp_lowcas = 'mg     '
             IF ( elements(j)%name == 'Si     ') sp_lowcas = 'si     '
             IF ( elements(j)%name == 'S      ') sp_lowcas = 's      '
             IF ( elements(j)%name == 'Fe     ') sp_lowcas = 'fe     '
             IF ( elements(j)%name == 'G      ') sp_lowcas = 'g      '
             template = TRIM(template)//TRIM(sp_lowcas)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE IF ( ( INDEX(metatab(i)%IDname, 'mass_of_*      ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%IDname, '*', .false.)
          template = metatab(i)%IDname(1:idx_star-1)//&
                     TRIM(speci(speidx)%SpecyID)//&
                     metatab(i)%IDname(idx_star+1:)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ELSE
          template = metatab(i)%IDname
          DO j = 1, nsubcols(i-ifil+1)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_IDname)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ENDIF
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 5 - Write the line containing all the dtypes in the output ASCII file
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       template = metatab(i)%dtype
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_dtype)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 6 - Write the line containing all the dtypes in the output ASCII file
    !     only radm_ini and radp_ini have units that depend on the other input parameters
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       template = metatab(i)%unit
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_unit)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 7 - Write the line containing all the SKOS in the output ASCII file
    !     SKOS need to be generated here only in a few cases
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       IF ( ( INDEX(metatab(i)%IDname, 'n_*           ', .false.) == 1 ).OR.&
            ( INDEX(metatab(i)%IDname, 'cd_prof_*     ', .false.) == 1 ).OR.&
            ( INDEX(metatab(i)%IDname, 'cd_*          ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%SKOS, '*', .false.)
          DO j = 1, Nspec
             template = metatab(i)%SKOS(1:idx_star-1)//&
                        TRIM(speci(j)%INCHI)//&
                        metatab(i)%SKOS(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_SKOS)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE
          template = metatab(i)%SKOS
          DO j = 1, nsubcols(i-ifil+1)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_SKOS)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ENDIF
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 8 - Write the line containing all the UCDs in the output ASCII file
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       template = metatab(i)%ucd
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_ucd)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 9 - Write the line containing all the utypes in the output ASCII file
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       template = metatab(i)%utype
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_utype)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 10 - Write the line containing all the groups in the output ASCII file
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       template = metatab(i)%group
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_group)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 11 - Write the line containing all the parents in the output ASCII file
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       template = metatab(i)%parent
       DO j = 1, nsubcols(i-ifil+1)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_parent)
          WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
       ENDDO
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 12 - Write the line containing all the descriptions in the output ASCII file
    !======================================================================================
    DO i = ifil, ifil + ndatafile - 1
       IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
       IF      ( ( INDEX(metatab(i)%IDname, 'elem_*        ', .false.) == 1 ) ) THEN
          idx_star = INDEX(metatab(i)%descrp, '*', .false.)
          DO j = 1, Nelements
             sp_uppcas = TRIM(ADJUSTL(elements(j)%name))
             template = metatab(i)%descrp(1:idx_star-1)//&
                        TRIM(sp_uppcas)//&
                        metatab(i)%descrp(idx_star+1:)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_descrp)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ELSE
          template = metatab(i)%descrp
          DO j = 1, nsubcols(i-ifil+1)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_descrp)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(template)
          ENDDO
       ENDIF
    ENDDO
    WRITE(iwrtmp,*)

    !======================================================================================
    ! 13 - Write the final lines containing the actual data in the output ASCII file
    !======================================================================================
    DO l = 1, nvals
       IF ( MOD(l-1, Nstep_hdf5) /= 0 ) CYCLE
       k = NINT(l * DBLE(counter) / DBLE(nvals))
       DO i = ifil, ifil + ndatafile - 1
          IF ( nsubcols(i-ifil+1) == 0 ) CYCLE
          key_temp = key_val
          IF (i == ifil) key_temp = ""
          IF      ( ( INDEX(metatab(i)%IDname, 'distance              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%distance
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'av                    ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Av
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'timeN                 ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%timeN
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'timeI                 ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%timeI
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'protdens              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nH
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rho_neut              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%rhoN
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rho_ions              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%rhoI
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rho_anio              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%rhoA
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rho_nega              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%rhoNeg
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'n_neut                ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%DensityN
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'n_ions                ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%DensityI
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'n_anio                ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%DensityA
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'n_nega                ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%DensityNeg
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mu_neut               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%muN
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mu_ions               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%muI
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mu_anio               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%muA
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mu_nega               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%muNeg
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'tgas                  ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Tn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'temp_neut             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Tn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'temp_ions             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Ti
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'temp_nega             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Te
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'iondegree             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%iondeg
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'ionfraction           ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%ionfrac
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'op_h2                 ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%oph2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'temp_grain            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%T_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'teff_grain            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Teff_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nlay_grain            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nlay_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'dens_grain            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%n_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mu_grain              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mu_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'much_grain            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%much_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'radius_grain          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%r_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mass_grain            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%m_gr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mass_core             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%m_grc
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mass_mantle           ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%m_grm
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'd_erosion             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%d_ero
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'd_adsorption          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%d_ads
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rsquare_core          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%rsq_grc
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rsquare_mantle        ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%rsq_grm
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'vsound                ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Vsound
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'vmagnet               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Vmagnet
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'valfven               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Valfven
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'velo_neut             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Vn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'velo_ions             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Vi
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'velo_nega             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%Vi
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'gradv_n               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%dVn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'gradv_i               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%dVi
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_tot          ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n            ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i            ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e            ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_molec      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_molec
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_H2         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_H2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_13CO       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_13CO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_OH         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_OH
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_NH3        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_NH3
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_rCO        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_rCO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_vCO        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_vCO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_CO         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_CO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_roH2O      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_roH2O
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_rpH2O      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_rpH2O
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_vH2O       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_vH2O
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_H2O        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_H2O
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_molec      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_molec
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_H2         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_H2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_molec      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_molec
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_H2         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_H2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_CO         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_CO
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_atoms      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_atoms
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Hat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Hat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Oat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Oat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Op         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Op
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Cat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Cat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Cp         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Cp
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Nat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Nat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Np         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Np
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Sat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Sat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Sp         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Sp
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Siat       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Siat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_n_Sip        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_n_Sip
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_atoms      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_atoms
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Hat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Hat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Oat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Oat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Op         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Op
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Cat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Cat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Cp         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Cp
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Nat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Nat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Np         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Np
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Sat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Sat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Sp         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Sp
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Siat       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Siat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_i_Sip        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_i_Sip
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_atoms      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_atoms
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Hat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Hat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Oat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Oat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Op         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Op
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Cat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Cat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Cp         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Cp
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Nat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Nat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Np         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Np
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Sat        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Sat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Sp         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Sp
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Siat       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Siat
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'line_rad_e_Sip        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%line_rad_e_Sip
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_tot           ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_n             ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_n
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_i             ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_i
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_e             ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_e
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_n_diss_H2     ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_n_diss_H2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_i_diss_H2     ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_i_diss_H2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_e_diss_H2     ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_e_diss_H2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_n_phgas       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_n_phgas
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_i_phgas       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_i_phgas
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_e_phgas       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_e_phgas
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_n_phgrn       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_n_phgrn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_i_phgrn       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_i_phgrn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_e_phgrn       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_e_phgrn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_n_cosmic      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_n_cosmic
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_i_cosmic      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_i_cosmic
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_e_cosmic      ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_e_cosmic
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_n_other       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_n_other
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_i_other       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_i_other
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chem_DE_e_other       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%chem_DE_e_other
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'elast_scat_tot        ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%elast_scat_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'elast_scat_n          ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%elast_scat_n
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'elast_scat_i          ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%elast_scat_i
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'elast_scat_e          ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%elast_scat_e
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'exch_eint_tot         ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%exch_eint_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'exch_eint_n           ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%exch_eint_n
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'exch_eint_i           ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%exch_eint_i
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'exch_eint_e           ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%exch_eint_e
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'therm_grain_tot       ', .false.) == 1 ) ) THEN 
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%therm_grain_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_tot         ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_n           ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_n
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_i           ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_i
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_e           ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_e
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_n_chem      ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_n_chem
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_i_chem      ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_i_chem
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_e_chem      ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_e_chem
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_n_compr     ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_n_compr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_i_compr     ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_i_compr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_e_compr     ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_e_compr
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_n_visc      ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_n_visc
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_i_visc      ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_i_visc
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mech_trsf_e_visc      ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mech_trsf_e_visc
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mom_flux_tot          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mom_flux_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mom_flux_kin          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mom_flux_kin
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mom_flux_the          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mom_flux_the
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mom_flux_mag          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mom_flux_mag
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mom_flux_vis          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%mom_flux_vis
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_tot          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_tot
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_kin          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_kin
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_the          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_the
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_mag          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_mag
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_vis          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_vis
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_int          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_int
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_src          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_src
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nrj_flux_pho          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(k)%nrj_flux_pho
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'n_*                   ', .false.) == 1 ) ) THEN
             DO j = 1, Nspec
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%abon(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'cd_prof_*             ', .false.) == 1 ) ) THEN
             DO j = 1,Nspec
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%coldens(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'pop_h2_*              ', .false.) == 1 ) ) THEN
             DO j = 1, NH2_lev
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%abon_h2lev(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h2_*              ', .false.) == 1 ) ) THEN
             DO j = 1, NH2_lines_out
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_h2lin(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_h_*               ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrhat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_h(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_c_*               ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrcat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_c(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_n_*               ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrnat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_n(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_o_*               ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntroat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_o(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_s_*               ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrsat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_s(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_si_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrsiat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_si(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_cp_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrcpl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_cp(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_np_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrnpl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_np(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_op_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntropl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_op(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sp_*              ', .false.) == 1 )) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrspl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_sp(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'emi_sip_*             ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrsipl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') traj_main(k)%emis_sip(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'cd_*                  ', .false.) == 1 ) ) THEN
             DO j = 1,Nspec
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') speci(j)%Col_dens
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'cd_h2_*               ', .false.) == 1 ) ) THEN
             DO j = 1, NH2_lev
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') H2_lev(j)%Col_dens
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_h2_*             ', .false.) == 1 ) ) THEN
             DO j = 1, NH2_lines_out
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') H2_lines(j)%intensity
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_h_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrhat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') inthat(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_c_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrcat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intcat(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_n_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrnat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intnat(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_o_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntroat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intoat(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_s_*              ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrsat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intsat(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_si_*             ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrsiat
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intsiat(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_cp_*             ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrcpl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intcpl(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_np_*             ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrnpl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intnpl(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_op_*             ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntropl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intopl(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_sp_*             ', .false.) == 1 )) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrspl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intspl(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'inta_sip_*            ', .false.) == 1 ) ) THEN
             idx_star = INDEX(metatab(i)%Hname, '*', .false.)
             DO j = 1, ntrsipl
                IF ( j /= 1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') intsipl(j)
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'uv_flux               ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') traj_main(k)%fluph
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'velo_drift            ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') ABS(traj_main(k)%Vn - traj_main(k)%Vi)
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'elecdens_gas          ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') traj_main(k)%abon(ind_e)
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'elecdens_grain        ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') traj_main(k)%abon(ind_Gminus) - traj_main(k)%abon(ind_Gplus) &
                                     + traj_main(k)%abon(ind_PAHm)   - traj_main(k)%abon(ind_PAHp)
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
         ELSE IF ( ( INDEX(metatab(i)%IDname, 'dens_spec_*           ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') traj_main(k)%abon(speidx)
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
         ELSE IF ( ( INDEX(metatab(i)%IDname, 'frate_*               ', .false.) == 1 ) ) THEN
            DO j = 1, nbf
               IF ( j /=1 ) key_temp = key_val
               WRITE(strf,'(1p,E16.8e3)') form_rate(j,l)
               WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
            ENDDO
         ELSE IF ( ( INDEX(metatab(i)%IDname, 'drate_*               ', .false.) == 1 ) ) THEN
            DO j = 1, nbd
               IF ( j /=1 ) key_temp = key_val
               WRITE(strf,'(1p,E16.8e3)') dest_rate(j,l)
               WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
            ENDDO
         ELSE IF ( ( INDEX(metatab(i)%IDname, 'tot_frate_*           ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') SUM(form_rate(1:nbf,l))
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
         ELSE IF ( ( INDEX(metatab(i)%IDname, 'tot_drate_*           ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') SUM(dest_rate(1:nbd,l))
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'species_file         ', .false.) == 1 ) ) THEN
             WRITE(strs,'(A)') TRIM(name_file_speci)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strs)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'chemistry_file        ', .false.) == 1 ) ) THEN
             WRITE(strs,'(A)') TRIM(name_file_chemistry)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strs)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'h2_file               ', .false.) == 1 ) ) THEN
             WRITE(strs,'(A)') TRIM(name_file_h2)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strs)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'shock_type_init       ', .false.) == 1 ) ) THEN
             WRITE(strs,'(A)') TRIM(shock_type_init)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strs)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'shock_type_fin        ', .false.) == 1 ) ) THEN
             WRITE(strs,'(A)') TRIM(shock_type_fin)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strs)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nfluids               ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') Nfluids
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'bbeta                 ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') Bbeta
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'vshock_km             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') Vs_km
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'delta_v_init          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') DeltaVmin
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'op_h2_init            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') op_H2_in
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'tn_init               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') traj_main(1)%Tn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nh_init               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') nH_init
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'tg_init               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') Tgrain
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'zeta                  ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') Zeta
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'f_isrf                ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') F_ISRF
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rad                   ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') RAD
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'av_init               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') Av0
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'f_coup_rad            ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') F_COUP_RAD
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'f_av                  ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') F_AV
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'f_invav               ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') F_invAv
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'conv_coldens          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') conv_coldens
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'n_h2_init             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') N_H2_0
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'n_co_init             ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') N_CO_0
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'vturb                 ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') vturb
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'f_tgr                 ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') F_TGR
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'amin_mrn              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') amin_mrn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'amax_mrn              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') amax_mrn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'alph_mrn              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') alph_mrn
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rho_grc               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') rho_grc
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'rho_grm               ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') rho_grm
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'nstep_max             ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') Nstep_max
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'xll_s_nh              ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') XLL
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'eps_v                 ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') Eps_V
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'timej                 ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') timeJ
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'duration_max          ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') duration_max
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'cool_kn               ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') Cool_KN
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'h_h2_flag             ', .false.) == 1 ) ) THEN
             WRITE(strs,'(A)') TRIM(H_H2_flag)
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strs)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'iforh2                ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') iforH2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'ikinh2                ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') ikinH2
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'elem_*                ', .false.) == 1 ) ) THEN
             DO j = 1, Nelements
                IF ( j /=1 ) key_temp = key_val
                WRITE(strf,'(1p,E16.8e3)') elements(j)%ab_init
                WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
             ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'compo_of_*             ', .false.) == 1 ) ) THEN
            DO j = 1, Nelements
               IF ( j /=1 ) key_temp = key_val
               WRITE(stri,'(I10)') speci(speidx)%formula(j)
               WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
            ENDDO
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'mass_of_*              ', .false.) == 1 ) ) THEN
            WRITE(strf,'(1p,E16.8e3)') speci(speidx)%mass / amu
            WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'stop_code             ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') stop_code
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'count_istat4          ', .false.) == 1 ) ) THEN
             WRITE(stri,'(I10)') count_istat4
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(stri)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'shock_size            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') size_shock
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ELSE IF ( ( INDEX(metatab(i)%IDname, 'shock_tsca            ', .false.) == 1 ) ) THEN
             WRITE(strf,'(1p,E16.8e3)') time_shock
             WRITE(iwrtmp,'(A)',ADVANCE='no') TRIM(key_temp)//TRIM(strf)
          ENDIF
       ENDDO
       WRITE(iwrtmp,*)
    ENDDO
 
    CLOSE(iwrtmp)
  
  END SUBROUTINE WRITE_METADATA


  SUBROUTINE FILL_RATES
    !---------------------------------------------------------------------------
    ! called by :
    !    WRITE_ASCII
    ! purpose :
    !    
    ! subroutine/function needed :
    ! input variables :
    ! ouput variables :
    ! results :
    !---------------------------------------------------------------------------
    USE MODULE_CHEM_REACT
    USE MODULE_PHYS_VAR,  ONLY : counter, Npthdf5
    USE MODULE_PROFIL_TABLES

    IMPLICIT none

    INTEGER                    :: npts
    INTEGER                    :: i, j, k, ir1, ir2, ir3, ip1, ip2, ip3, ip4                           !Tabone 3Body 12/2017
    INTEGER, DIMENSION(Nreact) :: multf, multd

    IF (counter <= Npthdf5) THEN
       npts = counter
    ELSE
       npts = Npthdf5
    ENDIF

    ! -------------------------------------------------------------------------------------
    ! Compute the number of reactions for the formation and destruction of species speidx
    ! -------------------------------------------------------------------------------------
    nbf   = 0
    nbd   = 0
    multf = 0
    multd = 0
    DO j = 1, Nreact
      ir1 = react(j)%R(1)
      ir2 = react(j)%R(2)
      ir3 = react(j)%R(3)                                                                        !Tabone 3Body 12/2017

      ip1 = react(j)%P(1)
      ip2 = react(j)%P(2)
      ip3 = react(j)%P(3)
      ip4 = react(j)%P(4)

      IF (ir1 == speidx) multd(j) = multd(j) + 1
      IF (ir2 == speidx) multd(j) = multd(j) + 1
      IF (ir3 == speidx) multd(j) = multd(j) + 1                                                    !Tabone 3Body 12/2017
      
      IF (ip1 == speidx) multf(j) = multf(j) + 1
      IF (ip2 == speidx) multf(j) = multf(j) + 1
      IF (ip3 == speidx) multf(j) = multf(j) + 1
      IF (ip4 == speidx) multf(j) = multf(j) + 1
      IF      (multf(j)-multd(j) > 0) THEN
         nbf = nbf + 1
      ELSE IF (multf(j)-multd(j) < 0) THEN
         nbd = nbd + 1
      ENDIF
    ENDDO

    ! -------------------------------------------------------------------------------------
    ! Allocate and fill the tables that contain the formation and destruction rates
    ! -------------------------------------------------------------------------------------
    ALLOCATE (form_rate(nbf,1:npts))
    ALLOCATE (dest_rate(nbd,1:npts))
    ALLOCATE (idx_freac(nbf))
    ALLOCATE (idx_dreac(nbd))

    nbf   = 0
    nbd   = 0
    DO j = 1, Nreact
       IF      (multf(j)-multd(j) > 0) THEN
          nbf = nbf + 1
          DO i = 1, npts
             k = NINT(i * DBLE(counter) / DBLE(npts))
             form_rate(nbf,i) = traj_main(k)%all_rate_chem(j) * ABS(multf(j)-multd(j))
          ENDDO
          idx_freac(nbf) = j
       ELSE IF (multf(j)-multd(j) < 0) THEN
          nbd = nbd + 1
          DO i = 1, npts
             k = NINT(i * DBLE(counter) / DBLE(npts))
             dest_rate(nbd,i) = traj_main(k)%all_rate_chem(j) * ABS(multf(j)-multd(j))
          ENDDO
          idx_dreac(nbd) = j
       ENDIF
    ENDDO

  END SUBROUTINE FILL_RATES

END MODULE OUTPUT_HDF5
