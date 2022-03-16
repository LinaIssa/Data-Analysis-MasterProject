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

MODULE MODULE_TECHCONFIG
  !***************************************************************************************
  !** This module contains all technical variable to produce HDF5 output files          **
  !***************************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! ======================================================================================
  ! List of flags
  ! - default values in parenthesis
  ! ======================================================================================
   INTEGER, PUBLIC            :: F_W_ASCII          ! (1) Write ascii files
   INTEGER, PUBLIC            :: F_W_HDF5_STD       ! (1) Write standard  HDF5 file
   INTEGER, PUBLIC            :: F_W_HDF5_CHE       ! (1) Write chemistry HDF5 file

  ! ======================================================================================
  ! Size of some character strings
  ! ======================================================================================
  INTEGER, PUBLIC, PARAMETER :: lenfilename   = 100 ! Max size of the name of a file
  INTEGER, PUBLIC, PARAMETER :: lenshcmnd     = 500 ! Max length of shell command called with SYSTEM

  INTEGER, PUBLIC, PARAMETER :: lenmetaname   = 100 ! Max size of the path / file / name / ID of any data (see METADATA TYPE)
  INTEGER, PUBLIC, PARAMETER :: lenmetatype   =  20 ! Max size of the type / unit / group / of any data (see METADATA TYPE)
  INTEGER, PUBLIC, PARAMETER :: lenmetadesc   = 500 ! Max size of the description of any data (see METADATA TYPE)
  INTEGER, PUBLIC, PARAMETER :: lenmetakey    =   8 ! Max size of the description of any data (see METADATA TYPE)

  INTEGER, PUBLIC, PARAMETER :: lenstrint     =  12 ! Max length of a string from an integer casted to string by INT_STR

  INTEGER, PUBLIC, PARAMETER :: lenIDLEV =  30      ! Size of a string to code quantum numbers identifications
  INTEGER, PUBLIC, PARAMETER :: lenIDLIN = 100      ! Size of a string to code line transitions

  ! ======================================================================================
  ! Input / output files & directories
  ! ======================================================================================
  CHARACTER (LEN=lenfilename), PUBLIC            :: modele
  CHARACTER (LEN=lenfilename), PUBLIC            :: specfile
  CHARACTER (LEN=lenfilename), PUBLIC            :: chemfile
  CHARACTER (LEN=lenfilename), PUBLIC            :: h2exfile
  CHARACTER (LEN=lenfilename), PUBLIC            :: gridfile
  CHARACTER (LEN=lenfilename), PUBLIC            :: inpfile
  CHARACTER (LEN=6),           PUBLIC, PARAMETER :: data_dir        = 'input/'
  CHARACTER (LEN=7),           PUBLIC, PARAMETER :: grains_dir      = 'Grains/'
  CHARACTER (LEN=15),          PUBLIC, PARAMETER :: phddat_dir      = 'cross_sections/'
  CHARACTER (LEN=11),          PUBLIC, PARAMETER :: metdat_dir      = 'TechConfig/'
  CHARACTER (LEN=17),          PUBLIC, PARAMETER :: metstd_fil      = 'meta_standard.dat'
  CHARACTER (LEN=17),          PUBLIC, PARAMETER :: metche_fil      = 'meta_chemistr.dat'
  CHARACTER (LEN=19),          PUBLIC, PARAMETER :: write_hdf5_file = 'shock_write_HDF5.py'
  CHARACTER (LEN=lenfilename), PUBLIC            :: out_dir         = 'output/'
  CHARACTER (LEN=16),          PUBLIC            :: ascii_std       = 'ASCII_files_dat/'
  CHARACTER (LEN=16),          PUBLIC            :: ascii_che       = 'ASCII_files_che/'
  CHARACTER (LEN=16),          PUBLIC            :: ascii_ana       = 'ASCII_files_ana/'
  CHARACTER (LEN=lenfilename), PUBLIC            :: ascii_dir_std ! Directory where standard ASCII files are stored
  CHARACTER (LEN=lenfilename), PUBLIC            :: ascii_dir_che ! Directory where chemical ASCII files are stored
  CHARACTER (LEN=lenfilename), PUBLIC            :: ascii_dir_ana ! Directory where chemical ASCII files are stored

  CHARACTER (LEN=255),         PUBLIC            :: fichier
  CHARACTER (LEN=9),           PUBLIC, PARAMETER :: python_version  = 'python2.7'
  INTEGER,                     PUBLIC            :: python_exists

  ! ======================================================================================
  ! Shell command
  ! ======================================================================================
  CHARACTER(LEN=lenshcmnd),    PUBLIC            :: sh_cmd     ! Shell command

  TYPE METADATA !------------------------------------------------------------------------!
    ! This type contains all necessary information (metadata) to completely describe the !
    ! outputs of the shock code. Metadata are used in the HDF5 file and the subsequent   !
    ! applications (extractor, online database services, ...)                            !
    ! -----------------------------------------------------------------------------------!
    CHARACTER(len=lenmetaname)                  :: path   ! path towards dataset         !
    CHARACTER(len=lenmetaname)                  :: file   ! ASCII file to store data     !
    CHARACTER(len=lenmetaname)                  :: datset ! name of the dataset          !
    CHARACTER(len=lenmetaname)                  :: Hname  ! human readable               !
    CHARACTER(len=lenmetaname)                  :: IDname ! param ID must be unique      !
    CHARACTER(len=lenmetatype)                  :: dtype  ! data type: real, int, string !
    CHARACTER(len=lenmetatype)                  :: unit   ! unit of the data             !
    CHARACTER(len=lenmetadesc)                  :: SKOS   ! documentary language         !
    CHARACTER(len=lenmetatype)                  :: ucd    ! UCD of the data              !
    CHARACTER(len=lenmetaname)                  :: utype  ! SimDM quantity... unclear    !
    CHARACTER(len=lenmetatype)                  :: group  ! gas, grain, ...              !
    CHARACTER(len=lenmetatype)                  :: parent ! parent object: meshcell, ... !
    CHARACTER(len=lenmetadesc)                  :: descrp ! human readable descr of data !
  END TYPE METADATA !--------------------------------------------------------------------!

   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_path   = TRIM('@path='  )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_datset = TRIM('@dset='  )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_Hname  = TRIM('@Name='  )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_IDname = TRIM('@PubID=' )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_dtype  = TRIM('@dtype=' )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_SKOS   = TRIM('@SKOS='  )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_unit   = TRIM('@unit='  )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_ucd    = TRIM('@UCD='   )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_utype  = TRIM('@utype=' )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_group  = TRIM('@Group=' )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_parent = TRIM('@ObjPar=')
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_descrp = TRIM('@Desc='  )
   CHARACTER(LEN=lenmetakey), PUBLIC, PARAMETER       :: key_val    = TRIM('val='    )

  ! ======================================================================================
  ! TechConfig files
  ! ======================================================================================
  INTEGER,          PUBLIC                            :: nhdf5_std   ! Number of output to write in the standard HDF5 file (number of lines in the file meta_standard.dat)
  INTEGER,          PUBLIC                            :: nhdf5_che   ! Number of output to write in the chemical HDF5 file (number of lines in the file meta_chemistr.dat)
  TYPE (METADATA),  PUBLIC, ALLOCATABLE, DIMENSION(:) :: metatab_std ! Store all the information written in the file metstd_fil
  TYPE (METADATA),  PUBLIC, ALLOCATABLE, DIMENSION(:) :: metatab_che ! Store all the information written in the file metche_fil

  ! ======================================================================================
  ! variables from "precision.f90" are private
  ! ======================================================================================
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

END MODULE MODULE_TECHCONFIG
