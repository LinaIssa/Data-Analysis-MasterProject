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

MODULE MODULE_SHIELD
  USE MODULE_TOOLS, ONLY : GET_FILE_NUMBER
  IMPLICIT none
  INCLUDE "precision.f90"
  private read_until,interpolate



CONTAINS



  REAL(KIND=DP) FUNCTION Av_guess(NH2)
     USE MODULE_PHYS_VAR, ONLY : F_AV, inv_Av_fac, Av0, N_H2_0
     REAL(KIND=DP) :: NH2

     IF (F_AV  == 0) THEN
        Av_guess = Av0
     ELSE
        Av_guess = Av0+(NH2-N_H2_0)*2_DP*abs(inv_Av_fac)
        Av_guess = MAX(0.0_DP,Av_guess)
     ENDIF
     RETURN
  END FUNCTION Av_guess


  REAL(KIND=DP) FUNCTION shield_H2(Av, NH2, key)
     REAL(KIND=DP) :: NH2, Av
     CHARACTER(5)  :: key

     SELECT CASE(key)
     CASE('lee  ')
        shield_H2 = shield_H2_lee(NH2)
     CASE('D&B  ')
        shield_H2 = shield_H2_DB(Av,NH2)
     CASE DEFAULT 
        shield_H2 = 1_DP
     END SELECT
  END FUNCTION shield_H2


  REAL(KIND=DP) FUNCTION shield_CO(Av,NH2,NCO,key)
     USE MODULE_TOOLS, ONLY : MINIMUM, MAXIMUM
     REAL(KIND=DP) :: Av, NH2, NCO
     CHARACTER(5)  :: key
     ! REAL(KIND=DP) :: theta1,theta2,theta3

     ! -----------------------------------------------------
     ! test - use analytical function instead of the table
     !        of Lee et al. (1996) -> in order to help dvode
     !        comment one of the block or the other
     ! -----------------------------------------------------
     ! theta1 = 6.0e-1_dp*exp(-8.0_dp*NCO/1e12_dp) &
     !        + 4.0e-1_dp*exp(-2.0_dp*NCO/1e16_dp) &
     !        + 2.0e-3_dp*exp(-2.0_dp*NCO/1e19_dp) &
     !        + 2.0e-1_dp*exp(-(log(NCO)-log(1e11_dp))**2.0_dp/3.0_dp**2.0_dp) &
     !        + 6.0e-1_dp*exp(-(log(NCO)-log(1e13_dp))**2.0_dp/5.0_dp**2.0_dp)
     ! theta2 = 1.0e-1_dp*exp(-1.0_dp*NH2/1e07_dp) &
     !        + 8.0e-1_dp*exp(-6.0_dp*NH2/1e22_dp) &
     !        + 1.2e-1_dp*exp(-1.0_dp*NH2/1e19_dp) &
     !       + 8.0e-2_dp*exp(-(log(NH2)-log(1e12_dp))**2.0_dp/13.0_dp**2.0_dp)
     ! theta3 = 1.0e+0_dp*exp(-6.0_dp*Av) &
     !       + 4.0e-2_dp*exp(-(log(Av)-log(0.55_dp))**2/0.58_dp**2.0_dp)

     ! theta1 = MINIMUM(theta1,1.0_dp)
     ! theta2 = MINIMUM(theta2,1.0_dp)
     ! theta3 = MINIMUM(theta3,1.0_dp)

     ! shield_CO = theta1*theta2*theta3
     ! -----------------------------------------------------
     SELECT CASE(key)
     CASE('lee  ')
        shield_CO=shield_CO_lee(Av,NH2,NCO)
     CASE('D&B  ')
        PRINT *, 'Draine&Bertoldi did not compute CO shielding...'
        STOP
     CASE DEFAULT
        shield_CO = 1_DP
     END SELECT
     ! -----------------------------------------------------
  END FUNCTION shield_CO


  REAL(KIND=DP) FUNCTION shield_H2_DB(Av,N)
     USE module_phys_var,  ONLY : Tn
     USE module_constants, ONLY : mp, kb
     IMPLICIT none
     REAL(KIND=DP), PARAMETER :: vturb=1e5
     REAL(KIND=DP)            :: Av, N
     REAL(KIND=DP)            :: b, x

     ! Draine&Bertoldi, ApJ 1996, equation (37)
     x = N / (5.e14_DP)
     b = sqrt(kb * Tn /mp + vturb**2) / 1e5 ! Doppler parameter in km/s
     shield_H2_DB = .965_DP / (1._DP + x / b)**2 &
                  + .035_DP / sqrt(1._DP + x) * exp(-8.5e-4_DP * sqrt(1._DP + x))

     ! Draine&Bertoldi, ApJ 1996, equation (40)
     ! interpreted by Guillaume in terms of Av
     shield_H2_DB = shield_H2_DB * exp(-6.3_DP * Av)
  END FUNCTION shield_H2_DB


  REAL(KIND=DP) FUNCTION shield_H2_lee(NH2)
     USE MODULE_TECHCONFIG
     CHARACTER(len=lenfilename)                     :: name_file_shield = 'shield_H2.dat'
     REAL(KIND=DP), SAVE, ALLOCATABLE, DIMENSION(:) :: xNH2, fNH2
     REAL(KIND=DP)                                  :: NH2, lNH2
     REAL(KIND=dp)                                  :: ofNH2
     LOGICAL, SAVE                                  :: first = .true.
     INTEGER, SAVE                                  :: iNH2 = 1 ! Positions in table from previous calls
     INTEGER                                        :: file_shield

     ! Read data at first call
     IF (first) THEN
        first = .false.
        file_shield=GET_FILE_NUMBER()
        name_file_shield = TRIM(data_dir) // TRIM(name_file_shield)
        OPEN(file_shield, file=name_file_shield, status='OLD')
        ! Read until hit keyword END
        CALL read_until('END', file_shield, xNH2, fNH2)
        CLOSE(file_shield)
        ! Some quantities are better behaved when logged
        xNH2 = log(xNH2)
        fNH2 = log(fNH2)
     ENDIF
      
     ! Convert input values to log, with security for small col_dens.
     IF (NH2.le.0_DP) THEN
        lNH2 = xNH2(1)
     ELSE
        lNH2 = max(log(NH2),xNH2(1))
     ENDIF

     ! Interpolate and return shielding value
     ofNH2 = interpolate(lNH2,iNH2,xNH2,fNH2)
     shield_H2_lee = exp(ofNH2)
  END FUNCTION shield_H2_lee


  REAL(KIND=DP) FUNCTION shield_CO_lee(Av,NH2,NCO)
     USE MODULE_TECHCONFIG
     CHARACTER(len=lenfilename)                     :: name_file_shield='shield_CO.dat'
     REAL(KIND=DP), SAVE, ALLOCATABLE, DIMENSION(:) :: xAv, fAv, xNH2, fNH2, xNCO, fNCO
     REAL(KIND=DP)                                  :: Av, NH2, NCO, lNH2, lNCO
     REAL(KIND=DP)                                  :: ofAv,ofNH2,ofNCO
     LOGICAL, SAVE                                  :: first=.true.
     INTEGER, SAVE                                  :: iAv = 1, iNH2 = 1, iNCO = 1 ! Positions in table from previous calls
     INTEGER                                        :: file_shield

     ! Read data at first call
     IF (first) THEN
        first = .false.
        file_shield = GET_FILE_NUMBER()
        name_file_shield = TRIM(data_dir) // TRIM(name_file_shield)
        OPEN(file_shield,file=name_file_shield,status='OLD')
        ! Read until hit keyword END
        CALL read_until('END',file_shield,xNCO,fNCO)
        CALL read_until('END',file_shield,xNH2,fNH2)
        CALL read_until('END',file_shield,xAv ,fAv)
        CLOSE(file_shield)
        ! Some quantities are better behaved when logged
        xNH2 = log(xNH2)
        xNCO = log(xNCO)
        fAv  = log(fAv)
     ENDIF

     ! Convert input values to log, with security for small col_dens.
     IF (NH2.le.0_DP) THEN
        lNH2 = xNH2(1)
     ELSE
        lNH2 = max(log(NH2),xNH2(1))
     ENDIF
     IF (NCO.le.0_DP) THEN
        lNCO = xNCO(1)
     ELSE
        lNCO = max(log(NCO),xNCO(1))
     ENDIF

     ! Interpolate 
     ofAv  = interpolate(Av,iAv,xAv,fAv)
     ofNH2 = interpolate(lNH2,iNH2,xNH2,fNH2)
     ofNCO = interpolate(lNCO,iNCO,xNCO,fNCO)

     ! Return shielding value
     shield_CO_lee = exp(ofAv) * ofNH2 * ofNCO
  END FUNCTION shield_CO_lee


  ! Read within a file a 2 column table until keyword is hit
  SUBROUTINE read_until(keyword,file_handle,x1,x2)
     INTEGER,       PARAMETER                 :: nmax = 1000
     REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: x1, x2
     REAL(KIND=DP), DIMENSION(nmax)           :: v1, v2
     CHARACTER(3)                             :: keyword
     CHARACTER(len=80)                        :: one_line
     INTEGER                                  :: file_handle
     INTEGER                                  :: size
     LOGICAL                                  :: ok

     size=0
     ok=.true.
     DO WHILE(ok)
        READ(file_handle,'(A80)')one_line
        IF (one_line(1:len(keyword))==keyword) THEN
           ! Break loop when keyword is hit
           ok=.false.
        ELSE
           IF (one_line(1:1).ne.'!') THEN ! Skip comments
              size=size+1
              READ(one_line,*)v1(size),v2(size)
           ENDIF
        ENDIF
     ENDDO
     ! allocate arrays and return values
     ALLOCATE(x1(size),x2(size))
     x1 = v1(:size)
     x2 = v2(:size)
  END SUBROUTINE read_until


  REAL(KIND=DP) FUNCTION interpolate(x,ix,xa,ya)
    REAL(KIND=DP)               :: x
    REAL(KIND=DP), DIMENSION(:) :: xa, ya
    INTEGER                     :: ix ! first guess for the position in the table
    logical                     :: ok
    integer                     :: n

    n = size(xa)
    ! Check bounds
    ! Smaller than lower bound is not allowed
    IF (x<xa(1)) THEN
       PRINT *, 'Input value too low... There must be a problem:'
       PRINT *, 'x=', x, ' Lower bound of array=', xa(1)
       STOP
    ENDIF
    ! Greater than last bound means we output the last value in the table
    IF (x.ge.xa(n)) THEN
       ix = n-1
       interpolate = ya(n)
       RETURN
    ENDIF
    ! Find position in table. ix is the first guess from a previous point, usually.
    ok=(x.ge.xa(ix)).and.(x<xa(ix+1))
    DO WHILE (.not.ok) 
       IF (x<xa(ix)) THEN
          ix = ix-1
       ELSE
          ix = ix+1
       ENDIF
       ok = (x.ge.xa(ix)).and.(x<xa(ix+1))
    ENDDO
    ! Finally interpolate
    interpolate = (ya(ix)*(xa(ix+1)-x)+ya(ix+1)*(x-xa(ix)))/(xa(ix+1)-xa(ix))
  END FUNCTION interpolate

END MODULE MODULE_SHIELD
