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

MODULE MODULE_GAMMA_NUMREC

   IMPLICIT NONE
   INCLUDE "precision.f90"

   ! variables from "precision.f90" are private
   PRIVATE :: DP, LONG, minus_infinity, plus_infinity



CONTAINS



   FUNCTION gammln(xx)
      !---------------------------------------------------------------------------
      ! called by :
      !     gser
      !     gcf
      ! purpose :
      !     Returns the value ln [Γ(xx)] for xx > 0.
      !     See Numerical recipes p 207
      ! subroutine/function needed :
      ! input variables :
      ! output variables :
      ! results :
      !     gammln -> ln [Γ(xx)]
      !---------------------------------------------------------------------------
      IMPLICIT none

      REAL(dp), PARAMETER     :: stp = 2.5066282746310005_dp
      REAL(dp)                :: gammln
      REAL(dp)                :: xx
      INTEGER                 :: j
      REAL(dp)                :: ser
      REAL(dp)                :: tmp
      REAL(dp)                :: x
      REAL(dp)                :: y
      REAL(dp), DIMENSION (6) :: cof
      DATA cof                 / 76.18009172947146e+0_dp,&
                                -86.50532032941677e+0_dp,&
                                 24.01409824083091e+0_dp,&
                                -1.231739572450155e+0_dp,&
                                 1.208650973866179e-3_dp,&
                                -5.395239384953000e-6_dp /

      x = xx
      y = x
      tmp = x + 5.5_dp
      tmp = ( x + 0.5_dp ) * LOG(tmp) - tmp
      ser = 1.000000000190015_dp
      DO j = 1, 6
         y   = y + 1.0_dp
         ser = ser + cof(j) / y
      ENDDO
      gammln = tmp + LOG( stp * ser / x )

   END FUNCTION gammln



   FUNCTION gammp(a,x)
      !---------------------------------------------------------------------------
      ! called by :
      !     SIMPLE_TRANSFER
      ! purpose :
      !     Returns the incomplete gamma function P(a,x).
      !     See Numerical recipes p 209-213
      ! subroutine/function needed :
      !     gser
      !     gcf
      ! input variables :
      ! output variables :
      ! results :
      !     gammp -> incomplete gamma function P
      !---------------------------------------------------------------------------
      IMPLICIT none

      REAL(dp) :: a
      REAL(dp) :: gammp
      REAL(dp) :: x
      REAL(dp) :: gammcf
      REAL(dp) :: gamser
      REAL(dp) :: gln

      IF ( x < 0.0_dp .OR. a < 0.0_dp ) THEN
         WRITE (*,*) "bad arguments in gammp"
         STOP
      ENDIF
      ! ---------------------------------------------
      ! use the series representation ... or ...
      ! use the continued fraction representation
      ! ---------------------------------------------
      IF (x < a+1.0_dp ) THEN 
         CALL gser(gamser,a,x,gln)
         gammp = gamser
      ELSE
         CALL gcf(gammcf,a,x,gln)
         gammp = 1.0_dp - gammcf
      ENDIF
 
   END FUNCTION gammp



   FUNCTION gammq(a,x)
      !---------------------------------------------------------------------------
      ! called by :
      !     nobody
      ! purpose :
      !     Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
      !     See Numerical recipes p 209-213
      ! subroutine/function needed :
      !     gser
      !     gcf
      ! input variables :
      ! output variables :
      ! results :
      !     gammq -> incomplete gamma function Q
      !---------------------------------------------------------------------------
      IMPLICIT none

      REAL(dp) :: a
      REAL(dp) :: gammq
      REAL(dp) :: x
      REAL(dp) :: gammcf
      REAL(dp) :: gamser
      REAL(dp) :: gln

      IF ( x < 0.0_dp .OR. a < 0.0_dp ) THEN
         WRITE (*,*) "bad arguments in gammq"
         STOP
      ENDIF
      ! ---------------------------------------------
      ! use the series representation ... or ...
      ! use the continued fraction representation
      ! ---------------------------------------------
      !IF (x < a+1.0_dp ) THEN 
         CALL gser(gamser,a,x,gln)
         gammq = 1.0_dp - gamser
      !ELSE
      !   CALL gcf(gammcf,a,x,gln)
      !   gammq = gammcf
      !ENDIF

   END FUNCTION gammq(a,x)



   SUBROUTINE gser(gamser,a,x,gln)
      !---------------------------------------------------------------------------
      ! called by :
      !     gammp
      !     gammq
      ! purpose :
      !     Returns the incomplete gamma function P(a,x) evaluated by its series
      !     representation as gamser. Also returns ln Γ(a) as gln.
      !     See Numerical recipes p 209-213
      ! subroutine/function needed :
      !     gammln
      ! input variables :
      ! output variables :
      ! results :
      !     gamser -> incomplete gamma function P(a,x)
      !---------------------------------------------------------------------------
      IMPLICIT none

      INTEGER,  PARAMETER :: ITMAX = 100
      REAL(dp), PARAMETER :: EPS   = 3.e-07_dp
      REAL(dp)            :: a
      REAL(dp)            :: gamser
      REAL(dp)            :: gln
      REAL(dp)            :: x
   
      INTEGER             :: n
      REAL(dp)            :: ap
      REAL(dp)            :: del
      REAL(dp)            :: sum
      REAL(dp)            :: gammln
      LOGICAL             :: success

      success = .FALSE.
      gln = gammln(a)
      IF ( x <= 0.0_dp ) THEN
         IF (x < 0.0_dp ) THEN
            WRITE (*,*) "bad arguments in gser"
            STOP
         ENDIF
         gamser = 0.0_dp
         RETURN
      ENDIF
      ap  = a
      sum = 1.0_dp / a
      del = sum
      DO n = 1, ITMAX
         ap  = ap + 1.0_dp
         del = del * x / ap
         sum = sum + del
         IF ( abs(del) < abs(sum) * EPS ) THEN
            success = .TRUE.
            EXIT
         ENDIF
      ENDDO
      gamser = sum * EXP(-x+a*log(x)-gln)
      IF ( .NOT.success ) THEN
         WRITE(*,*) "Warning ... a too large, ITMAX too small in gser"
      ENDIF

   END SUBROUTINE gser



   SUBROUTINE gcf(gammcf,a,x,gln)
      !---------------------------------------------------------------------------
      ! called by :
      !     gammp
      !     gammq
      ! purpose :
      !     Returns the incomplete gamma function Q(a,x) evaluated by its continued 
      !     fraction representation as gammcf. Also returns ln Γ(a) as gln.
      !     Parameters: ITMAX is the maximum allowed number of iterations
      !                 EPS is the relative accuracy
      !                 FPMIN is a number near the smallest representable floating-point number.
      !                 Set up for evaluating continued fraction by modified
      !                 Lentz’s method (§5.2) with b0 = 0.
      !                 Iterate to convergence.
      !     See Numerical recipes p 209-213
      ! subroutine/function needed :
      !     gammln
      ! input variables :
      ! output variables :
      ! results :
      !     gammcf -> incomplete gamma function Q(a,x)
      !---------------------------------------------------------------------------
      IMPLICIT none

      INTEGER,  PARAMETER :: ITMAX = 100
      REAL(dp), PARAMETER :: EPS   = 3.e-07_dp
      REAL(dp), PARAMETER :: FPMIN = 1.e-30_dp
      REAL(dp)            :: a
      REAL(dp)            :: gammcf
      REAL(dp)            :: gln
      REAL(dp)            :: x
 
      INTEGER             :: i
      INTEGER             :: an
      REAL(dp)            :: b
      REAL(dp)            :: c
      REAL(dp)            :: d
      REAL(dp)            :: del
      REAL(dp)            :: h
      REAL(dp)            :: gammln
      LOGICAL             :: success

      success = .FALSE.
      gln = gammln(a)
      b   = x + 1.0_dp - a
      c   = 1.0_dp / FPMIN
      d   = 1.0_dp / b
      h   = d
      DO i = 1, ITMAX
         an = -i * ( i - a )
         b  = b + 2.0_dp
         d  = an * d + b
         IF ( abs(d) < FPMIN ) THEN
            d=FPMIN
         ENDIF
         c = b + an / c
         IF ( abs(c) < FPMIN ) THEN
            c = FPMIN
         ENDIF
         d   = 1.0_dp / d
         del = d * c
         h   = h * del
         IF ( abs( del - 1.) < EPS ) THEN
            success = .TRUE.
            EXIT
         ENDIF
      ENDDO
      gammcf = exp ( -x + a * log(x) - gln ) * h
      IF ( .NOT.success ) THEN
         WRITE(*,*) "Warning ... a too large, ITMAX too small in gcf"
      ENDIF

   END SUBROUTINE gcf

END MODULE MODULE_GAMMA_NUMREC