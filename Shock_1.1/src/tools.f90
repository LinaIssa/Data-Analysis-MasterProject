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

! Tabulated exponentials and powers
MODULE tex
  IMPLICIT NONE
  INCLUDE "precision.f90"
  ! Tabulated exponential and power
  integer,parameter::n=10001
  real(kind=DP),parameter::xmin=-50._DP,xmax=50._DP
  real(kind=DP)::exp_max
  real(kind=DP),dimension(n)::table_exp
  real(kind=DP)::pow_max
  real(kind=DP),dimension(n)::table_pow
  private
  public texp,tpow
contains

! Tabulated exponential
  real(kind=dp) function texp(x)
    real(kind=dp)::x
    real(kind=dp)::xi
    logical,save::first=.true.
    integer::i
! Use system exp
    texp=exp(x)
    return
! Skip the table
    if (first) then
       first=.false.
       exp_max=exp(xmax)
       ! build table
       do i=1,n
          table_exp(i)=exp(dble(i-1)*(xmax-xmin)/dble(n-1)+xmin)
       enddo
    endif

    
    if (x<xmin) then
       texp=0_dp
       return
    else if(x>=xmax) then
       texp=exp_max
       return
    else
       xi=(x-xmin)/(xmax-xmin)*dble(n-1)+1._DP
       i=int(xi)
       xi=xi-dble(i)
       texp=table_exp(i)*(1._DP-xi)+xi*table_exp(i+1)
       return
    endif
  end function texp
  
! Tabulated powers of 10
  real(kind=dp) function tpow(x)
    real(kind=dp)::x
    real(kind=dp)::xi
    logical,save::first=.true.
    integer::i
! Use system power
    tpow=10._DP**(x)
    return
! Skip the table
    if (first) then
       first=.false.
       pow_max=10._DP**(xmax)
       ! build table
       do i=1,n
          table_pow(i)=10._DP**(dble(i-1)*(xmax-xmin)/dble(n-1)+xmin)
       enddo
    endif

    
    if (x<xmin) then
       tpow=0_dp
       return
    else if(x>=xmax) then
       tpow=pow_max
       return
    else
       xi=(x-xmin)/(xmax-xmin)*dble(n-1)+1._DP
       i=int(xi)
       xi=xi-dble(i)
       tpow=table_pow(i)*(1._DP-xi)+xi*table_pow(i+1)
       return
    endif
  end function tpow  
end MODULE tex


MODULE MODULE_TOOLS
  !*****************************************************************************
  !** The module 'MODULE_TOOL' contains useful non-numerical procedures       **
  !** used by severeal subroutine of the MHD shock code                       **
  !*****************************************************************************
  IMPLICIT NONE
  INCLUDE "precision.f90"

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity


CONTAINS


  FUNCTION GET_FILE_NUMBER() RESULT (num)
    !---------------------------------------------------------------------------
    ! purpose :
    !     gives a file number available for opening, i.e. a free logical unit
    !     example :
    !         n=GET_FILE_NUMBER()
    !         OPEN(n,file='sample.txt')
    ! subroutine/function needed :
    ! input variables :
    ! output variables :
    ! results :
    !     num -> 'LONG' integer : the first free logical unit
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG) :: num
    LOGICAL :: opened

    opened=.TRUE.
    num=9 ! start file number after 10
    DO WHILE (opened .AND. num < 50)
       num=num+1
       INQUIRE(UNIT=num, opened=opened)
    END DO
    IF (num == 50 .and. opened) STOP "*** WARNING : no file number available ***"
  END FUNCTION GET_FILE_NUMBER


  FUNCTION INTEGRATE_SCALAR(dist_step, value1, value2) RESULT(res)
    !---------------------------------------------------------------------------
    ! purpose :
    !     Integrates the function funct using trapezian rule, the independant
    !     variable is the distance (dist_step =distance-distance_old).
    !     Useful to calculate column densities or flow times.
    ! subroutine/function needed :
    ! input variables :
    !     funct          -> 'function' to integrate (scalar)
    !     dist_step  -> = distance2-distance1 (cm)
    !     value1, value2 -> values to integrate at distance1 and distance2
    ! output variables :
    ! results :
    !     res = 0.5_DP*dist_step*(value1+value2)
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(in) :: dist_step, value1, value2
    REAL(KIND=DP) :: res

    res = 0.5_DP*dist_step*(value1+value2)

  END FUNCTION INTEGRATE_SCALAR



  SUBROUTINE SORT_INCREASING(dim,vect1,vect2,vect3,vect4)
    !---------------------------------------------------------------------------
    ! called by :
    ! purpose :
    !     sorts 4 vectors in the increasing order of the first one
    !     method=shell
    ! subroutine/function needed :
    ! input variables :
    !     dim -> dimension of the 4 vectors, or number of elements to sort
    !     vect1, vect2, vect3, vect4 : vectors to sort
    ! output variables :
    !     vect1, vect2, vect3, vect4 : sorted vectors
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER(KIND=LONG), INTENT(in) :: dim
    REAL(KIND=DP),DIMENSION(:),INTENT(inout) :: vect1,vect2,vect3,vect4
    REAL(KIND=DP) :: v1,v2,v3,v4
    INTEGER(KIND=LONG) :: i,j,inc

    ! dim has to be <= dimension of the vectors
    IF (dim > SIZE(vect1)) STOP "*** WARNING, the vector is to small"
    ! the four vectors must have the same dimension
    IF ( SIZE(vect1)==SIZE(vect2) .AND. &
         SIZE(vect2)==SIZE(vect3) .AND. &
         SIZE(vect3)==SIZE(vect4)) THEN
       ! initial value of the increment
       inc=1
       DO WHILE (inc <= dim)
          inc=3*inc+1
       END DO

       ! sub-vectors distants by inc are sorted first, and then inc is decreased
       DO WHILE (inc > 1)
          inc=inc/3
          DO i=inc+1,dim
             v1=vect1(i)
             v2=vect2(i)
             v3=vect3(i)
             v4=vect4(i)
             j=i
             DO WHILE ((vect1(j-inc) > v1) .AND. (j > inc))
                vect1(j)=vect1(j-inc)
                vect2(j)=vect2(j-inc)
                vect3(j)=vect3(j-inc)
                vect4(j)=vect4(j-inc)
                j=j-inc
             END DO
             vect1(j)=v1
             vect2(j)=v2
             vect3(j)=v3
             vect4(j)=v4
          END DO
       END DO
    ELSE
       STOP "***, WARNING, the four vectors must have the same dimension in SORT_INCREASING"
    END IF

  END SUBROUTINE SORT_INCREASING


  REAL(KIND=DP) FUNCTION MINIMUM(a,b)
    !---------------------------------------------------------------------------
    ! called by :
    ! purpose :
    !     gives an approximation for the min(a,b) so that the resulting function
    !     is infinitely derivable and continuous with respect to a, b
    ! subroutine/function needed :
    ! input variables :
    !     a -> real number
    !     b -> real number
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(in) :: a
    REAL(KIND=DP), INTENT(in) :: b

    IF      ( ( a > 0.0_dp ).AND.( b > 0.0_dp ) ) THEN
       MINIMUM = 1.0_dp / ( 1.0_dp/a**8.0_dp + 1.0_dp/b**8_dp )**(1.0_dp/8.0_dp)
    ELSE IF ( ( a < 0.0_dp ).AND.( b < 0.0_dp ) ) THEN
       MINIMUM = -MAXIMUM(-a,-b)
    ELSE IF ( ( a <= 0.0_dp ).AND.( b >= 0.0_dp ) ) THEN
       MINIMUM = a
    ELSE
       MINIMUM = b
    ENDIF
  END FUNCTION MINIMUM


  REAL(KIND=DP) FUNCTION MAXIMUM(a,b)
    !---------------------------------------------------------------------------
    ! called by :
    ! purpose :
    !     gives an approximation for the min(a,b) so that the resulting function
    !     is infinitely derivable and continuous with respect to a, b
    ! subroutine/function needed :
    ! input variables :
    !     a -> real number
    !     b -> real number
    ! output variables :
    ! results :
    !---------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(in) :: a
    REAL(KIND=DP), INTENT(in) :: b

    IF      ( ( a > 0.0_dp ).AND.( b > 0.0_dp ) ) THEN
       MAXIMUM = ( a**8.0_dp + b**8_dp )**(1.0_dp/8.0_dp)
    ELSE IF ( ( a < 0.0_dp ).AND.( b < 0.0_dp ) ) THEN
       MAXIMUM = -MINIMUM(-a,-b)
    ELSE IF ( ( a <= 0.0_dp ).AND.( b >= 0.0_dp ) ) THEN
       MAXIMUM = b
    ELSE
       MAXIMUM = a
    ENDIF
  END FUNCTION MAXIMUM


  subroutine JACO
  end subroutine JACO

END MODULE MODULE_TOOLS
