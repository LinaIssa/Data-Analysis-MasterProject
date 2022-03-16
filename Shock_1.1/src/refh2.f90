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

MODULE refH2
  IMPLICIT none

  INCLUDE "precision.f90"

  PUBLIC coolH2

  ! variables from "precision.f90" are private
  PRIVATE :: DP, LONG, minus_infinity, plus_infinity

  INTEGER,       PARAMETER                  :: npt = 23
  INTEGER,       PARAMETER                  :: npd = 49 
  INTEGER,       PARAMETER                  :: nph = 36
  INTEGER,       PARAMETER                  :: npr =  6

  REAL(KIND=dp), DIMENSION(npt)             :: T_N
  REAL(KIND=dp), DIMENSION(npd)             :: DN
  REAL(KIND=dp), DIMENSION(nph)             :: RH
  REAL(KIND=dp), DIMENSION(npr)             :: ROP

  REAL(KIND=dp), DIMENSION(npt,npd,nph,npr) :: WG

  REAL(KIND=dp), PARAMETER                  :: T_Nmin = 100
  REAL(KIND=dp), PARAMETER                  :: T_Nmax = 1.e4
  REAL(KIND=dp), PARAMETER                  :: DNmin = 1.
  REAL(KIND=dp), PARAMETER                  :: DNmax = 1.0d+08

CONTAINS

!!$  SUBROUTINE setrefH2
!!$    logical,save::first=.true.
!!$    real WH2,tt,dd,rr,hh
!!$    integer i
!!$
!!$    ! Au premier appel, charge les données
!!$    if (first) then
!!$       first=.false.
!!$       call initWH2
!!$    endif
!!$
!!$    ! Interpole dans le cube 
!!$    do i=1,nz
!!$       tt=log10(Tn(i))
!!$       dd=log10(nH(i)*unitNchim)
!!$       hh=log10(dens(i,H)/n_H2(i))
!!$       rr=3.
!!$       CALL coolH2(tt,dd,hh,rr,WH2)
!!$       H2_cooling(i)=-n_H2(i)*10.**WH2*(unitNchim*erg/cm3/seconde)
!!$    enddo
!!$
!!$  end SUBROUTINE setrefH2

  SUBROUTINE coolH2(tt, dd, hh, rr, WH2)

    ! - Compute interpolated value in 4D-grid
    REAL(KIND=dp) :: WH2
    REAL(KIND=dp) :: tt,dd,hh,rr
    REAL(KIND=dp) :: t,u,v,w
    INTEGER       ::  i,it,id,ir,ih

    REAL(KIND=dp) :: y0000
    REAL(KIND=dp) :: y0001
    REAL(KIND=dp) :: y0010
    REAL(KIND=dp) :: y0011
    REAL(KIND=dp) :: y0100
    REAL(KIND=dp) :: y0101
    REAL(KIND=dp) :: y0110
    REAL(KIND=dp) :: y0111
    REAL(KIND=dp) :: y1000
    REAL(KIND=dp) :: y1001
    REAL(KIND=dp) :: y1010
    REAL(KIND=dp) :: y1011
    REAL(KIND=dp) :: y1100
    REAL(KIND=dp) :: y1101
    REAL(KIND=dp) :: y1110
    REAL(KIND=dp) :: y1111

    logical,save::first=.true.

    ! Au premier appel, charge les données
    if (first) then
       first=.false.
       call initWH2
    endif


!!$    ! - Check that the point is inside our grid
!!$
!!$    print *, tt, dd, hh, rr
!!$
!!$    if (tt .lt. T_N(1) .or. tt .gt. T_N(npt)) then
!!$       print *, " Value out of range (tt)!"
!!$       print *, " T should be in the range:", 10.**T_N(1), " - ", 10.**T_N(npt)
!!$    endif
!!$
!!$    if (dd .lt. DN(1) .or. dd .gt. DN(npd)) then
!!$       print *, " Value out of range (dd)!"
!!$       print *, " nH should be in the range:", 10.**DN(1), " - ", 10.**DN(npd)
!!$    endif
!!$
!!$    if (hh .lt. RH(1) .or. hh .gt. RH(nph)) then
!!$       print *, " Value out of range (hh)"
!!$       print *, " H/H2 should be in the range:", 10.**RH(1), " - ", 10.**RH(nph)
!!$    endif
!!$
!!$    if (rr .lt. ROP(1) .or. rr .gt. ROP(npr)) then
!!$       print *, " Value out of range (rr)"
!!$       print *, " O/P should be in the range:", ROP(1), " - ", ROP(npr)
!!$    endif

    !  Find place in grid

    it=1
    do i=1,npt-1
       if (T_N(i)<tt) it=i
    enddo

    id=1
    do i=1,npd-1
       if (DN(i)<dd) id=i
    enddo

    ih=1
    do i=1,nph-1
       if (RH(i)<hh) ih=i
    enddo

    ir=1
    do i=1,npr-1
       if (ROP(i)<rr) ir=i
    enddo


    !  Do a simple linear interpolation

    y0000 = WG(it  ,id  ,ih  ,ir  )
    y1000 = WG(it+1,id  ,ih  ,ir  )
    y0100 = WG(it  ,id+1,ih  ,ir  )
    y1100 = WG(it+1,id+1,ih  ,ir  )
    y0010 = WG(it  ,id  ,ih+1,ir  )
    y1010 = WG(it+1,id  ,ih+1,ir  )
    y0110 = WG(it  ,id+1,ih+1,ir  )
    y1110 = WG(it+1,id+1,ih+1,ir  )
    y0001 = WG(it  ,id  ,ih  ,ir+1)
    y1001 = WG(it+1,id  ,ih  ,ir+1)
    y0101 = WG(it  ,id+1,ih  ,ir+1)
    y1101 = WG(it+1,id+1,ih  ,ir+1)
    y0011 = WG(it  ,id  ,ih+1,ir+1)
    y1011 = WG(it+1,id  ,ih+1,ir+1)
    y0111 = WG(it  ,id+1,ih+1,ir+1)
    y1111 = WG(it+1,id+1,ih+1,ir+1)

    t = (tt - T_N(it)) / (T_N(it+1) - T_N(it))
    u = (dd - DN(id)) / (DN(id+1) - DN(id))
    v = (hh - RH(ih)) / (RH(ih+1) - RH(ih))
    w = (rr - ROP(ir)) / (ROP(ir+1) - ROP(ir))

    WH2 =  (1.-t) * (1.-u) * (1.-v) * (1.-w) * y0000&
         +    t   * (1.-u) * (1.-v) * (1.-w) * y1000&
         + (1.-t) *    u   * (1.-v) * (1.-w) * y0100&
         +    t   *    u   * (1.-v) * (1.-w) * y1100&
         + (1.-t) * (1.-u) *    v   * (1.-w) * y0010&
         +    t   * (1.-u) *    v   * (1.-w) * y1010&
         + (1.-t) *    u   *    v   * (1.-w) * y0110&
         +    t   *    u   *    v   * (1.-w) * y1110&
         + (1.-t) * (1.-u) * (1.-v) *    w   * y0001&
         +    t   * (1.-u) * (1.-v) *    w   * y1001&
         + (1.-t) *    u   * (1.-v) *    w   * y0101&
         +    t   *    u   * (1.-v) *    w   * y1101&
         + (1.-t) * (1.-u) *    v   *    w   * y0011&
         +    t   * (1.-u) *    v   *    w   * y1011&
         + (1.-t) *    u   *    v   *    w   * y0111&
         +    t   *    u   *    v   *    w   * y1111

    return
  end subroutine coolH2

  subroutine initWH2

    implicit none
    integer it,id,ih,io
    !  Read Cooling file (binary)

    open (10, file="input/le_cube", status="old")

    read(10,*) RH
    read(10,*) ROP
    read(10,*) T_N
    read(10,*) DN
    read(10,*) ((((WG(it,id,ih,io), id=1,npd), it=1,npt), io=1,npr), ih=1,nph)

    close (10)

    return
  end subroutine initWH2

end MODULE refH2


