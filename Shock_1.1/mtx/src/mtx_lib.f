c ***********************************************************************
!
!   Copyright (C) 2006, 2007  Bill Paxton
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
c ***********************************************************************


		module mtx_lib
		
		! mtx_lib has some utility routines for working with banded jacobians.
		
		! most of the mtx_lib routines come directly from LAPACK
		! and do not have definitions here
		
		implicit none
		
		
		contains
		
		subroutine write_jacobian_eqn_var_info(
     >			jac, sz, ldj, n, m, ml, mu, diag,
     >			ml_cells, mu_cells, var, eqn, jac_filename,
     >			form, out_of_bounds, io_unit, ierr)
			use alert_lib
     		integer, intent(in) :: sz, ldj, n, m, ml, mu, diag, 
     >			ml_cells, mu_cells, var, eqn, form, io_unit
			! diag is the row number for the diagonal elements of the jacobian
			! i.e., jacobian matrix diagonal (i,i) is stored in array element jac(diag,i)
			! diag is usually mu+ml+1 for numerical jacobian, mu+1 for analytic
			! ml_cells == number of cells in stencil with lower index
			! mu_cells == number of cells in stencil with higher index
			! total number of cells in stencil is thus = ml_cells + 1 + mu_cells
			! var == the variable number (1 to m)
			! eqn == the equation number (1 to m)
			! form == 0 means "raw" data
			! form == 1 means log10(max(1d-99,data)))
			! form == 2 means sign(1d0,data)*log10(1+data**2)
     		double precision, intent(in) :: jac(ldj,sz), out_of_bounds
			character (len=256), intent(in) :: jac_filename
			integer, intent(out) :: ierr
			
			
			integer :: k, kk, j, i
			double precision :: v
			character (len=256) :: message
			
			ierr = 0
			
			if (sz /= m*n) then
				message = 'sz /= m*n in write_banded_jacobian'
				call alert(ierr,message)
				return
			end if
			
			open(unit = io_unit, file = trim(jac_filename), iostat=ierr)
			if (ierr /= 0) then
				write(message,*) ierr, ' write_jacobian_eqn_var_info failed to open ', trim(jac_filename)
				call alert(ierr,message)
				return
			end if
			do kk = -ml_cells, mu_cells ! cell offset in stencil
				do k = 1, n ! cell
					j = (k-1)*m + var
					i = (k+kk-1)*m + eqn
					if (i >= j-mu .and. i <= j+ml) then
						v = jac(diag+i-j,j)
						if (form == 1) then
							v = log10(max(1d-99,v))
						else if (form == 2) then
							v = sign(1d0,v)*log10(1+v**2)
						else if (form /= 0) then
							message = 'invalid value for form in write_jacobian_eqn_var_info'
							ierr = -1
							call alert(ierr,message)
							return
						end if
					else
						v = out_of_bounds
					end if
					write(io_unit, '(e18.10)') v
				end do
			end do
			close(io_unit)
		
		end subroutine write_jacobian_eqn_var_info
		
		
		subroutine write_jacobian_bands(
     >			jac, sz, ldj, n, m, ml, mu, diag, names, 
     >			ml_cells, mu_cells,
     >			jac_filename, rows_filename, cols_filename, 
     >			form, out_of_bounds, ierr)
			use alert_lib
     		integer, intent(in) :: sz, ldj, n, m, ml, mu, diag, 
     >			ml_cells, mu_cells, form
			! diag is the row number for the diagonal elements of the jacobian
			! i.e., jacobian matrix diagonal (i,i) is stored in array element jac(diag,i)
			! diag is usually mu+ml+1 for numerical jacobian, mu+1 for analytic
			! ml_cells == number of cells in stencil with lower index
			! mu_cells == number of cells in stencil with higher index
			! total number of cells in stencil is thus = ml_cells + 1 + mu_cells
			! form == 0 means "raw" data
			! form == 1 means log10(max(1d-99,data)))
			! form == 2 means sign(1d0,data)*log10(1+data**2)
     		double precision, intent(in) :: jac(ldj,sz), out_of_bounds
			character (len=256), intent(in) :: names(m), jac_filename, rows_filename, cols_filename
			integer, intent(out) :: ierr
			
			integer :: var, eqn, k, kk, j, i, qcnt, io_unit
			double precision :: v
			character (len=256) :: message
			
			io_unit = 56
			ierr = 0
			
			if (sz /= m*n) then
				message = 'sz /= m*n in write_banded_jacobian'
				call alert(ierr,message)
				return
			end if
			
			do var = 1, m
				qcnt = 0
				open(unit = io_unit, file = trim(jac_filename), iostat=ierr)
				if (ierr /= 0) then
					write(message,*) 'write_jacobian_bands failed to open ', trim(jac_filename)
					call alert(ierr,message)
					return
				end if
				
				do eqn = 1, m
					qcnt = qcnt + ml_cells + mu_cells + 1
					call write_jacobian_eqn_var_info(
     >				jac, sz, ldj, n, m, ml, mu, diag,
     >				ml_cells, mu_cells, var, eqn, jac_filename,
     >				form, out_of_bounds, io_unit, ierr)
				end do
				
				close(io_unit)
				
			end do
			
			open(unit = io_unit, file = trim(rows_filename), iostat=ierr)
			if (ierr /= 0) then
				write(message,*) 'write_jacobian_bands failed to open ', trim(rows_filename)
				call alert(ierr,message)
				return
			end if
			do i = 1, qcnt
				write(io_unit,'(e18.10)') dble(i) - 0.5
			end do
			close(io_unit)

			open(unit = io_unit, file = trim(cols_filename), iostat=ierr)
			if (ierr /= 0) then
				write(message,*) 'write_jacobian_bands failed to open ', trim(cols_filename)
				call alert(ierr,message)
				return
			end if
			do i = 1, n
				write(io_unit,'(i6)') i
			end do
			close(io_unit)
		
		end subroutine write_jacobian_bands


		subroutine show_banded_jacobian(jac,sz,ldj,mu,ml,m,n)
     		integer, intent(in) :: sz, ldj, mu, ml, m, n
     		double precision, intent(inout) :: jac(ldj,sz)
			
			integer :: i, j, im, in
			double precision :: v
			
			if (m*n /= sz) then
				write(*,*) 'incorrect args for show_banded_jacobian: expect m*n = sz.'
				write(*,*) '		m = number of components per point'
				write(*,*) '		n = number of points'
				write(*,*) '	  sz = number of equations'
				return
			end if
			
			im = 1; in = 1
			write(*,fmt='(12x)',advance='no')
			do j = 1, sz
				write(*,fmt='(3x,"(",i2,",",i2,")",4x)',advance='no') im, in
				im = im + 1
				if (im > m) then
					im = 1; in = in + 1
				end if
			end do
			write(*,*)
			im = 1; in = 1
			do i = 1, sz
				write(*,fmt='("(",i2,",",i2,")",4x)',advance='no') im, in
				im = im + 1
				if (im > m) then
					im = 1; in = in + 1
				end if
				do j = 1, sz
					if (i >= max(1,j-mu) .and. i <= min(sz,j+ml)) then
						v = jac(ml+mu+1+i-j,j)
						if (v == 0d0) then
							write(*,fmt='(6x,"0",7x)',advance='no')
						else
							write(*,fmt='(e14.4)',advance='no') v
						end if
					else
						write(*,fmt='(14x)',advance='no')
					end if
				end do
				write(*,*)
			end do
			write(*,*)
			
		end subroutine show_banded_jacobian


      subroutine write_eqn_var_partials(
     >         n, m, io_unit, var, eqn, em2, em1, e00, ep1, ep2, filename, ierr)
         use alert_lib
         integer, intent(in) :: n, m, io_unit, var, eqn
         double precision, dimension(m,m,n), intent(in) :: em2, em1, e00, ep1, ep2
         character (len=256), intent(in) :: filename
         integer, intent(out) :: ierr
         
         character (len=256) :: message

         open(unit = io_unit, file = trim(filename), iostat=ierr)
         if (ierr /= 0) then
            write(message,*) ierr, ' write_eqn_var_partials failed to open ', trim(filename)
            call alert(ierr,message)
            return
         end if
         
         write(io_unit, '(e18.10)') 0d0, 0d0, em2(eqn,var,3:n) 
         write(io_unit, '(e18.10)') 0d0, em1(eqn,var,2:n)   
         write(io_unit, '(e18.10)') e00(eqn,var,1:n)  
         write(io_unit, '(e18.10)') ep1(eqn,var,1:n-1), 0d0
         write(io_unit, '(e18.10)') ep2(eqn,var,1:n-2), 0d0, 0d0

         close(io_unit)
         
         contains
         
         subroutine write_blanks
            integer :: i
            do i = 1, n
               write(io_unit, '(e18.10)') -100d0
            end do
         end subroutine write_blanks

      end subroutine write_eqn_var_partials

      
      subroutine copy_to_jacobian(jac, sz, ldj, n, m, ml, mu, diag, em2, em1, e00, ep1, ep2, use_OpenMP)
         integer, intent(in) :: sz, ldj, n, m, ml, mu, diag
         double precision, dimension(m,m,n), intent(in) :: em2, em1, e00, ep1, ep2
         double precision, intent(out) :: jac(ldj,sz)
         logical, intent(in) :: use_OpenMP
         
         integer :: eqn

         jac = 0
         if (use_OpenMP) then
!$OMP PARALLEL DO PRIVATE(eqn)
            do eqn = 1, m
               call copy_eqn_to_jacobian(jac, sz, ldj, n, m, ml, mu, diag, eqn, em2, em1, e00, ep1, ep2)
            end do
!$OMP END PARALLEL DO
         else
            do eqn = 1, m
               call copy_eqn_to_jacobian(jac, sz, ldj, n, m, ml, mu, diag, eqn, em2, em1, e00, ep1, ep2)
            end do
         end if

      end subroutine copy_to_jacobian
      
      
      subroutine copy_eqn_to_jacobian(jac, sz, ldj, n, m, ml, mu, diag, eqn, em2, em1, e00, ep1, ep2)
         integer, intent(in) :: sz, ldj, n, m, ml, mu, diag
         ! diag is the row number for the diagonal elements of the jacobian
         ! i.e., jacobian matrix diagonal (i,i) is stored in array element jac(diag,i)
         ! diag is usually mu+ml+1 for numerical jacobian, mu+1 for analytic
         double precision, dimension(m,m,n), intent(in) :: em2, em1, e00, ep1, ep2
         double precision, intent(out) :: jac(ldj,sz)
         integer, intent(in) :: eqn
         
         integer :: var

         do var = 1, m
            call copy_QV_to_jacobian(jac, sz, ldj, n, m, ml, mu, diag, eqn, var, em2, em1, e00, ep1, ep2)
         end do

      end subroutine copy_eqn_to_jacobian
      
      
      subroutine copy_QV_to_jacobian(jac, sz, ldj, n, m, ml, mu, diag, eqn, var, em2, em1, e00, ep1, ep2)
         integer, intent(in) :: sz, ldj, n, m, ml, mu, diag
         double precision, dimension(m,m,n), intent(in) :: em2, em1, e00, ep1, ep2
         double precision, intent(out) :: jac(ldj,sz)
         integer, intent(in) :: eqn, var
         
         integer :: k, dk, ii, jj

         do dk = -2, 2
            ii = eqn - var - m*dk + diag
            jj = var + m*(dk-1)
            select case(dk)
               case(-2) 
                  forall (k=3:n)   jac(ii,jj+m*k) = em2(eqn,var,k)
               case(-1) 
                  forall (k=2:n)   jac(ii,jj+m*k) = em1(eqn,var,k)
               case(0) 
                  forall (k=1:n)   jac(ii,jj+m*k) = e00(eqn,var,k)
               case(1) 
                  forall (k=1:n-1) jac(ii,jj+m*k) = ep1(eqn,var,k)
               case(2) 
                  forall (k=1:n-2) jac(ii,jj+m*k) = ep2(eqn,var,k)
            end select
         end do

      end subroutine copy_QV_to_jacobian

      
      subroutine get_from_jacobian(
     >         jac, sz, ldj, n, m, ml, mu, diag, em2, em1, e00, ep1, ep2, use_OpenMP)
         integer, intent(in) :: sz, ldj, n, m, ml, mu, diag
         ! diag is the row number for the diagonal elements of the jacobian
         ! i.e., jacobian matrix diagonal (i,i) is stored in array element jac(diag,i)
         ! diag is usually mu+ml+1 for numerical jacobian, mu+1 for analytic
         ! code assumes ml_cells = 2 and mu_cells = 2
         double precision, intent(in) :: jac(ldj,sz)
         double precision, dimension(m,m,n), intent(out) :: em2, em1, e00, ep1, ep2
         logical, intent(in) :: use_OpenMP
         
         integer :: eqn

         if (use_OpenMP) then
!$OMP PARALLEL DO PRIVATE(eqn)
            do eqn = 1, m
               call do_eqn(eqn)
            end do
!$OMP END PARALLEL DO
         else
            do eqn = 1, m
               call do_eqn(eqn)
            end do
         end if

      contains
      
      subroutine do_eqn(eqn)
         integer, intent(in) :: eqn
         
         integer :: k, dk, var, ii, jj

         do dk = -2, 2
            do var = 1, m
               ii = eqn - var - m*dk + diag
               jj = var + m*(dk-1)
               select case(dk)
                  case(-2) 
                     forall (k=3:n)   em2(eqn,var,k) = jac(ii,jj+m*k)
                     em2(eqn,var,1:2) = 0
                  case(-1) 
                     forall (k=2:n)   em1(eqn,var,k) = jac(ii,jj+m*k)
                     em1(eqn,var,1) = 0
                  case(0) 
                     forall (k=1:n)   e00(eqn,var,k) = jac(ii,jj+m*k)
                  case(1) 
                     forall (k=1:n-1) ep1(eqn,var,k) = jac(ii,jj+m*k)
                     ep1(eqn,var,n) = 0
                  case(2) 
                     forall (k=1:n-2) ep2(eqn,var,k) = jac(ii,jj+m*k)
                     ep2(eqn,var,n-1:n) = 0
               end select
            end do
         end do

      end subroutine do_eqn
      
         
      end subroutine get_from_jacobian
		

		end module mtx_lib
