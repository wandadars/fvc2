c ---------------------------------------------------------------------
      program project
c
      include 'MYDATA'

      call rzero(elnorm(1,1),11*lm)
      call rzero(nodel(1,1),2*lm)
      call izero(elml(1,1),9*lm)
      call izero(edgl(1,1),6*lm)

c     -------------------- setup routines -----------------------------
      call read_user_input_file
      call read_sorted_grid
      call user_ic
c     -------------------- end setup routines -------------------------

c     -------------------- solver routines ----------------------------
      time = 0.
      do istep=0,nstepmax
         if (modulo(istep,iostep).eq.0) then
            call write_soln ! output soln
         endif
c         !King, modified for last time only
c!         if (istep .eq. nstepmax) then
c!            call write_soln
c!         endif
c!         if (modulo(istep,iostep) .eq. 0) then
c!            print*, istep, nstepmax
c!         endif

c        rcmax = 0.200
c        rumax = 300.
c        rvmax = 300.
c        rdxmin = 100.
c        do i =1,nel
c           if (abs(u(i,2)/u(i,1)) .gt. rumax) rumax = u(i,2)/u(i,1)
c           if (abs(u(i,3)/u(i,1)) .gt. rvmax) rvmax = u(i,3)/u(i,1)
c           if (sqrt(elnorm(3,i)) .lt. rdxmin) rdxmin= sqrt(elnorm(3,i))
c        enddo
c        dt = rcmax/(rumax/rdxmin + rvmax/rdxmin)
         time = time + dt

         if (time_scheme .eq. 1) then ! forward euler

            call zero_matricies(b,istep)

            do ied=1,ned ! loop over edges, main loop
               call compute_edge_values
               
               if (edgl(5,ied).eq.0) then ! interior edges
                  call convective_flux_interior(b)
               else
                  call convective_flux_boundary(b)
               endif
      
c              call compute_sources(A,b)
            enddo

            call solve_equations_fwdE(b,1)

         elseif (time_scheme .eq. 2) then ! rk4

            ! RK4 - code here

         endif




      enddo
c     -------------------- end solver routines ------------------------

      stop
      end
c ---------------------------------------------------------------------
      subroutine rzero(rdum,n)

      real*8 rdum(n)

      do i=1,n
         rdum(i) = 0.
      enddo

      return
      end
c ---------------------------------------------------------------------
      subroutine izero(idum,n)

      integer idum(n)

      do i=1,n
         idum(i) = 0
      enddo

      return
      end
c ---------------------------------------------------------------------
