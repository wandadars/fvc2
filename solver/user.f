c ---------------------------------------------------------------------
      subroutine user_ic
c
      include 'MYDATA'

      if (ic_type .eq. 2) then
         call user_ic_fwd_step
      elseif (ic_type .eq. 3) then
         call user_ic_shocktube
      else
         call user_ic_uniform
      endif

      return
      end
c ---------------------------------------------------------------------
      subroutine user_ic_uniform
c
      include 'MYDATA'

      real*8 us(4),rpres,rdens,ru,rv,ret

      real*8 MixtPerf_et
      external MixtPerf_et

      rdens = ic_uniform(1)
      rpres = ic_uniform(2)
      ru    = ic_uniform(3)
      rv    = ic_uniform(4)

      us(1) = rdens
      us(2) = rdens*ru
      us(3) = rdens*rv
      us(4) = -1. ! dummy, set after
      ret   = MixtPerf_et(us(1),us(2),us(3),us(4),rpres)

      do i=1,nel
         u(i,1) = rdens     ! rho
         u(i,2) = rdens*ru  ! rho*u
         u(i,3) = rdens*rv  ! rho*v
         u(i,4) = rdens*ret ! rho*e_t
      enddo

      return
      end
c ---------------------------------------------------------------------
      subroutine user_ic_fwd_step
c
      include 'MYDATA'

      real*8 us(4),rpres,rdens,ru,rv,ret

      real*8 MixtPerf_et
      external MixtPerf_et

c     legacy forward-step initial condition
      rdens = 1.4
      rpres = 1.
      ru = 3.
      rv = 0.

      us(1) = rdens
      us(2) = rdens*ru
      us(3) = rdens*rv
      us(4) = -1. ! dummy, set after
      ret   = MixtPerf_et(us(1),us(2),us(3),us(4),rpres)

      do i=1,nel
         u(i,1) = rdens     ! rho
         u(i,2) = rdens*ru  ! rho*u
         u(i,3) = rdens*rv  ! rho*v
         u(i,4) = rdens*ret ! rho*e_t
      enddo

      return
      end
c ---------------------------------------------------------------------
      subroutine user_ic_shocktube
c
      include 'MYDATA'

      real*8 us(4),rpres,rdens,ru,rv,ret
      real*8 rx

      real*8 MixtPerf_et
      external MixtPerf_et

      do i=1,nel
         rx = elnorm(1,i) ! x centroid loc

         if (rx .le. ic_split_x) then
            rdens = ic_left(1)
            rpres = ic_left(2)
            ru    = ic_left(3)
            rv    = ic_left(4)
         else
            rdens = ic_right(1)
            rpres = ic_right(2)
            ru    = ic_right(3)
            rv    = ic_right(4)
         endif

         us(1) = rdens
         us(2) = rdens*ru
         us(3) = rdens*rv
         us(4) = -1. ! dummy, set after
         ret   = MixtPerf_et(us(1),us(2),us(3),us(4),rpres)

         u(i,1) = rdens     ! rho
         u(i,2) = rdens*ru  ! rho*u
         u(i,3) = rdens*rv  ! rho*v
         u(i,4) = rdens*ret ! rho*e_t
      enddo

      return
      end
c ---------------------------------------------------------------------
      subroutine user_source(rx,ry,rsource)
c
      include 'MYDATA'

      real*8 rsource(2)

      rsource(1) = 0. ! x vel for now source
      rsource(2) = 0. ! y vel for now source

      return
      end
c ---------------------------------------------------------------------
c    Evaluate gas quantities
c ---------------------------------------------------------------------
      FUNCTION MixtPerf_Ht(u1,u2,u3,u4,p)
      IMPLICIT NONE
      real*8       rgam,rspec,rcp,rcv,k
      common /rgasprops/ rgam,rspec,rcp,rcv,k
      REAL*8  p,u1,u2,u3,u4
      REAL*8  MixtPerf_Ht
      MixtPerf_Ht = (u4 + p)/u1
      return
      END
c ---------------------------------------------------------------------
      FUNCTION MixtPerf_T(u1,u2,u3,u4,p)
      IMPLICIT NONE
      real*8       rgam,rspec,rcp,rcv,k
      common /rgasprops/ rgam,rspec,rcp,rcv,k
      REAL*8  p,u1,u2,u3,u4
      REAL*8  MixtPerf_T
      MixtPerf_T = p/rspec/u1
      return
      END
c ---------------------------------------------------------------------
      FUNCTION MixtPerf_P(u1,u2,u3,u4)
      IMPLICIT NONE
      real*8       rgam,rspec,rcp,rcv,k
      common /rgasprops/ rgam,rspec,rcp,rcv,k
      REAL*8  u1,u2,u3,u4
      REAL*8  MixtPerf_P
      MixtPerf_P = (rgam-1.)*(u4-0.5*(u2*u2+u3*u3)/u1)
      return
      END
c ---------------------------------------------------------------------
      FUNCTION MixtPerf_et(u1,u2,u3,u4,p)
      IMPLICIT NONE
      real*8       rgam,rspec,rcp,rcv,k
      common /rgasprops/ rgam,rspec,rcp,rcv,k
      REAL*8  u1,u2,u3,u4,p
      REAL*8 MixtPerf_et
      MixtPerf_et = (p/(rgam-1.) + 0.5*(u2*u2 + u3*u3)/u1)/u1
      return
      END
c ---------------------------------------------------------------------
      FUNCTION MixtPerf_A(u1,u2,u3,u4,t)
      IMPLICIT NONE
      real*8       rgam,rspec,rcp,rcv,k
      common /rgasprops/ rgam,rspec,rcp,rcv,k
      REAL*8  t,u1,u2,u3,u4
      REAL*8  MixtPerf_A
      MixtPerf_A = sqrt(rgam*rspec*t)
      return
      end
c ---------------------------------------------------------------------
