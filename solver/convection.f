c ---------------------------------------------------------------------
      subroutine convective_flux_interior(bdum)
c
      include 'MYDATA'

c     real*8 bdum(nmm),Adum(nmm,nmm),flux(4),ul1(4),ur1(4)
      real*8 bdum(nmm),flux(4),ul1(4),ur1(4),nhatx,nhaty

      ul1(1) = u(iel1,1)
      ul1(2) = u(iel1,2)
      ul1(3) = u(iel1,3)
      ul1(4) = u(iel1,4)

      ur1(1) = u(iel2,1)
      ur1(2) = u(iel2,2)
      ur1(3) = u(iel2,3)
      ur1(4) = u(iel2,4)

      nhatx = elnorm(4+is1*2,iel1)
      nhaty = elnorm(5+is1*2,iel1)

      call convective_flux(ul1,ur1,nhatx,nhaty,flux)

      bdum(iems1r)  = bdum(iems1r)  - flux(1)*deta
      bdum(iems1ru) = bdum(iems1ru) - flux(2)*deta
      bdum(iems1rv) = bdum(iems1rv) - flux(3)*deta
      bdum(iems1re) = bdum(iems1re) - flux(4)*deta

      bdum(iems2r)  = bdum(iems2r)  + flux(1)*deta
      bdum(iems2ru) = bdum(iems2ru) + flux(2)*deta
      bdum(iems2rv) = bdum(iems2rv) + flux(3)*deta
      bdum(iems2re) = bdum(iems2re) + flux(4)*deta

      return
      end
c ---------------------------------------------------------------------
      subroutine convective_flux_boundary(bdum)
c
      include 'MYDATA'

c     real*8 bdum(nmm),Adum(nmm,nmm),flux(4),ul1(4), ur1(4)
      real*8 bdum(nmm),flux(4),ul1(4),ur1(4),rnhatx,rnhaty

      idum = edgl(5,ied) ! location in bcvals(1:4,idum) with bc data

      ! supersonic outflow
      if (idum .eq. 3) then
c        rnhatx = elnorm(4+is1*2,iel1)
c        rnhaty = elnorm(5+is1*2,iel1)

c        rdens = u(iel1,1)
c        ru    = u(iel1,2)
c        rv    = u(iel1,3)
c        ret   = u(iel1,4)

c        rpres = MixtPerf_P(rdens,ru,rv,ret)

c        rdens = u(iel1,1)
c        ru    = u(iel1,2)/rdens
c        rv    = u(iel1,3)/rdens
c        ret   = u(iel1,4)/rdens

c        write(6,*) rdens,ru,rv,ret,elnorm(1,iel1),elnorm(2,iel1)

c        flux(1) = rdens*(ru*rnhatx + rv*rnhaty)
c        flux(2) = rdens*ru*(ru*rnhatx + rv*rnhaty) + rpres*rnhatx
c        flux(3) = rdens*rv*(ru*rnhatx + rv*rnhaty) + rpres*rnhaty
c        flux(4) = (ret + rpres/rdens)*rdens*(ru*rnhatx + rv*rnhaty)

c        bdum(iems1r)  = bdum(iems1r)  - flux(1)*deta
c        bdum(iems1ru) = bdum(iems1ru) - flux(2)*deta
c        bdum(iems1rv) = bdum(iems1rv) - flux(3)*deta
c        bdum(iems1re) = bdum(iems1re) - flux(4)*deta

         ! left state (i.e. this cell center)
         ul1(1) = u(iel1,1)
         ul1(2) = u(iel1,2)
         ul1(3) = u(iel1,3)
         ul1(4) = u(iel1,4)

         rnhatx = elnorm(4+is1*2,iel1)
         rnhaty = elnorm(5+is1*2,iel1)

         ur1(1) = u(iel1,1)
         ur1(2) = u(iel1,2)
         ur1(3) = u(iel1,3)
         ur1(4) = u(iel1,4)


         call convective_flux(ul1,ur1,rnhatx,rnhaty,flux)
c        rpres = (rgam-1.)*(u(iel1,4)- 0.)
c        flux(1) = 0.
c        flux(2) = rpres*rnhatx
c        flux(3) = rpres*rnhaty
c        flux(4) = 0.

c        flux(1) = 0.
c        flux(2) = 0.
c        flux(3) = 0.
c        flux(4) = 0.

         bdum(iems1r)  = bdum(iems1r)  - flux(1)*deta
         bdum(iems1ru) = bdum(iems1ru) - flux(2)*deta
         bdum(iems1rv) = bdum(iems1rv) - flux(3)*deta
         bdum(iems1re) = bdum(iems1re) - flux(4)*deta

      ! supersonic inflow
      elseif (idum .eq. 2) then
         rnhatx = elnorm(4+is1*2,iel1)
         rnhaty = elnorm(5+is1*2,iel1)

c        rdens0 = 1.2
c        rpres0 = 101000.
c        rtemp0 = 293.263

c        rmach_inlet = 3.
c        rdens_rdens0 = 0.07622631
c        rpres_rpres0 = 0.02722368
c        rtemp_rtemp0 = 0.35714285

c        rdens_rdens0 = 1.
c        rpres_rpres0 = 1.
c        rtemp_rtemp0 = 1.

         
c        rtemp = rtemp_rtemp0*rtemp0
c        rdens = rdens_rdens0*rdens0
c        rpres = rpres_rpres0*rpres0
c        rsoun = sqrt(rgam*rspec*rtemp)


         rdens = bcvals(1,idum)
         rpres = bcvals(2,idum)
         ru    = bcvals(3,idum)
         rv    = bcvals(4,idum)
         if (rdens .le. 0. .or. rpres .le. 0.) then
            rdens = ic_uniform(1)
            rpres = ic_uniform(2)
            ru    = ic_uniform(3)
            rv    = ic_uniform(4)
         endif
c        ru = rsoun*rmach_inlet
c        rv = 0.
         ret = (rpres/(rgam-1.) + 0.5*rdens*(ru*ru + rv*rv))/rdens

         flux(1) = rdens*(ru*rnhatx + rv*rnhaty)
         flux(2) = rdens*ru*(ru*rnhatx + rv*rnhaty) + rpres*rnhatx
         flux(3) = rdens*rv*(ru*rnhatx + rv*rnhaty) + rpres*rnhaty
         flux(4) = (ret + rpres/rdens)*rdens*(ru*rnhatx + rv*rnhaty)

         bdum(iems1r)  = bdum(iems1r)  - flux(1)*deta
         bdum(iems1ru) = bdum(iems1ru) - flux(2)*deta
         bdum(iems1rv) = bdum(iems1rv) - flux(3)*deta
         bdum(iems1re) = bdum(iems1re) - flux(4)*deta

c     ! no penetration wall bc
      elseif (idum .eq. 1) then

         ! left state (i.e. this cell center)
         ul1(1) = u(iel1,1)
         ul1(2) = u(iel1,2)
         ul1(3) = u(iel1,3)
         ul1(4) = u(iel1,4)

         rnhatx = elnorm(4+is1*2,iel1)
         rnhaty = elnorm(5+is1*2,iel1)

         ! right state (i.e. wall mirrored)
         ur1(1) = u(iel1,1)    ! mirror
         if (abs(rnhaty) .lt. 1E-9) then
            ur1(2) = 0.
            ur1(3) = u(iel1,3)
         elseif (abs(rnhatx) .lt. 1E-9) then
            ur1(2) = u(iel1,2)    ! no flow
            ur1(3) = 0.
         else
            ur1(2) = u(iel1,2)*0.    ! no flow
            ur1(3) = -ur1(2)*rnhatx/rnhaty
         endif
c        ur1(3) = u(iel1,3)    ! no flow
         ur1(4) = u(iel1,4)    ! mirror

         call convective_flux(ul1,ur1,rnhatx,rnhaty,flux)
c        rpres = (rgam-1.)*(u(iel1,4)- 0.)
c        flux(1) = 0.
c        flux(2) = rpres*rnhatx
c        flux(3) = rpres*rnhaty
c        flux(4) = 0.

c        flux(1) = 0.
c        flux(2) = 0.
c        flux(3) = 0.
c        flux(4) = 0.

         bdum(iems1r)  = bdum(iems1r)  - flux(1)*deta
         bdum(iems1ru) = bdum(iems1ru) - flux(2)*deta
         bdum(iems1rv) = bdum(iems1rv) - flux(3)*deta
         bdum(iems1re) = bdum(iems1re) - flux(4)*deta

      endif

      return
      end
c ---------------------------------------------------------------------
      subroutine convective_flux(ul1,ur1,nhatx,nhaty,flux)
c
      include 'MYDATA' 
c
      real*8 ul1(4),ur1(4),nhatx,nhaty,flux(4)

      if (conv_scheme .eq. 0) then
         ! do nothing for convective term
         flux(1) = 0.
         flux(2) = 0.
         flux(3) = 0.
         flux(4) = 0.

      elseif (conv_scheme .eq. 1) then
         ! AUSM scheme
         call ausm_plus(ul1,ur1,nhatx,nhaty,flux)

      elseif (conv_scheme .eq. 2) then
         ! Roe scheme
         call Roe(ul1,ur1,nhatx,nhaty,flux)

      endif

      return
      end
c ---------------------------------------------------------------------
!*********************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Diff.
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*********************************************************************
      subroutine Roe(uL, uR, nx, ny,flux)
      real*8 uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
      real*8 nx, ny       !  Input: face normal vector, [nx, ny] 
                           !         (Left-to-Right)
      real*8 flux(4)      ! Output: Roe flux function (upwind)
!Local constants
      real*8 gamma                          ! Ratio of specific heat.
      real*8 zero, fifth, half, one, two    ! Numbers
!Local variables
      real*8 tx, ty       ! Tangent vector (perpendicular 
                           ! to the face normal)
      real*8 vxL, vxR, vyL, vyR             ! Velocity components.
      real*8 rhoL, rhoR, pL, pR             ! Primitive variables.
      real*8 vnL, vnR, vtL, vtR             ! Normal and tangent vel
      real*8 aL, aR, HL, HR, tL, tR         ! Speeds of sound and temps.
      real*8 RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
      real*8 drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
      real*8 ws(4),dws(4), Rv(4,4)          ! Wave speeds and 
                                             ! right-eigevectors
      real*8 fL(4), fR(4), diss(4)          ! Fluxes ad dissipation 

      real*8 MixtPerf_Ht,MixtPerf_P,MixtPerf_T,MixtPerf_A
      external MixtPerf_Ht,MixtPerf_P,MixtPerf_T,MixtPerf_A

      real*8       rgam,rspec,rcp,rcv,k
      common /rgasprops/ rgam,rspec,rcp,rcv,k

!Constants.
      gamma = rgam
      zero = 0.0
      fifth = 0.2
      half = 0.5
      one = 1.0
      two = 2.0

!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
      tx = -ny
      ty = nx

!Primitive and other variables.
!  Left state
      rhoL = uL(1)
      vxL = uL(2)/uL(1)
      vyL = uL(3)/uL(1)
      vnL = vxL*nx+vyL*ny
      vtL = vxL*tx+vyL*ty

      pL = MixtPerf_P(uL(1),uL(2),uL(3),uL(4))
      tL = MixtPerf_T(uL(1),uL(2),uL(3),uL(4),pL)
      aL = MixtPerf_A(uL(1),uL(2),uL(3),uL(4),tL)
      hL = MixtPerf_Ht(uL(1),uL(2),uL(3),uL(4),pL)
c     pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
c     aL = sqrt(gamma*pL/rhoL)
c     HL = ( uL(4) + pL ) / rhoL
!  Right state
      rhoR = uR(1)
      vxR = uR(2)/uR(1)
      vyR = uR(3)/uR(1)
      vnR = vxR*nx+vyR*ny
      vtR = vxR*tx+vyR*ty
c     pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
c     aR = sqrt(gamma*pR/rhoR)
c     HR = ( uR(4) + pR ) / rhoR
      pR = MixtPerf_P(uR(1),uR(2),uR(3),uR(4))
      tR = MixtPerf_T(uR(1),uR(2),uR(3),uR(4),pR)
      aR = MixtPerf_A(uR(1),uR(2),uR(3),uR(4),tR)
      hR = MixtPerf_Ht(uR(1),uR(2),uR(3),uR(4),pR)

!First compute the Roe Averages
      RT = sqrt(rhoR/rhoL)
      rho = RT*rhoL
      vx = (vxL+RT*vxR)/(one+RT)
      vy = (vyL+RT*vyR)/(one+RT)
      H = ( HL+RT* HR)/(one+RT)
      a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
      vn = vx*nx+vy*ny
      vt = vx*tx+vy*ty

!Wave Strengths
      drho = rhoR - rhoL 
      dp =   pR - pL
      dvn =  vnR - vnL
      dvt =  vtR - vtL

      dV(1) = (dp - rho*a*dvn )/(two*a*a)
      dV(2) = rho*dvt/a
      dV(3) =  drho - dp/(a*a)
      dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speed
      ws(1) = abs(vn-a)
      ws(2) = abs(vn)
      ws(3) = abs(vn)
      ws(4) = abs(vn+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
      dws(1) = fifth
      if ( ws(1) .lt. dws(1) ) ws(1) = half * 
     >                     ( ws(1)*ws(1)/dws(1)+dws(1) )
      dws(4) = fifth
      if ( ws(4) .gt. dws(4) ) ws(4) = half * 
     >                     ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
      Rv(1,1) = one    
      Rv(2,1) = vx - a*nx
      Rv(3,1) = vy - a*ny
      Rv(4,1) =  H - vn*a

      Rv(1,2) = zero
      Rv(2,2) = a*tx
      Rv(3,2) = a*ty
      Rv(4,2) = vt*a

      Rv(1,3) = one
      Rv(2,3) = vx
      Rv(3,3) = vy 
      Rv(4,3) = half*(vx*vx+vy*vy)

      Rv(1,4) = one
      Rv(2,4) = vx + a*nx
      Rv(3,4) = vy + a*ny
      Rv(4,4) =  H + vn*a

!Dissipation Term
      diss = zero
      do i=1,4
         do j=1,4
            diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
         enddo
      enddo

!Compute the flux.
      fL(1) = rhoL*vnL
      fL(2) = rhoL*vnL * vxL + pL*nx
      fL(3) = rhoL*vnL * vyL + pL*ny
      fL(4) = rhoL*vnL *  HL

      fR(1) = rhoR*vnR
      fR(2) = rhoR*vnR * vxR + pR*nx
      fR(3) = rhoR*vnR * vyR + pR*ny
      fR(4) = rhoR*vnR *  HR

      flux(1) = half * (fL(1) + fR(1) - diss(1))
      flux(2) = half * (fL(2) + fR(2) - diss(2))
      flux(3) = half * (fL(3) + fR(3) - diss(3))
      flux(4) = half * (fL(4) + fR(4) - diss(4))

      return
      end
! ******************************************************************************
!
! Purpose: Compute convective fluxes using AUSM+ scheme.
!
! Description: None.
!
! Input: 
!   nx          x-component of face normal
!   ny          y-component of face normal
!   nz          z-component of face normal
!   nm          Magnitude of face normal
!   fs          Face speed
!   rl          Density of left state
!   ul          x-component of velocity of left state
!   vl          y-component of velocity of left state
!   wl          z-component of velocity of left state   
!   Hl		Total enthalpy of left state
!   al		Speed of sound of left state
!   pl          Pressure of left state
!   rr          Density of right state
!   ur          x-component of velocity of right state
!   vr          y-component of velocity of right state
!   wr          z-component of velocity of right state  
!   pr          Pressure of right state
!   Hr		Total enthalpy of right state
!   ar		Speed of sound of right state
!
! Output: 
!   flx         Fluxes
!   vf          Face velocities ! NOT USED IN CMT-NEK YET
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSM+, AIAA Paper
!      95-1701, 1995
!   2. Do not use computation of face speed of sound which leads to exact 
!      capturing of isolated normal shock waves because of robustness problems
!      for unsteady flows and because that formulation is not applicable to 
!      anything but calorically and thermally perfect gases.
!
! ******************************************************************************

      SUBROUTINE ausm_plus(rul,rur,nx,ny,flx)

!     IMPLICIT NONE ! HAHAHHAHHAHA
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************
      real*8 MixtPerf_Ht,MixtPerf_P,MixtPerf_T,MixtPerf_A
      external MixtPerf_Ht,MixtPerf_P,MixtPerf_T,MixtPerf_A

! ==============================================================================
! Arguments
! ==============================================================================
      REAL*8 al,ar,fs,nx,ny,
     >     nz,pl,pr,rl,rr,ul,
     >     ur,vl,vr,wl,wr,cpl,
     >     cpr,tl,tr! INTENT(IN) ::
      REAL*8 rul(4),rur(4),flx(4)!,vf(3) ! INTENT(OUT) ::

! ==============================================================================
! Locals
! ==============================================================================

      REAL*8 af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,pf,ql,qr,vml,vmr,
     >        wtl,wtr,Hl,Hr

      real*8       rgam,rspec,rcp,rcv,k
      common /rgasprops/ rgam,rspec,rcp,rcv,k

      rl   = rul(1)
      ul   = rul(2)/rul(1)
      vl   = rul(3)/rul(1)
      retl = rul(4)/rul(1)

      rr   = rur(1)
      ur   = rur(2)/rur(1)
      vr   = rur(3)/rur(1)
      retr = rur(4)/rur(1)

      pl = MixtPerf_P(rul(1),rul(2),rul(3),rul(4))
      pr = MixtPerf_P(rur(1),rur(2),rur(3),rur(4))
c     pl = rl*(rgam-1.)*(retl - (ul*ul + vl*vl)/2.)
c     pr = rr*(rgam-1.)*(retr - (ur*ur + vr*vr)/2.)

      tl = MixtPerf_T(rul(1),rul(2),rul(3),rul(4),pl)
      tr = MixtPerf_T(rur(1),rur(2),rur(3),rur(4),pr)
c     tl = pl/rl/rspec
c     tr = pr/rr/rspec

      al = MixtPerf_A(rul(1),rul(2),rul(3),rul(4),tl)
      ar = MixtPerf_A(rur(1),rur(2),rur(3),rur(4),tr)
c     al = sqrt(rgam*rspec*tl)
c     ar = sqrt(rgam*rspec*tr)

      hl = MixtPerf_Ht(rul(1),rul(2),rul(3),rul(4),pl)
      hr = MixtPerf_Ht(rur(1),rur(2),rur(3),rur(4),pr)
c     hl = rcp*tl + 0.5*(ul**2 + vl**2)
c     hr = rcp*tr + 0.5*(ur**2 + vr**2)
         
c        ql = ul*nx + vl*ny + wl*nz - fs
c        qr = ur*nx + vr*ny + wr*nz - fs
         ql = ul*nx + vl*ny 
         qr = ur*nx + vr*ny 

         af = 0.5*(al+ar) ! NOTE not using original formulation, see note
         ml  = ql/af
         mla = ABS(ml)

         mr  = qr/af
         mra = ABS(mr)    

         IF ( mla .le. 1.0 ) THEN 
            mlp = 0.25*(ml+1.0)*(ml+1.0) + 0.125*(ml*ml-1.0)*(ml*ml-1.0)
            wtl = 0.25*(ml+1.0)*(ml+1.0)*(2.0-ml) +
     >            0.1875*ml*(ml*ml-1.0)*(ml*ml-1.0)
         ELSE
            mlp = 0.5*(ml+mla)
            wtl = 0.5*(1.0+ml/mla)
         END IF ! mla

         IF ( mra .le. 1.0 ) THEN 
            mrm = -0.25*(mr-1.0)*(mr-1.0)-0.125*(mr*mr-1.0)*(mr*mr-1.0)
            wtr = 0.25*(mr-1.0)*(mr-1.0)*(2.0+mr) -
     >            0.1875*mr*(mr*mr-1.0)*(mr*mr-1.0)
         ELSE
            mrm = 0.5*(mr-mra)
            wtr = 0.5*(1.0-mr/mra)
         END IF ! mla

         mf  = mlp + mrm
         mfa = ABS(mf)
         mfp = 0.5*(mf+mfa)
         mfm = 0.5*(mf-mfa)

         pf = wtl*pl + wtr*pr

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

!        vf(1) = mfp*ul + mfm*ur ! I'm sure we'll need this someday
!        vf(2) = mfp*vl + mfm*vr
!        vf(3) = mfp*wl + mfm*wr

         flx(1)=(af*(mfp*rl      +mfm*rr   )        )
         flx(2)=(af*(mfp*rl*ul+mfm*rr*ur)+pf*nx)
         flx(3)=(af*(mfp*rl*vl+mfm*rr*vr)+pf*ny)
c        flx(4)=(af*(mfp*rl*wl+mfm*rr*wr)+pf*nz)
         flx(4)=(af*(mfp*rl*Hl   +mfm*rr*Hr))
      return
      END
