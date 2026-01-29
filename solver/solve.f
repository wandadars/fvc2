c ---------------------------------------------------------------------
      subroutine zero_matricies(bdum,istep)
c
      include 'MYDATA'

c     real*8 bdum(nmm),Adum(nmm,nmm)
      real*8 bdum(nmm)

      do j=1,nmm
         bdum(j) = 0.
      enddo

c     ! need to zero everything once...
c     if (istep .eq. 0) then
c     do j=1,nmm
c     do i=1,nmm
c        Adum(i,j) = 0.
c     enddo
c     enddo
c     
c     ! but after, only zero where you will update...
c     else
c     
c     do ied=1,ned
c        iel1 = edgl(3,ied)

c        nsd1 = elml(9,iel1)
c        do is=0,nsd1-1
c            iel1_s = elml(5+is,iel1)

c            Adum(iel1,iel1_s) = 0.
c        enddo

c        if (edgl(5,ied).eq.0) then
c           iel2 = edgl(4,ied)
c           nsd2 = elml(9,iel2)
c           do is=0,nsd2-1
c               iel2_s = elml(5+is,iel2)
c          
c               Adum(iel2,iel2_s) = 0.
c           enddo
c        endif
c     enddo

c     endif
         

      return
      end

c ---------------------------------------------------------------------
      subroutine solve_equations_fwdE(bdum,istage)
c
      include 'MYDATA'

c     real*8 bdum(nmm),Adum(nmm,nmm)
      real*8 bdum(nmm)
      integer ipiv(nmm),info

      ! forward euler time scheme
      do j=1,nmm
         jel = int((j-1)/neq) + 1

         bdum(j) = bdum(j)/elnorm(3,jel)*dt
      enddo

c     set solution
      ic = 0
      do i=1,nel
         do ii=1,neq
            ic = ic + 1
            u(i,ii) = u(i,ii) + bdum(ic)
         enddo
      enddo

      return
      end
c ---------------------------------------------------------------------
      subroutine compute_sources(bdum)
c
      include 'MYDATA'

      real*8 bdum(nmm),rsource1(neq),rsource2(neq)

      ! user source term
      do j=1,neq
         rsource1(j) = 0.
         rsource2(j) = 0.
      enddo
      call user_source(elnorm(1,iel1),elnorm(2,iel1),rsource1)
      do i=0,neq-1
         iems1 = iem1 + i
         rsource1(i+1) = rsource1(i+1)/4. ! divide by number of edges since this
         bdum(iems1) = bdum(iems1) + rsource1(i+1)*elnorm(3,iel1) ! by area
      enddo

      if (edgl(5,ied).eq.0) then ! opposite element if exist
         call user_source(elnorm(1,iel2),elnorm(2,iel2),rsource2)
         do i=0,neq-1
            iems2 = iem2 + i
            rsource2(i+1) = rsource2(i+1)/4. !    ideally would be callled on elm basis
            bdum(iems2) = bdum(iems2) + rsource2(i+1)*elnorm(3,iel2) ! by area
         enddo
      endif

      return
      end
c ---------------------------------------------------------------------
      subroutine compute_edge_values
c
      include 'MYDATA'

      iel1 = edgl(3,ied)
      iem1 = 1 + neq*(iel1-1)

      iems1r  = iem1 + 0
      iems1ru = iem1 + 1
      iems1rv = iem1 + 2
      iems1re = iem1 + 3

      nsd1 = elml(9,iel1)
      do is=0,nsd1-1
         if (ied .eq. elml(5+is,iel1)) is1 = is
      enddo


      ! interior edges
      if (edgl(5,ied).eq.0) then


         ! remember, there is another element to connect to...
         iel2 = edgl(4,ied)
         iem2 = 1 + neq*(iel2-1)
       
         iems2r  = iem2 + 0
         iems2ru = iem2 + 1
         iems2rv = iem2 + 2
         iems2re = iem2 + 3
       
         nsd2 = elml(9,iel2)
         do is=0,nsd2-1
            if (ied .eq. elml(5+is,iel2)) is2 = is
         enddo

         rx1 = elnorm(1,iel1)
         ry1 = elnorm(2,iel1)
         rx2 = elnorm(1,iel2)
         ry2 = elnorm(2,iel2)
         dxi = sqrt((ry2-ry1)**2 + (rx2-rx1)**2)
         exi1(1) = (rx2-rx1)/dxi
         exi1(2) = (ry2-ry1)/dxi
         exi2(1) = -exi1(1)
         exi2(2) = -exi1(2)

         rx1 = nodel(1,edgl(1,ied))
         ry1 = nodel(2,edgl(1,ied))
         rx2 = nodel(1,edgl(2,ied))
         ry2 = nodel(2,edgl(2,ied))
         deta = sqrt((ry2-ry1)**2 + (rx2-rx1)**2)
         eeta(1) = (rx2-rx1)/deta
         eeta(2) = (ry2-ry1)/deta

         n_dot_xi1 = elnorm(4+is1*2,iel1)*exi1(1)
     >             + elnorm(5+is1*2,iel1)*exi1(2)
         n_dot_xi2 = elnorm(4+is2*2,iel2)*exi2(1)
     >             + elnorm(5+is2*2,iel2)*exi2(2)

      elseif (edgl(5,ied).ne.0) then ! boundaries


         rx1 = elnorm(1,iel1)
         ry1 = elnorm(2,iel1)
         rx2 = (nodel(1,edgl(1,ied)) +  ! just average face
     >          nodel(1,edgl(2,ied)))/2.
         ry2 = (nodel(2,edgl(1,ied)) + 
     >          nodel(2,edgl(2,ied)))/2.
         dxi = sqrt((ry2-ry1)**2 + (rx2-rx1)**2)
         exi1(1) = (rx2-rx1)/dxi
         exi1(2) = (ry2-ry1)/dxi
         exi2(1) = 0.
         exi2(2) = 0.

         rx1 = nodel(1,edgl(1,ied))
         ry1 = nodel(2,edgl(1,ied))
         rx2 = nodel(1,edgl(2,ied))
         ry2 = nodel(2,edgl(2,ied))
         deta = sqrt((ry2-ry1)**2 + (rx2-rx1)**2)
         eeta(1) = (rx2-rx1)/deta
         eeta(2) = (ry2-ry1)/deta

         n_dot_xi1 = elnorm(4+is1*2,iel1)*exi1(1)
     >             + elnorm(5+is1*2,iel1)*exi1(2)
         n_dot_xi2 = 0.

      endif

      return
      end
