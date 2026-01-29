c ---------------------------------------------------------------------
      program setup
c
      include 'MYDATA'

      call rzero(elnorm,11*lm)
      call rzero(nodel,2*lm)
      call izero(elml,9*lm)
      call izero(edgl,6*lm)
      call izero(edgb,6*lm)

c     -------------------- setup routines -----------------------------
      call read_user_input_file
      call read_grid
      call sort_grid
      call write_new_grid
c     -------------------- end setup routines -------------------------

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
