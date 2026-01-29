c ---------------------------------------------------------------------
      subroutine write_soln
c
      include 'MYDATA'

      character(len=32) :: fnamec
      character(len=256) :: outdir, outfile
      character(len=32) :: nn_str, nel_str
      character(len=200) :: line

      integer icalld
      integer cells_size
      integer nbad_cfl
      real*8 cfl_max, cfl_local, rho, uvel, vvel, p, a, speed, area, dx
      save    icalld
      data    icalld  /-1/


      icalld = icalld + 1
      cfl_max = 0.d0
      nbad_cfl = 0
      do i=1,nel
         rho = u(i,1)
         if (rho .le. 0.d0) then
            nbad_cfl = nbad_cfl + 1
            goto 10
         endif
         uvel = u(i,2)/rho
         vvel = u(i,3)/rho
         p = (rgam-1.d0)*
     >       (u(i,4) - 0.5d0*(u(i,2)*u(i,2) + u(i,3)*u(i,3))/rho)
         if (p .le. 0.d0) then
            nbad_cfl = nbad_cfl + 1
            goto 10
         endif
         area = elnorm(3,i)
         if (area .le. 0.d0) then
            nbad_cfl = nbad_cfl + 1
            goto 10
         endif
         dx = sqrt(area)
         a = sqrt(rgam*p/rho)
         speed = sqrt(uvel*uvel + vvel*vvel)
         cfl_local = dt*(speed + a)/dx
         if (cfl_local .gt. cfl_max) cfl_max = cfl_local
 10      continue
      enddo

      write (6,'(A,I6,2X,ES12.6,2X,ES12.6,2X,A,ES12.6)') 'iostep: ',
     >    icalld,time,dt,'CFLmax=',cfl_max
      if (nbad_cfl .gt. 0) then
         write (6,'(A,I6)') '  CFL skipped cells: ', nbad_cfl
      endif

      call ensure_output_dir
      outdir = 'output/'

      cells_size = 0
      do i=1,nel
         cells_size = cells_size + elml(9,i) + 1
      enddo

      if (output_mode .eq. 1 .or. output_mode .eq. 2) then
         write(fnamec,'(A6,I5.5,A4)') 'soln_c', icalld, '.vtk'
         outfile = trim(outdir)//trim(fnamec)
         open(unit=72,file=outfile)

         write(72,'(A26)') '# vtk DataFile Version 3.0'
         write(72,*) '2D scalar data'
         write(72,*) 'ASCII'
         write(72,*) ''

         write(72,*) 'DATASET UNSTRUCTURED_GRID'
         write(72,*) 'POINTS',' ', nn,' ',' float'
         do i=1,nn
            write(72,*) nodel(1,i),nodel(2,i),0.
         enddo
         write(72,*) ''

         write(72,*) 'CELLS ', nel, cells_size
         do i=1,nel
            if (elml(9,i) .eq. 3) then
               write(72,*) 3,elml(1,i)-1,elml(2,i)-1,elml(3,i)-1
            elseif (elml(9,i) .eq. 4) then
               write(72,*) 4,elml(1,i)-1,elml(2,i)-1,
     >                    elml(3,i)-1,elml(4,i)-1
            endif
         enddo
         write(72,*) ''

         write(72,*) 'CELL_TYPES', nel
         do i=1,nel
            if (elml(9,i) .eq. 3) write(72,*) 5
            if (elml(9,i) .eq. 4) write(72,*) 9
         enddo
         write(72,*) ''

         write(72,*) 'CELL_DATA',nel
         write(72,*) 'SCALARS density float 1'
         write(72,*) 'LOOKUP_TABLE default'
         do i=1,nel
            write(72,*) u(i,1)
         enddo
         write(72,*) 'SCALARS ru float 1'
         write(72,*) 'LOOKUP_TABLE default'
         do i=1,nel
            write(72,*) u(i,2)
         enddo
         write(72,*) 'SCALARS rv float 1'
         write(72,*) 'LOOKUP_TABLE default'
         do i=1,nel
            write(72,*) u(i,3)
         enddo
         write(72,*) 'SCALARS re float 1'
         write(72,*) 'LOOKUP_TABLE default'
         do i=1,nel
            write(72,*) u(i,4)
         enddo
         write(72,*) ''

         close(72)
      endif

      if (output_mode .eq. 0 .or. output_mode .eq. 2) then
         write(fnamec,'(A6,I5.5,A4)') 'soln_c', icalld, '.dat'
         outfile = trim(outdir)//trim(fnamec)
         open(unit=72,file=outfile)
      
      write(72,*) 'TITLE = "Your Title Here"'
      !write(72,*) 'VARIABLES = "X", "Y", "Z", "Density", "RU", "RV", "RE"'
      !write(72,*) 'VARIABLES = ', '"X",', '"Y",', '"Z",', '"Density",',&
      !'"RU",', '"RV",', '"RE"'
      !write(72,*) 'VARIABLES = ' // '"X", "Y", "Z", "Density", "RU", "RV", "RE"'
      !write(72,*) "VARIABLES = 'X', 'Y', 'Z', 'Density', 'RU', 'RV', 'RE'"
      write(72,*) 'VARIABLES = "X", "Y", "Z", "Density", 
     & "RU", "RV", "RE"' 
      write(72,*) 'ZONE N=',nn,', E=',nel,', DATAPACKING=POINT, 
     & ZONETYPE=FETRIANGLE'
!      ! Write nodal coordinates
!      do i=1,nn
!         write(72,*) nodel(1,i),nodel(2,i),0.0
!      enddo
!      
!      ! Write connectivity
!      do i=1,nel
!         if (elml(9,i) .eq. 3) then
!            write(72,*) elml(1,i),elml(2,i),elml(3,i)
!         elseif (elml(9,i) .eq. 4) then
!            write(72,*) elml(1,i),elml(2,i),elml(3,i),elml(4,i)
!         endif
!      enddo
!      
!      ! Write solution data
!      do j=1,4
!         write(72,*) 'ZONE VARLOCATION=([4-7]=CELLCENTERED)'
!         do i=1,nel
!            write(72,*) u(i,j)
!         enddo
!      enddo

      ! Write nodal coordinates
      do i=1,nn
         write(72,*) nodel(1,i),nodel(2,i),0.0, u(i,1), u(i,2),
     &       u(i,3), u(i,4)
      enddo

      ! Write connectivity
      do i=1,nel
         if (elml(9,i) .eq. 3) then
            write(72,*) elml(1,i),elml(2,i),elml(3,i)
         elseif (elml(9,i) .eq. 4) then
            write(72,*) elml(1,i),elml(2,i),elml(3,i),elml(4,i)
         endif
      enddo

      ! Write solution data
      ! Assuming the solution is at the cell centers and u is a 2D array
      ! with nel rows (elements) and 4 columns (variables)
      !write(72,*) 'ZONE VARLOCATION=([4-7]=CELLCENTERED)'
      !do i=1,nel
      !   write(72,*) u(i,1), u(i,2), u(i,3), u(i,4)
      !enddo
      !
         close(72)
      endif
      

      return
      end
c ---------------------------------------------------------------------
      subroutine read_sorted_grid
c
      include 'MYDATA'

      character(len=256) :: mesh_dir, path
      logical use_grid
      integer ios

      mesh_dir = ''
      inquire(file='grid/nodel', exist=use_grid)
      if (use_grid) then
         mesh_dir = 'grid/'
      endif

      write(6,*) 'Reading nodel file'
      path = trim(mesh_dir)//'nodel'
      open(unit=76,file=path,status='old',action='read',iostat=ios)
      if (ios .ne. 0) then
         write(6,*) 'Error: cannot open mesh file: ', trim(path)
         stop
      endif
      read(76,*) nn
      do i = 1,nn
         read(76,*) nodel(1,i)
         read(76,*) nodel(2,i)
      enddo
      close(76)
      write(6,*) 'Finished reading nodel file'

      write(6,*) 'Reading elnorm file'
      path = trim(mesh_dir)//'elnorm'
      open(unit=76,file=path,status='old',action='read',iostat=ios)
      if (ios .ne. 0) then
         write(6,*) 'Error: cannot open mesh file: ', trim(path)
         stop
      endif
      read(76,*) nel
      do i = 1,nel
         read(76,*) elnorm(1,i)
         read(76,*) elnorm(2,i)
         read(76,*) elnorm(3,i)
         read(76,*) elnorm(4,i)
         read(76,*) elnorm(5,i)
         read(76,*) elnorm(6,i)
         read(76,*) elnorm(7,i)
         read(76,*) elnorm(8,i)
         read(76,*) elnorm(9,i)
         read(76,*) elnorm(10,i)
         read(76,*) elnorm(11,i)
      enddo
      close(76)
      write(6,*) 'Finished reading elnorm file'

      write(6,*) 'Reading elml file'
      path = trim(mesh_dir)//'elml'
      open(unit=76,file=path,status='old',action='read',iostat=ios)
      if (ios .ne. 0) then
         write(6,*) 'Error: cannot open mesh file: ', trim(path)
         stop
      endif
      read(76,*) nel
      do i = 1,nel
         read(76,*) elml(1,i)
         read(76,*) elml(2,i)
         read(76,*) elml(3,i)
         read(76,*) elml(4,i)
         read(76,*) elml(5,i)
         read(76,*) elml(6,i)
         read(76,*) elml(7,i)
         read(76,*) elml(8,i)
         read(76,*) elml(9,i)
      enddo
      close(76)
      write(6,*) 'Finished reading elml file'

      write(6,*) 'Reading edgl file'
      path = trim(mesh_dir)//'edgl'
      open(unit=76,file=path,status='old',action='read',iostat=ios)
      if (ios .ne. 0) then
         write(6,*) 'Error: cannot open mesh file: ', trim(path)
         stop
      endif
      read(76,*) ned
      do i = 1,ned
         read(76,*) edgl(1,i)
         read(76,*) edgl(2,i)
         read(76,*) edgl(3,i)
         read(76,*) edgl(4,i)
         read(76,*) edgl(5,i)
         read(76,*) edgl(6,i)
      enddo
      close(76)
      write(6,*) 'Finished reading edgl file'

      nmm = neq*nel ! matrix actual size

      return
      end
c ---------------------------------------------------------------------
      subroutine read_user_input_file
c
      include 'MYDATA'

      real*8 mmass
      integer ios
      character(len=256) :: line

      open(unit=81,file="prj.ini",form="formatted")
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) dt
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) iostep
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) nstepmax
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) rgam
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) mmass
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) k
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) time_scheme
      call read_next_line(81,line,ios)
      if (ios .ne. 0) goto 910
      read (line,*) conv_scheme
      do i=1,5 ! dirichlet
         call read_next_line(81,line,ios)
         if (ios .ne. 0) goto 910
         read (line,*) bcvals(1,i),bcvals(2,i),bcvals(3,i),bcvals(4,i)
      enddo
      do i=6,10 ! flux
         call read_next_line(81,line,ios)
         if (ios .ne. 0) goto 910
         read (line,*) bcvals(1,i),bcvals(2,i)
      enddo

      ! Initial condition defaults (used if IC block is missing)
      ic_type = 1
      ic_uniform(1) = 1.4
      ic_uniform(2) = 1.0
      ic_uniform(3) = 3.0
      ic_uniform(4) = 0.0

      ic_split_x = 0.5
      ic_left(1) = 1.0
      ic_left(2) = 1.0
      ic_left(3) = 0.0
      ic_left(4) = 0.0
      ic_right(1) = 0.125
      ic_right(2) = 0.1
      ic_right(3) = 0.0
      ic_right(4) = 0.0

      output_mode = 2

      ! Optional IC block (all lines required if ic_type is provided)
      call read_next_line(81,line,ios)
      if (ios .eq. 0) then
         read (line,*) ic_type

         call read_next_line(81,line,ios)
         if (ios .ne. 0) goto 920
         read (line,*) ic_uniform(1),ic_uniform(2),
     >        ic_uniform(3),ic_uniform(4)

         call read_next_line(81,line,ios)
         if (ios .ne. 0) goto 930
         read (line,*) ic_split_x

         call read_next_line(81,line,ios)
         if (ios .ne. 0) goto 940
         read (line,*) ic_left(1),ic_left(2),
     >        ic_left(3),ic_left(4)

         call read_next_line(81,line,ios)
         if (ios .ne. 0) goto 950
         read (line,*) ic_right(1),ic_right(2),
     >        ic_right(3),ic_right(4)
      endif

      call read_next_line(81,line,ios)
      if (ios .eq. 0) then
         read (line,*) output_mode
      endif
      close(81)


      rideal = 8.3144598
      rspec  = rideal/mmass
      rcp    = rspec/(1.-1./rgam)
      rcv    = rcp/rgam


      return
 910  continue
      write(6,*) 'Error: unexpected end of prj.ini'
      stop
 920  continue
      write(6,*) 'Error: expected ic_uniform line in prj.ini'
      stop
 930  continue
      write(6,*) 'Error: expected ic_split_x line in prj.ini'
      stop
 940  continue
      write(6,*) 'Error: expected ic_left line in prj.ini'
      stop
 950  continue
      write(6,*) 'Error: expected ic_right line in prj.ini'
      stop
      contains
      function writeToString(i)
         integer, intent(in) :: i
         character(len=32) :: writeToString
     
         write(writeToString, '(I32)') i
         writeToString = adjustl(writeToString)
      end function writeToString
     
      end
c ---------------------------------------------------------------------
      subroutine read_next_line(unit, line, ios)
c
      integer unit, ios
      character*(*) line
      integer p1, p2, cut

 10   continue
      read(unit,'(A)',iostat=ios) line
      if (ios .ne. 0) return

      p1 = index(line,'!')
      p2 = index(line,'#')
      cut = 0
      if (p1 .gt. 0) cut = p1
      if (p2 .gt. 0 .and. (cut .eq. 0 .or. p2 .lt. cut)) cut = p2
      if (cut .gt. 0) line = line(1:cut-1)

      if (len_trim(line) .eq. 0) goto 10

      return
      end
c ---------------------------------------------------------------------
      subroutine ensure_output_dir
c
      call system('mkdir -p output')
      return
      end
c ---------------------------------------------------------------------
