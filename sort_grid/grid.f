c ---------------------------------------------------------------------
      subroutine sort_grid
c
      include 'MYDATA'

      logical ifadded
      integer in1,in2,jn1,jn2,kn1,kn2
      integer n_bdry, n_bdry_unset, n_bdry_other, n_int_with_bc
      integer bc1, bc2, bc3

c     compute edges
      ned = 0
      do i=1,nel
         in1 = elml(1,i)
         in2 = elml(2,i)

         jn1 = elml(2,i)
         jn2 = elml(3,i)

         kn1 = elml(3,i)
         kn2 = elml(1,i)

         if (elml(9,i) .eq. 4) then
            kn1 = elml(3,i)
            kn2 = elml(4,i)

            nn1 = elml(4,i)
            nn2 = elml(1,i)
         endif


         if (i .eq. 1) then
            call add_edge(5,in1,in2,ned,i)
            call add_edge(6,jn1,jn2,ned,i)
            call add_edge(7,kn1,kn2,ned,i)
            if (elml(9,i) .eq. 4) call add_edge(8,nn1,nn2,ned,i)
         else
            call check_if_added(in1,in2,ned,ifadded,je) 
            if (ifadded)      call add_edge_element_only(5,je,i)
            if (.not.ifadded) call add_edge(5,in1,in2,ned,i)

            call check_if_added(jn1,jn2,ned,ifadded,je) 
            if (ifadded)      call add_edge_element_only(6,je,i)
            if (.not.ifadded) call add_edge(6,jn1,jn2,ned,i)

            call check_if_added(kn1,kn2,ned,ifadded,je) 
            if (ifadded)      call add_edge_element_only(7,je,i)
            if (.not.ifadded) call add_edge(7,kn1,kn2,ned,i)

            if (elml(9,i) .eq. 4) then
               call check_if_added(nn1,nn2,ned,ifadded,je) 
               if (ifadded)      call add_edge_element_only(8,je,i)
               if (.not.ifadded) call add_edge(8,nn1,nn2,ned,i)
            endif
         endif
         if (verbosity .ge. 2) then
            write(6,*) 'Sorted element: ',i,nel
         endif
      enddo

c     compute centroids
      do i=1,nel
        inode1 = elml(1,i)
        inode2 = elml(2,i)
        inode3 = elml(3,i)

        if (elml(9,i) .eq. 4) then
           inode4 = elml(4,i)
           rdum4x = nodel(1,inode4)
           rdum4y = nodel(2,inode4)
        else
           rdum4x = 0.
           rdum4y = 0.
        endif

        elnorm(1,i)=(nodel(1,inode1)+nodel(1,inode2)+
     >               nodel(1,inode3)+rdum4x)/elml(9,i)
        elnorm(2,i)=(nodel(2,inode1)+nodel(2,inode2)+
     >               nodel(2,inode3)+rdum4y)/elml(9,i)

         if (verbosity .ge. 2) then
            write(6,*) 'Centroid element: ',i,nel
         endif
      enddo

c     compute normals to edges
      do i=1,nel
         in1 = elml(1,i)
         in2 = elml(2,i)

         jn1 = elml(2,i)
         jn2 = elml(3,i)

         kn1 = elml(3,i)
         kn2 = elml(1,i)

         if (elml(9,i) .eq. 4) then
            kn1 = elml(3,i)
            kn2 = elml(4,i)
            
            nn1 = elml(4,i)
            nn2 = elml(1,i)
         endif

c        for first edge i
         rxc = (nodel(1,in1)-nodel(1,in2))
         ryc = (nodel(2,in1)-nodel(2,in2))
         rmx = ryc
         rmy = rxc
         rlen = sqrt(rxc**2 + ryc**2)

         xx = (nodel(1,in1)+nodel(1,in2))/2.
         yy = (nodel(2,in1)+nodel(2,in2))/2.

         rs = 1E-6
         xx1= xx + (-rmx/rlen)*rs
         yy1= yy + (rmy/rlen)*rs

         xx2= xx + (rmx/rlen)*rs
         yy2= yy + (-rmy/rlen)*rs

         rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
         rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

         if (rd1 .gt. rd2) then
            elnorm(4,i) = -rmx/rlen
            elnorm(5,i) = rmy/rlen
         else
            elnorm(4,i) =   rmx/rlen
            elnorm(5,i) = -rmy/rlen
         endif



c        for second edge j
         rxc = (nodel(1,jn1)-nodel(1,jn2))
         ryc = (nodel(2,jn1)-nodel(2,jn2))
         rmx = ryc
         rmy = rxc
         rlen = sqrt(rxc**2 + ryc**2)

         xx = (nodel(1,jn1)+nodel(1,jn2))/2.
         yy = (nodel(2,jn1)+nodel(2,jn2))/2.

         rs = 1E-6
         xx1= xx + (-rmx/rlen)*rs
         yy1= yy + (rmy/rlen)*rs

         xx2= xx + (rmx/rlen)*rs
         yy2= yy + (-rmy/rlen)*rs

         rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
         rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

         if (rd1 .gt. rd2) then
            elnorm(6,i) = -rmx/rlen
            elnorm(7,i) = rmy/rlen
         else
            elnorm(6,i) =   rmx/rlen
            elnorm(7,i) = -rmy/rlen
         endif


c        for third edge k
         rxc = (nodel(1,kn1)-nodel(1,kn2))
         ryc = (nodel(2,kn1)-nodel(2,kn2))
         rmx = ryc
         rmy = rxc
         rlen = sqrt(rxc**2 + ryc**2)

         xx = (nodel(1,kn1)+nodel(1,kn2))/2.
         yy = (nodel(2,kn1)+nodel(2,kn2))/2.

         rs = 1E-6
         xx1= xx + (-rmx/rlen)*rs
         yy1= yy + (rmy/rlen)*rs

         xx2= xx + (rmx/rlen)*rs
         yy2= yy + (-rmy/rlen)*rs

         rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
         rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

         if (rd1 .gt. rd2) then
            elnorm(8,i) = -rmx/rlen
            elnorm(9,i) = rmy/rlen
         else
            elnorm(8,i) =   rmx/rlen
            elnorm(9,i) = -rmy/rlen
         endif


         if (elml(9,i) .eq. 4) then
c           for fourth edge n
         rxc = (nodel(1,nn1)-nodel(1,nn2))
         ryc = (nodel(2,nn1)-nodel(2,nn2))
         rmx = ryc
         rmy = rxc
         rlen = sqrt(rxc**2 + ryc**2)

         xx = (nodel(1,nn1)+nodel(1,nn2))/2.
         yy = (nodel(2,nn1)+nodel(2,nn2))/2.

         rs = 1E-6
         xx1= xx + (-rmx/rlen)*rs
         yy1= yy + (rmy/rlen)*rs

         xx2= xx + (rmx/rlen)*rs
         yy2= yy + (-rmy/rlen)*rs

         rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
         rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

         if (rd1 .gt. rd2) then
            elnorm(10,i) = -rmx/rlen
            elnorm(11,i) = rmy/rlen
         else
            elnorm(10,i) =   rmx/rlen
            elnorm(11,i) = -rmy/rlen
         endif

         endif

         if (verbosity .ge. 2) then
            write(6,*) 'Normals element: ',i,nel
         endif
      enddo


c     connect edges to elements
      do i=1,nel
         in1 = elml(1,i)
         in2 = elml(2,i)

         jn1 = elml(2,i)
         jn2 = elml(3,i)

         kn1 = elml(3,i)
         kn2 = elml(1,i)

         if (elml(9,i) .eq. 4) then
            kn1 = elml(3,i)
            kn2 = elml(4,i)
          
            nn1 = elml(4,i)
            nn2 = elml(1,i)
         endif

         do j=1,ned ! now check edges
            in1e = edgl(1,j)
            in2e = edgl(2,j)

            if (in1 .eq. in1e) then ! first edge
            if (in2 .eq. in2e) then
               elml(5,i) = j
            endif
            endif
            if (in2 .eq. in1e) then
            if (in1 .eq. in2e) then
               elml(5,i) = j
            endif
            endif

            if (jn1 .eq. in1e) then ! second edge
            if (jn2 .eq. in2e) then
               elml(6,i) = j
            endif
            endif
            if (jn2 .eq. in1e) then
            if (jn1 .eq. in2e) then
               elml(6,i) = j
            endif
            endif

            if (kn1 .eq. in1e) then ! third edge
            if (kn2 .eq. in2e) then
               elml(7,i) = j
            endif
            endif
            if (kn2 .eq. in1e) then
            if (kn1 .eq. in2e) then
               elml(7,i) = j
            endif
            endif

            if (elml(9,i) .eq. 4) then
            if (nn1 .eq. in1e) then ! fourth edge
            if (nn2 .eq. in2e) then
               elml(8,i) = j
            endif
            endif
            if (nn2 .eq. in1e) then
            if (nn1 .eq. in2e) then
               elml(8,i) = j
            endif
            endif
            endif
         enddo

         if (verbosity .ge. 2) then
            write(6,*) 'Tagged element: ',i,nel
         endif
      enddo

c     compute areas of elements
      do ie=1,nel

         if (elml(9,ie) .eq. 4) then
            rxd1 = nodel(1,elml(1,ie)) - nodel(1,elml(3,ie))
            rxd2 = nodel(1,elml(2,ie)) - nodel(1,elml(4,ie))
            ryd1 = nodel(2,elml(1,ie)) - nodel(2,elml(3,ie))
            ryd2 = nodel(2,elml(2,ie)) - nodel(2,elml(4,ie))
            rd1 = sqrt(rxd1**2 + ryd1**2)
            rd2 = sqrt(rxd2**2 + ryd2**2)
            rxd1 = rxd1/rd1
            rxd2 = rxd2/rd2
            ryd1 = ryd1/rd1
            ryd2 = ryd2/rd2
            rth = acos(rxd1*rxd2+ryd1*ryd2)
            rarea = rd1*rd2*sin(rth)/2.
          
            elnorm(3,ie) = rarea
         else
            rxa = nodel(1,elml(1,ie))
            rya = nodel(2,elml(1,ie))
            rxb = nodel(1,elml(2,ie))
            ryb = nodel(2,elml(2,ie))
            rxc = nodel(1,elml(3,ie))
            ryc = nodel(2,elml(3,ie))

            elnorm(3,ie) = 0.5*((rxa-rxc)*(ryb-rya)-(rxa-rxb)*(ryc-rya))
            elnorm(3,ie) = abs(elnorm(3,ie))
         endif
         if (verbosity .ge. 2) then
            write(6,*) 'Area element: ',ie,nel
         endif
      enddo

c     sanity checks on boundary tagging
      n_bdry = 0
      n_bdry_unset = 0
      n_bdry_other = 0
      n_int_with_bc = 0
      bc1 = 0
      bc2 = 0
      bc3 = 0
      do i=1,ned
         if (edgl(4,i) .eq. 0) then
            n_bdry = n_bdry + 1
            if (edgl(5,i) .eq. 1) then
               bc1 = bc1 + 1
            elseif (edgl(5,i) .eq. 2) then
               bc2 = bc2 + 1
            elseif (edgl(5,i) .eq. 3) then
               bc3 = bc3 + 1
            elseif (edgl(5,i) .eq. 0) then
               n_bdry_unset = n_bdry_unset + 1
            else
               n_bdry_other = n_bdry_other + 1
            endif
         else
            if (edgl(5,i) .ne. 0) n_int_with_bc = n_int_with_bc + 1
         endif
      enddo

      if (verbosity .ge. 1) then
         write(6,*) 'Boundary edge counts (bc=1/2/3/0/other):',
     >        bc1, bc2, bc3, n_bdry_unset, n_bdry_other
      endif
      if (n_bdry_unset .gt. 0) then
         write(6,*) 'Warning: boundary edges with bc=0 (unassigned):',
     >        n_bdry_unset
      endif
      if (n_bdry_other .gt. 0) then
         write(6,*) 'Warning: boundary edges with unknown bc tag:',
     >        n_bdry_other
      endif
      if (n_int_with_bc .gt. 0) then
         write(6,*) 'Warning: interior edges with bc tag:',
     >        n_int_with_bc
      endif

c     do i=1,ned
c        if (edgl(5,i).eq. 0) then

c        iel1 = edgl(3,i)
c        iel2 = edgl(4,i)

c        nsd1 = elml(9,iel1)
c        do is = 0,nsd1-1
c           if (i .eq. elml(5+is,iel1)) is1 = is
c        enddo

c        nsd2 = elml(9,iel2)
c        do is = 0,nsd2-1
c           if (i .eq. elml(5+is,iel2)) is2 = is
c        enddo

c        rnx1 = elnorm(4+is1*2,iel1)
c        rny1 = elnorm(5+is1*2,iel1)

c        rnx2 = elnorm(4+is2*2,iel2)
c        rny2 = elnorm(5+is2*2,iel2)

c        rdot = rnx1*rnx2 + rny1*rny2


c        if (rdot .gt. 0) then
c           elnorm(4+is2*2,iel2) = -rnx2
c           elnorm(5+is2*2,iel2) = -rny2
c        endif

c        endif
c     enddo




      return
      end
c ---------------------------------------------------------------------
      subroutine check_if_added(in1,in2,ntmp,ifadded,j)
c
      include 'MYDATA'
      logical ifadded

      ifadded = .false.
      do j=1,ntmp
         if (edgl(1,j) .eq. in1) then
         if (edgl(2,j) .eq. in2) then
            ifadded = .true.
            goto 1511
         endif
         endif
         if (edgl(1,j) .eq. in2) then
         if (edgl(2,j) .eq. in1) then
            ifadded = .true.
            goto 1511
         endif
         endif
      enddo

 1511 continue
      return
      end
c ---------------------------------------------------------------------
      subroutine add_edge_element_only(neled,je,i)
c
      include 'MYDATA'

       elml(neled,i) = je
       edgl(4,je) = i

      return
      end
c ---------------------------------------------------------------------
      subroutine add_edge(neled,in1,in2,ntmp,i)
c
      include 'MYDATA'

      ntmp = ntmp + 1
      elml(neled,i) = ntmp

      edgl(1,ntmp) = in1
      edgl(2,ntmp) = in2
      edgl(3,ntmp) = i

      ! check if a boundary edge
      edgl(5,ntmp) = 0.
      do j=1,nbc
         if (edgb(1,j) .eq. in1) then
         if (edgb(2,j) .eq. in2) then
            edgl(5,ntmp) = edgb(5,j)
            goto 1511
         endif
         endif
         if (edgb(1,j) .eq. in2) then
         if (edgb(2,j) .eq. in1) then
            edgl(5,ntmp) = edgb(5,j)
            goto 1511
         endif
         endif
      enddo
       
 1511 continue
      return
      end
c ---------------------------------------------------------------------
      subroutine write_new_grid
c
      include 'MYDATA'


      if (write_vtk_mesh .ne. 0) then
         if (verbosity .ge. 1) then
            write(6,*) 'Writing grid vtk file'
         endif
         call write_grid_vtk
         if (verbosity .ge. 1) then
            write(6,*) 'Finished writing grid vtk file'
         endif
      endif

      if (verbosity .ge. 1) then
         write(6,*) 'Writing nodel file'
      endif
      open(unit=76,file='nodel')
      write(76,*) nn
      do i = 1,2*nn
         write(76,*) nodel(i,1)
      enddo
      close(76)
      if (verbosity .ge. 1) then
         write(6,*) 'Finished writing nodel file'
         write(6,*) 'Writing elnorm file'
      endif
      open(unit=76,file='elnorm')
      write(76,*) nel
      do i = 1,11*nel
         write(76,*) elnorm(i,1)
      enddo
      close(76)
      if (verbosity .ge. 1) then
         write(6,*) 'Finished writing elnorm file'
         write(6,*) 'Writing elml file'
      endif
      open(unit=76,file='elml')
      write(76,*) nel
      do i = 1,9*nel
         write(76,*) elml(i,1)
      enddo
      close(76)
      if (verbosity .ge. 1) then
         write(6,*) 'Finished writing elml file'
         write(6,*) 'Writing edgl file'
      endif
      open(unit=76,file='edgl')
      write(76,*) ned
      do i = 1,6*ned
         write(76,*) edgl(i,1)
      enddo
      close(76)
      if (verbosity .ge. 1) then
         write(6,*) 'Finished writing edgl file'
      endif

      return
      end
c ---------------------------------------------------------------------
      subroutine read_grid
c
      include 'MYDATA'

      character*64 dum_str, phys_name
      character*1024 line
      real*8 dum_real, mesh_version
      real*8 rx, ry, rz, rpar1, rpar2, rpar3
      integer dum_int, dum_tag
      integer ios
      integer npoints, ncurves, nsurfs, nvols
      integer nblocks, nblock, minTag, maxTag
      integer ent_dim, ent_tag, parametric
      integer elem_tag, elem_type, nnodes
      integer idx, node_tag
      integer gmsh_elem_nnodes
      integer block_tags(lm)
      integer nodemap(lm)
      integer curve_phys(lm)

      open(unit=82,file=gmsh_filename,form="formatted",
     >     status="old",action="read",iostat=ios)
      if (ios .ne. 0) then
         write(6,*) 'Error: cannot open mesh file: ', gmsh_filename
         stop
      endif

      ! $MeshFormat
      read (82,*) dum_str
      read (82,*) mesh_version, dum_int, dum_int
      read (82,*) dum_str

      ! $PhysicalNames (optional)
      read (82,*) dum_str
      ntags = 0
      if (dum_str .eq. '$PhysicalNames') then
         read (82,*) ntags
         do i=1,ntags
            read (82,*) ntags_dim(i), ntags_num(i), phys_name

            ntags_type(i) = 0                               ! no boundary
            if ( phys_name(1:4) .eq. 'wall' .or.
     >           phys_name(1:4) .eq. 'bump' .or.
     >           phys_name(1:4) .eq. 'top_' .or.
     >           phys_name(1:4) .eq. 'symm') then
               ntags_type(i) = 1
            elseif ( phys_name(1:4) .eq. 'infl' .or.
     >               phys_name(1:4) .eq. 'inle') then
               ntags_type(i) = 2
            elseif ( phys_name(1:4) .eq. 'outf' .or.
     >               phys_name(1:4) .eq. 'outl') then
               ntags_type(i) = 3
            endif
         enddo
         read (82,*) dum_str ! $EndPhysicalNames
         read (82,*) dum_str ! next section
      endif

      if (mesh_version .ge. 4.0d0) then
         do i=1,lm
            nodemap(i) = 0
            curve_phys(i) = 0
         enddo

         if (dum_str .ne. '$Entities') then
            write(6,*) 'Error: expected $Entities, found: ', dum_str
            stop
         endif
         read (82,*) npoints, ncurves, nsurfs, nvols

         do i=1,npoints
            read (82,'(A)') line
            call parse_entity_line(line, 0, ent_tag, dum_int)
         enddo
         do i=1,ncurves
            read (82,'(A)') line
            call parse_entity_line(line, 1, ent_tag, dum_int)
            if (ent_tag .ge. 1 .and. ent_tag .le. lm) then
               curve_phys(ent_tag) = dum_int
            endif
         enddo
         do i=1,nsurfs
            read (82,'(A)') line
            call parse_entity_line(line, 2, ent_tag, dum_int)
         enddo
         do i=1,nvols
            read (82,'(A)') line
            call parse_entity_line(line, 3, ent_tag, dum_int)
         enddo
         read (82,*) dum_str ! $EndEntities

         ! $Nodes
         read (82,*) dum_str
         if (dum_str .ne. '$Nodes') then
            write(6,*) 'Error: expected $Nodes, found: ', dum_str
            stop
         endif
         read (82,*) nblocks, nn, minTag, maxTag
         if (nn .gt. lm) then
            write(6,*) 'Error: too many nodes: ', nn, ' > ', lm
            stop
         endif

         idx = 0
         do ib=1,nblocks
            read (82,*) ent_dim, ent_tag, parametric, nblock
            if (nblock .gt. lm) then
               write(6,*) 'Error: node block too large: ', nblock
               stop
            endif
            call read_int_list(82, block_tags, nblock)
            do i=1,nblock
               if (parametric .eq. 0) then
                  read (82,*) rx, ry, rz
               else
                  if (ent_dim .eq. 1) then
                     read (82,*) rx, ry, rz, rpar1
                  elseif (ent_dim .eq. 2) then
                     read (82,*) rx, ry, rz, rpar1, rpar2
                  else
                     read (82,*) rx, ry, rz, rpar1, rpar2, rpar3
                  endif
               endif
               idx = idx + 1
               node_tag = block_tags(i)
               if (node_tag .gt. lm) then
                  write(6,*) 'Error: node tag out of range: ', node_tag
                  stop
               endif
               nodemap(node_tag) = idx
               nodel(1,idx) = rx
               nodel(2,idx) = ry
            enddo
         enddo
         read (82,*) dum_str ! $EndNodes

         ! $Elements
         read (82,*) dum_str
         if (dum_str .ne. '$Elements') then
            write(6,*) 'Error: expected $Elements, found: ', dum_str
            stop
         endif
         read (82,*) nblocks, dum_int, minTag, maxTag
         nbc = 0
         nel = 0
         do ib=1,nblocks
            read (82,*) ent_dim, ent_tag, elem_type, nblock
            nnodes = gmsh_elem_nnodes(elem_type)
            if (nnodes .le. 0) then
               write(6,*) 'Error: unsupported element type: ',
     >                    elem_type
               stop
            endif

            if (nnodes .eq. 1) then
               do i=1,nblock
                  read (82,*) elem_tag, node_tag
               enddo
            elseif (nnodes .eq. 2) then
               do i=1,nblock
                  read (82,*) elem_tag, in1, in2
                  nbc = nbc + 1
                  if (nodemap(in1) .eq. 0 .or. nodemap(in2) .eq. 0) then
                     write(6,*) 'Error: element node tag not found'
                     stop
                  endif
                  edgb(1,nbc) = nodemap(in1)
                  edgb(2,nbc) = nodemap(in2)
                  edgb(5,nbc) = curve_phys(ent_tag)
                  ndum1 = edgb(5,nbc)
                  do j=1,ntags ! mark boundary condition
                     if (ndum1 .eq. ntags_num(j)) then
                        edgb(5,nbc)=ntags_type(j) 
                     endif
                  enddo
               enddo
            elseif (nnodes .eq. 3) then
               do i=1,nblock
                  read (82,*) elem_tag, in1, in2, in3
                  nel = nel + 1
                  elml(1,nel) = nodemap(in1)
                  elml(2,nel) = nodemap(in2)
                  elml(3,nel) = nodemap(in3)
                  elml(9,nel) = 3 !number of sides
               enddo
            elseif (nnodes .eq. 4) then
               do i=1,nblock
                  read (82,*) elem_tag, in1, in2, in3, in4
                  nel = nel + 1
                  elml(1,nel) = nodemap(in1)
                  elml(2,nel) = nodemap(in2)
                  elml(3,nel) = nodemap(in3)
                  elml(4,nel) = nodemap(in4)
                  elml(9,nel) = 4 !number of sides
               enddo
            endif
         enddo
         read (82,*) dum_str ! $EndElements
         close(82)

         nmm = neq*nel ! matrix actual size
         return
      endif

      ! ---- MSH 2.x legacy format ----
      if (dum_str .ne. '$Nodes') then
         write(6,*) 'Error: expected $Nodes, found: ', dum_str
         stop
      endif
      read (82,*) nn
      do i=1,nn
         read (82,*) dum_int, nodel(1,i), nodel(2,i), dum_real
      enddo
      read (82,*) dum_str

      nbc = 0
      nel = 0
      if (dum_str .ne. '$EndNodes') then
         write(6,*) 'Error: expected $EndNodes, found: ', dum_str
         stop
      endif
      read (82,*) dum_str
      if (dum_str .ne. '$Elements') then
         write(6,*) 'Error: expected $Elements, found: ', dum_str
         stop
      endif
      read (82,*) nelt
      do i=1,nelt
         read (82,*) dum_int,dum_int,dum_tag

         if (dum_int .eq. 1) then     ! line element
            ndum = 2
            nbc = nbc + 1
            backspace 82
            ! assume only 2 tags for lines ...
            read (82,*) dum_int,dum_int,dum_int,edgb(5,nbc),dum_int
     >                 ,edgb(1,nbc),edgb(2,nbc)
            ndum1 = edgb(5,nbc)
            do j=1,ntags ! mark boundary condition
               if (ndum1 .eq. ntags_num(j)) then
                  edgb(5,nbc)=ntags_type(j) 
               endif
            enddo

         elseif (dum_int .eq. 2) then ! triangle
            ndum = 3
            nel = nel + 1
            backspace 82
            ! assume only 2 tags for triangles ...?
            read (82,*) dum_int,dum_int,dum_int,dum_int,dum_int
     >                 ,elml(1,nel),elml(2,nel),elml(3,nel)
            elml(9,nel) = 3 !number of sides
         elseif (dum_int .eq. 3) then ! quad
            ndum = 4
            nel = nel + 1
            backspace 82
            read (82,*) dum_int,dum_int,dum_int,dum_int,dum_int
     >               ,elml(1,nel),elml(2,nel),elml(3,nel),elml(4,nel)
            elml(9,nel) = 4 !number of sides
         endif
      enddo
      read (82,*) dum_str
      close(82)

      nmm = neq*nel ! matrix actual size

      return
      end
c ---------------------------------------------------------------------
      integer function gmsh_elem_nnodes(elem_type)
c
      integer elem_type

      if (elem_type .eq. 1) then
         gmsh_elem_nnodes = 2
      elseif (elem_type .eq. 2) then
         gmsh_elem_nnodes = 3
      elseif (elem_type .eq. 3) then
         gmsh_elem_nnodes = 4
      elseif (elem_type .eq. 15) then
         gmsh_elem_nnodes = 1
      else
         gmsh_elem_nnodes = 0
      endif

      return
      end
c ---------------------------------------------------------------------
      subroutine parse_entity_line(line, dim, ent_tag, phys_tag)
c
      character*(*) line
      integer dim, ent_tag, phys_tag
      integer num_phys, num_bnd
      integer ntok, p

      character*32 tokens(512)

      call tokenize(line, tokens, ntok)
      p = 1
      if (ntok .le. 0) goto 900

      read(tokens(p),*) ent_tag
      p = p + 1

      if (dim .eq. 0) then
         p = p + 3
      else
         p = p + 6
      endif
      if (p .gt. ntok) goto 900

      read(tokens(p),*) num_phys
      p = p + 1

      phys_tag = 0
      if (num_phys .gt. 0) then
         read(tokens(p),*) phys_tag
      endif
      p = p + num_phys
      if (p .gt. ntok) return

      read(tokens(p),*) num_bnd
      p = p + 1 + num_bnd

      return

 900  continue
      write(6,*) 'Error: failed to parse $Entities line:'
      write(6,*) line
      stop
      end
c ---------------------------------------------------------------------
      subroutine tokenize(line, tokens, ntok)
c
      character*(*) line
      character*32 tokens(*)
      integer ntok
      integer i, linelen, start, tlen
      integer max_tok
      parameter (max_tok=512)

      ntok = 0
      linelen = len(line)
      do while (linelen .gt. 0)
         if (line(linelen:linelen) .ne. ' ' .and.
     >       line(linelen:linelen) .ne. char(9)) goto 10
         linelen = linelen - 1
      enddo
 10   continue

      i = 1
      do while (i .le. linelen)
         do while (i .le. linelen .and.
     >            (line(i:i) .eq. ' ' .or. line(i:i) .eq. char(9)))
            i = i + 1
         enddo
         if (i .gt. linelen) goto 20

         start = i
         do while (i .le. linelen .and.
     >            (line(i:i) .ne. ' ' .and. line(i:i) .ne. char(9)))
            i = i + 1
         enddo

         tlen = i - start
         ntok = ntok + 1
         if (ntok .gt. max_tok) then
            write(6,*) 'Error: too many tokens in entity line'
            stop
         endif
         tokens(ntok) = ' '
         if (tlen .gt. 32) tlen = 32
         tokens(ntok)(1:tlen) = line(start:start+tlen-1)
      enddo
 20   continue

      return
      end
c ---------------------------------------------------------------------
      subroutine read_int_list(unit, arr, n)
c
      integer unit, n
      integer arr(n)
      integer idx, i, ntok
      character*1024 line
      character*32 tokens(512)

      idx = 0
      do while (idx .lt. n)
         read(unit,'(A)') line
         call tokenize(line, tokens, ntok)
         do i=1,ntok
            idx = idx + 1
            read(tokens(i),*) arr(idx)
            if (idx .ge. n) goto 20
         enddo
      enddo
 20   continue

      return
      end
c ---------------------------------------------------------------------
      subroutine write_grid_vtk
c
      include 'MYDATA'

      character*15 fnamec,fnamen

      integer icalld
      save    icalld
      data    icalld  /-1/


      icalld = icalld + 1

      write(fnamec,'(A6,I5.5,A4)') 'grid_c', icalld, '.vtk'
      open(unit=72,file=fnamec)

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

      write(72,*) 'CELLS ', nel, nel*5
      do i=1,nel
         if (elml(9,i) .eq. 3) then
            write(72,*) 
     >        elml(9,i),elml(1,i)-1,elml(2,i)-1,elml(3,i)-1
         elseif (elml(9,i) .eq. 4) then
            write(72,*) 
     >        elml(9,i),elml(1,i)-1,elml(2,i)-1,elml(3,i)-1,elml(4,i)-1
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
      write(72,*) 'SCALARS  dumdata  float 1'
      write(72,*) 'LOOKUP_TABLE default'
      do i=1,nel
         write(72,*) 0.
      enddo
      write(72,*) ''

      close(72)

      return
      end
c ---------------------------------------------------------------------
      subroutine read_user_input_file
c
      include 'MYDATA'

      integer ios

      open(unit=81,file="prj.ini",form="formatted")
      read (81,*) gmsh_filename
      write_vtk_mesh = 1
      verbosity = 1
      read (81,*,iostat=ios) write_vtk_mesh
      if (ios .eq. 0) then
         read (81,*,iostat=ios) verbosity
         if (ios .ne. 0) verbosity = 1
      endif
      close(81)

      return
      end
