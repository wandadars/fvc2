module msh2fvc_data
  ! Global mesh state shared by all routines (modern replacement for COMMON blocks).
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: lm = 200000   ! max nodes/elements
  integer, parameter :: neq = 1       ! number equns solving
  integer, parameter :: lmm = lm * neq
  integer, parameter :: lbc = 20      ! total number of boundary conditions

  integer :: nel = 0, nn = 0, ned = 0
  integer :: nmm = 0

  real(real64) :: nodel(2, lm) = 0.0_real64
  real(real64) :: elnorm(11, lm) = 0.0_real64

  integer :: elml(9, lm) = 0
  integer :: edgl(6, lm) = 0
  integer :: edgb(6, lm) = 0

  integer :: ntags = 0
  integer :: ntags_dim(lbc) = 0
  integer :: ntags_num(lbc) = 0
  integer :: ntags_type(lbc) = 0
  character(len=64) :: ntags_name(lbc) = ''

  integer :: nbc = 0

  character(len=255) :: gmsh_filename = ''
  real(real64) :: coord_scale = 1.0_real64

  integer :: write_vtk_mesh = 0
  integer :: verbosity = 1
contains
  subroutine reset_arrays()
    ! Ensure all mesh buffers are cleared before reading a new grid.
    nodel = 0.0_real64
    elnorm = 0.0_real64
    elml = 0
    edgl = 0
    edgb = 0
  end subroutine reset_arrays
end module msh2fvc_data

module msh2fvc_ops
  use iso_fortran_env, only: real64, int64, output_unit, error_unit
  use msh2fvc_data
  implicit none
contains
  subroutine read_cmdline()
    ! Parse CLI args: mesh filename + optional unit scale flag.
    character(len=256) :: arg
    character(len=256) :: unit_arg
    integer :: nargs
    integer :: i
    logical :: have_mesh, have_units

    nargs = command_argument_count()
    if (nargs < 1 .or. nargs > 2) then
      write(output_unit,*) 'Usage: msh2fvc <meshfile> [-m|-mm|-ft|-in]'
      error stop
    end if

    coord_scale = 1.0_real64
    gmsh_filename = ''
    have_mesh = .false.
    have_units = .false.
    do i = 1, nargs
      call get_command_argument(i, arg)
      unit_arg = adjustl(arg)
      select case (trim(unit_arg))
      case ('-m')
        if (have_units) goto 900
        coord_scale = 1.0_real64
        have_units = .true.
      case ('-mm')
        if (have_units) goto 900
        coord_scale = 1.0e-3_real64
        have_units = .true.
      case ('-ft')
        if (have_units) goto 900
        coord_scale = 0.3048_real64
        have_units = .true.
      case ('-in')
        if (have_units) goto 900
        coord_scale = 0.0254_real64
        have_units = .true.
      case default
        if (have_mesh) goto 901
        gmsh_filename = trim(arg)
        have_mesh = .true.
      end select
    end do

    if (.not. have_mesh) then
      write(output_unit,*) 'Error: missing mesh file'
      write(output_unit,*) 'Usage: msh2fvc <meshfile> [-m|-mm|-ft|-in]'
      error stop
    end if

    write_vtk_mesh = 0
    verbosity = 1
    return

900 continue
    write(output_unit,*) 'Error: multiple units flags provided'
    write(output_unit,*) 'Usage: msh2fvc <meshfile> [-m|-mm|-ft|-in]'
    error stop
901 continue
    write(output_unit,*) 'Error: multiple mesh files provided'
    write(output_unit,*) 'Usage: msh2fvc <meshfile> [-m|-mm|-ft|-in]'
    error stop
  end subroutine read_cmdline

  subroutine sort_grid()
    ! Build edge list, centroids, normals, and element areas (edge lookup via hash).
    integer :: in1, in2, jn1, jn2, kn1, kn2, nn1, nn2
    integer :: n_bdry, n_bdry_unset, n_bdry_other, n_int_with_bc
    integer :: bc1, bc2, bc3
    integer :: i, ie
    integer :: inode1, inode2, inode3, inode4
    real(real64) :: rdum4x, rdum4y
    real(real64) :: rxc, ryc, rmx, rmy, rlen
    real(real64) :: xx, yy, rs, xx1, yy1, xx2, yy2
    real(real64) :: rd1, rd2
    real(real64) :: rxd1, rxd2, ryd1, ryd2, rth, rarea
    real(real64) :: rxa, rya, rxb, ryb, rxc_tri, ryc_tri
    integer, allocatable :: h_key1(:), h_key2(:), h_val(:)
    integer :: hsize, target

    ned = 0
    target = max(1024, 8 * max(1, nel))
    hsize = 1
    do while (hsize < target)
      hsize = hsize * 2
    end do
    allocate(h_key1(hsize), h_key2(hsize), h_val(hsize))
    h_key1 = 0
    h_key2 = 0
    h_val = 0
    do i = 1, nel
      in1 = elml(1,i)
      in2 = elml(2,i)

      jn1 = elml(2,i)
      jn2 = elml(3,i)

      kn1 = elml(3,i)
      kn2 = elml(1,i)

      if (elml(9,i) == 4) then
        kn1 = elml(3,i)
        kn2 = elml(4,i)

        nn1 = elml(4,i)
        nn2 = elml(1,i)
      end if

      call register_edge(5, in1, in2, i)
      call register_edge(6, jn1, jn2, i)
      call register_edge(7, kn1, kn2, i)
      if (elml(9,i) == 4) call register_edge(8, nn1, nn2, i)

      if (verbosity >= 2) then
        write(output_unit,*) 'Sorted element: ', i, nel
      end if
    end do

    do i = 1, nel
      inode1 = elml(1,i)
      inode2 = elml(2,i)
      inode3 = elml(3,i)

      if (elml(9,i) == 4) then
        inode4 = elml(4,i)
        rdum4x = nodel(1,inode4)
        rdum4y = nodel(2,inode4)
      else
        rdum4x = 0.0_real64
        rdum4y = 0.0_real64
      end if

      elnorm(1,i) = (nodel(1,inode1) + nodel(1,inode2) + &
                     nodel(1,inode3) + rdum4x) / real(elml(9,i), real64)
      elnorm(2,i) = (nodel(2,inode1) + nodel(2,inode2) + &
                     nodel(2,inode3) + rdum4y) / real(elml(9,i), real64)

      if (verbosity >= 2) then
        write(output_unit,*) 'Centroid element: ', i, nel
      end if
    end do

    do i = 1, nel
      in1 = elml(1,i)
      in2 = elml(2,i)

      jn1 = elml(2,i)
      jn2 = elml(3,i)

      kn1 = elml(3,i)
      kn2 = elml(1,i)

      if (elml(9,i) == 4) then
        kn1 = elml(3,i)
        kn2 = elml(4,i)

        nn1 = elml(4,i)
        nn2 = elml(1,i)
      end if

      rxc = (nodel(1,in1) - nodel(1,in2))
      ryc = (nodel(2,in1) - nodel(2,in2))
      rmx = ryc
      rmy = rxc
      rlen = sqrt(rxc**2 + ryc**2)

      xx = (nodel(1,in1) + nodel(1,in2)) / 2.0_real64
      yy = (nodel(2,in1) + nodel(2,in2)) / 2.0_real64

      rs = 1.0e-6_real64
      xx1 = xx + (-rmx / rlen) * rs
      yy1 = yy + (rmy / rlen) * rs

      xx2 = xx + (rmx / rlen) * rs
      yy2 = yy + (-rmy / rlen) * rs

      rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
      rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

      if (rd1 > rd2) then
        elnorm(4,i) = -rmx / rlen
        elnorm(5,i) = rmy / rlen
      else
        elnorm(4,i) = rmx / rlen
        elnorm(5,i) = -rmy / rlen
      end if

      rxc = (nodel(1,jn1) - nodel(1,jn2))
      ryc = (nodel(2,jn1) - nodel(2,jn2))
      rmx = ryc
      rmy = rxc
      rlen = sqrt(rxc**2 + ryc**2)

      xx = (nodel(1,jn1) + nodel(1,jn2)) / 2.0_real64
      yy = (nodel(2,jn1) + nodel(2,jn2)) / 2.0_real64

      rs = 1.0e-6_real64
      xx1 = xx + (-rmx / rlen) * rs
      yy1 = yy + (rmy / rlen) * rs

      xx2 = xx + (rmx / rlen) * rs
      yy2 = yy + (-rmy / rlen) * rs

      rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
      rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

      if (rd1 > rd2) then
        elnorm(6,i) = -rmx / rlen
        elnorm(7,i) = rmy / rlen
      else
        elnorm(6,i) = rmx / rlen
        elnorm(7,i) = -rmy / rlen
      end if

      rxc = (nodel(1,kn1) - nodel(1,kn2))
      ryc = (nodel(2,kn1) - nodel(2,kn2))
      rmx = ryc
      rmy = rxc
      rlen = sqrt(rxc**2 + ryc**2)

      xx = (nodel(1,kn1) + nodel(1,kn2)) / 2.0_real64
      yy = (nodel(2,kn1) + nodel(2,kn2)) / 2.0_real64

      rs = 1.0e-6_real64
      xx1 = xx + (-rmx / rlen) * rs
      yy1 = yy + (rmy / rlen) * rs

      xx2 = xx + (rmx / rlen) * rs
      yy2 = yy + (-rmy / rlen) * rs

      rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
      rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

      if (rd1 > rd2) then
        elnorm(8,i) = -rmx / rlen
        elnorm(9,i) = rmy / rlen
      else
        elnorm(8,i) = rmx / rlen
        elnorm(9,i) = -rmy / rlen
      end if

      if (elml(9,i) == 4) then
        rxc = (nodel(1,nn1) - nodel(1,nn2))
        ryc = (nodel(2,nn1) - nodel(2,nn2))
        rmx = ryc
        rmy = rxc
        rlen = sqrt(rxc**2 + ryc**2)

        xx = (nodel(1,nn1) + nodel(1,nn2)) / 2.0_real64
        yy = (nodel(2,nn1) + nodel(2,nn2)) / 2.0_real64

        rs = 1.0e-6_real64
        xx1 = xx + (-rmx / rlen) * rs
        yy1 = yy + (rmy / rlen) * rs

        xx2 = xx + (rmx / rlen) * rs
        yy2 = yy + (-rmy / rlen) * rs

        rd1 = sqrt((elnorm(1,i) - xx1)**2 + (elnorm(2,i) - yy1)**2)
        rd2 = sqrt((elnorm(1,i) - xx2)**2 + (elnorm(2,i) - yy2)**2)

        if (rd1 > rd2) then
          elnorm(10,i) = -rmx / rlen
          elnorm(11,i) = rmy / rlen
        else
          elnorm(10,i) = rmx / rlen
          elnorm(11,i) = -rmy / rlen
        end if
      end if

      if (verbosity >= 2) then
        write(output_unit,*) 'Normals element: ', i, nel
      end if
    end do

    do ie = 1, nel
      if (elml(9,ie) == 4) then
        rxd1 = nodel(1,elml(1,ie)) - nodel(1,elml(3,ie))
        rxd2 = nodel(1,elml(2,ie)) - nodel(1,elml(4,ie))
        ryd1 = nodel(2,elml(1,ie)) - nodel(2,elml(3,ie))
        ryd2 = nodel(2,elml(2,ie)) - nodel(2,elml(4,ie))
        rd1 = sqrt(rxd1**2 + ryd1**2)
        rd2 = sqrt(rxd2**2 + ryd2**2)
        rxd1 = rxd1 / rd1
        rxd2 = rxd2 / rd2
        ryd1 = ryd1 / rd1
        ryd2 = ryd2 / rd2
        rth = acos(rxd1 * rxd2 + ryd1 * ryd2)
        rarea = rd1 * rd2 * sin(rth) / 2.0_real64

        elnorm(3,ie) = rarea
      else
        rxa = nodel(1,elml(1,ie))
        rya = nodel(2,elml(1,ie))
        rxb = nodel(1,elml(2,ie))
        ryb = nodel(2,elml(2,ie))
        rxc_tri = nodel(1,elml(3,ie))
        ryc_tri = nodel(2,elml(3,ie))

        elnorm(3,ie) = 0.5_real64 * ((rxa - rxc_tri) * (ryb - rya) - &
                                    (rxa - rxb) * (ryc_tri - rya))
        elnorm(3,ie) = abs(elnorm(3,ie))
      end if

      if (verbosity >= 2) then
        write(output_unit,*) 'Area element: ', ie, nel
      end if
    end do

    n_bdry = 0
    n_bdry_unset = 0
    n_bdry_other = 0
    n_int_with_bc = 0
    bc1 = 0
    bc2 = 0
    bc3 = 0
    do i = 1, ned
      if (edgl(4,i) == 0) then
        n_bdry = n_bdry + 1
        if (edgl(5,i) == 1) then
          bc1 = bc1 + 1
        elseif (edgl(5,i) == 2) then
          bc2 = bc2 + 1
        elseif (edgl(5,i) == 3) then
          bc3 = bc3 + 1
        elseif (edgl(5,i) == 0) then
          n_bdry_unset = n_bdry_unset + 1
        else
          n_bdry_other = n_bdry_other + 1
        end if
      else
        if (edgl(5,i) /= 0) n_int_with_bc = n_int_with_bc + 1
      end if
    end do

    if (verbosity >= 1) then
      write(output_unit,*) 'Boundary edge counts (bc=1/2/3/0/other):', &
        bc1, bc2, bc3, n_bdry_unset, n_bdry_other
    end if
    if (n_bdry_unset > 0) then
      write(output_unit,*) 'Warning: boundary edges with bc=0 (unassigned):', &
        n_bdry_unset
    end if
    if (n_bdry_other > 0) then
      write(output_unit,*) 'Warning: boundary edges with unknown bc tag:', &
        n_bdry_other
    end if
    if (n_int_with_bc > 0) then
      write(output_unit,*) 'Warning: interior edges with bc tag:', n_int_with_bc
    end if

    deallocate(h_key1, h_key2, h_val)

  contains
    subroutine register_edge(neled, n1, n2, elem)
      integer, intent(in) :: neled, n1, n2, elem
      integer :: a, b, pos, probes
      integer(int64) :: h

      if (n1 <= n2) then
        a = n1
        b = n2
      else
        a = n2
        b = n1
      end if

      h = int(a, int64) * 73856093_int64 + int(b, int64) * 19349663_int64
      pos = int(mod(h, int(hsize, int64))) + 1
      probes = 0
      do
        if (h_val(pos) == 0) then
          call add_edge(neled, n1, n2, ned, elem)
          h_val(pos) = ned
          h_key1(pos) = a
          h_key2(pos) = b
          exit
        elseif (h_key1(pos) == a .and. h_key2(pos) == b) then
          call add_edge_element_only(neled, h_val(pos), elem)
          exit
        else
          pos = pos + 1
          if (pos > hsize) pos = 1
          probes = probes + 1
          if (probes >= hsize) then
            write(error_unit,*) 'Error: edge hash table full'
            error stop
          end if
        end if
      end do
    end subroutine register_edge
  end subroutine sort_grid

  subroutine check_if_added(in1, in2, ntmp, ifadded, j)
    ! Find an edge (in1,in2) in the existing list, independent of direction.
    integer, intent(in) :: in1, in2, ntmp
    logical, intent(out) :: ifadded
    integer, intent(out) :: j
    integer :: jj

    ifadded = .false.
    j = 0
    do jj = 1, ntmp
      if ((edgl(1,jj) == in1 .and. edgl(2,jj) == in2) .or. &
          (edgl(1,jj) == in2 .and. edgl(2,jj) == in1)) then
        ifadded = .true.
        j = jj
        exit
      end if
    end do
  end subroutine check_if_added

  subroutine add_edge_element_only(neled, je, i)
    ! Attach an existing edge index to an element and tag the second neighbor.
    integer, intent(in) :: neled, je, i

    elml(neled,i) = je
    edgl(4,je) = i
  end subroutine add_edge_element_only

  subroutine add_edge(neled, in1, in2, ntmp, i)
    ! Create a new edge and apply any boundary condition tags if present.
    integer, intent(in) :: neled, in1, in2, i
    integer, intent(inout) :: ntmp
    integer :: j

    ntmp = ntmp + 1
    elml(neled,i) = ntmp

    edgl(1,ntmp) = in1
    edgl(2,ntmp) = in2
    edgl(3,ntmp) = i
    edgl(4,ntmp) = 0

    edgl(5,ntmp) = 0
    edgl(6,ntmp) = 0
    do j = 1, nbc
      if (edgb(1,j) == in1 .and. edgb(2,j) == in2) then
        edgl(5,ntmp) = edgb(5,j)
        edgl(6,ntmp) = edgb(6,j)
        exit
      end if
      if (edgb(1,j) == in2 .and. edgb(2,j) == in1) then
        edgl(5,ntmp) = edgb(5,j)
        edgl(6,ntmp) = edgb(6,j)
        exit
      end if
    end do

  end subroutine add_edge

  subroutine write_new_grid()
    ! Write optional VTK mesh and the case.fvc output.
    if (write_vtk_mesh /= 0) then
      if (verbosity >= 1) then
        write(output_unit,*) 'Writing grid vtk file'
      end if
      call write_grid_vtk()
      if (verbosity >= 1) then
        write(output_unit,*) 'Finished writing grid vtk file'
      end if
    end if

    if (verbosity >= 1) then
      write(output_unit,*) 'Writing case.fvc'
    end if
    call write_case_fvc()
    if (verbosity >= 1) then
      write(output_unit,*) 'Finished writing case.fvc'
    end if
  end subroutine write_new_grid

  subroutine write_case_fvc()
    ! Emit the FVC2_CASE file in the solver's expected ASCII layout.
    integer :: unit
    integer :: i, j

    open(newunit=unit, file='case.fvc', status='replace', action='write')
    write(unit,*) 'FVC2_CASE 1'

    write(unit,*) '[phys_names]'
    write(unit,*) ntags
    do i = 1, ntags
      write(unit,*) ntags_dim(i), ntags_num(i), trim(ntags_name(i))
    end do
    write(unit,*) '[end_phys_names]'

    write(unit,*) '[mesh]'
    write(unit,*) nn, nel, ned

    write(unit,*) 'node_coords'
    do i = 1, nn
      write(unit,*) nodel(1,i), nodel(2,i)
    end do

    write(unit,*) 'element_geom'
    do i = 1, nel
      write(unit,*) (elnorm(j,i), j = 1, 11)
    end do

    write(unit,*) 'element_connectivity'
    do i = 1, nel
      write(unit,*) (elml(j,i), j = 1, 9)
    end do

    write(unit,*) 'edge_connectivity'
    do i = 1, ned
      write(unit,*) (edgl(j,i), j = 1, 6)
    end do

    write(unit,*) '[end_mesh]'
    close(unit)
  end subroutine write_case_fvc

  subroutine read_grid()
    ! Read Gmsh 4.x or legacy 2.x mesh formats into the internal arrays.
    character(len=64) :: dum_str, phys_name
    character(len=1024) :: line
    real(real64) :: dum_real, mesh_version
    real(real64) :: rx, ry, rz, rpar1, rpar2, rpar3
    integer :: dum_int, dum_tag
    integer :: ios
    integer :: npoints, ncurves, nsurfs, nvols
    integer :: nblocks, nblock, minTag, maxTag
    integer :: ent_dim, ent_tag, parametric
    integer :: elem_tag, elem_type, nnodes
    integer :: idx, node_tag
    integer, allocatable :: block_tags(:)
    integer, allocatable :: nodemap(:)
    integer, allocatable :: curve_phys(:)
    integer :: ndum1, nelt
    integer :: in1, in2, in3, in4
    integer :: i, j, ib
    integer :: unit
    integer :: istat

    open(newunit=unit, file=trim(gmsh_filename), form='formatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot open mesh file: ', trim(gmsh_filename)
      error stop
    end if

    read(unit,*) dum_str
    read(unit,*) mesh_version, dum_int, dum_int
    read(unit,*) dum_str
    read(unit,*) dum_str

    ntags = 0
    if (dum_str == '$PhysicalNames') then
      read(unit,*) ntags
      do i = 1, ntags
        read(unit,*) ntags_dim(i), ntags_num(i), phys_name
        ntags_name(i) = phys_name

        ntags_type(i) = 0
        if (phys_name(1:4) == 'wall' .or. &
            phys_name(1:4) == 'bump' .or. &
            phys_name(1:4) == 'top_' .or. &
            phys_name(1:4) == 'symm') then
          ntags_type(i) = 1
        elseif (phys_name(1:4) == 'infl' .or. &
                phys_name(1:4) == 'inle') then
          ntags_type(i) = 2
        elseif (phys_name(1:4) == 'outf' .or. &
                phys_name(1:4) == 'outl') then
          ntags_type(i) = 3
        end if
      end do
      read(unit,*) dum_str
      read(unit,*) dum_str
    end if

    if (mesh_version >= 4.0_real64) then
      allocate(block_tags(lm), nodemap(lm), curve_phys(lm), stat=istat)
      if (istat /= 0) then
        write(error_unit,*) 'Error: failed to allocate mesh buffers'
        error stop
      end if

      block_tags = 0
      nodemap = 0
      curve_phys = 0

      if (dum_str /= '$Entities') then
        write(error_unit,*) 'Error: expected $Entities, found: ', trim(dum_str)
        error stop
      end if
      read(unit,*) npoints, ncurves, nsurfs, nvols

      do i = 1, npoints
        read(unit,'(A)') line
        call parse_entity_line(line, 0, ent_tag, dum_int)
      end do
      do i = 1, ncurves
        read(unit,'(A)') line
        call parse_entity_line(line, 1, ent_tag, dum_int)
        if (ent_tag >= 1 .and. ent_tag <= lm) then
          curve_phys(ent_tag) = dum_int
        end if
      end do
      do i = 1, nsurfs
        read(unit,'(A)') line
        call parse_entity_line(line, 2, ent_tag, dum_int)
      end do
      do i = 1, nvols
        read(unit,'(A)') line
        call parse_entity_line(line, 3, ent_tag, dum_int)
      end do
      read(unit,*) dum_str

      read(unit,*) dum_str
      if (dum_str /= '$Nodes') then
        write(error_unit,*) 'Error: expected $Nodes, found: ', trim(dum_str)
        error stop
      end if
      read(unit,*) nblocks, nn, minTag, maxTag
      if (nn > lm) then
        write(error_unit,*) 'Error: too many nodes: ', nn, ' > ', lm
        error stop
      end if

      idx = 0
      do ib = 1, nblocks
        read(unit,*) ent_dim, ent_tag, parametric, nblock
        if (nblock > lm) then
          write(error_unit,*) 'Error: node block too large: ', nblock
          error stop
        end if
        call read_int_list(unit, block_tags, nblock)
        do i = 1, nblock
          if (parametric == 0) then
            read(unit,*) rx, ry, rz
          else
            if (ent_dim == 1) then
              read(unit,*) rx, ry, rz, rpar1
            elseif (ent_dim == 2) then
              read(unit,*) rx, ry, rz, rpar1, rpar2
            else
              read(unit,*) rx, ry, rz, rpar1, rpar2, rpar3
            end if
          end if
          idx = idx + 1
          node_tag = block_tags(i)
          if (node_tag > lm) then
            write(error_unit,*) 'Error: node tag out of range: ', node_tag
            error stop
          end if
          nodemap(node_tag) = idx
          nodel(1,idx) = rx * coord_scale
          nodel(2,idx) = ry * coord_scale
        end do
      end do
      read(unit,*) dum_str

      read(unit,*) dum_str
      if (dum_str /= '$Elements') then
        write(error_unit,*) 'Error: expected $Elements, found: ', trim(dum_str)
        error stop
      end if
      read(unit,*) nblocks, dum_int, minTag, maxTag
      nbc = 0
      nel = 0
      do ib = 1, nblocks
        read(unit,*) ent_dim, ent_tag, elem_type, nblock
        nnodes = gmsh_elem_nnodes(elem_type)
        if (nnodes <= 0) then
          write(error_unit,*) 'Error: unsupported element type: ', elem_type
          error stop
        end if

        if (nnodes == 1) then
          do i = 1, nblock
            read(unit,*) elem_tag, node_tag
          end do
        elseif (nnodes == 2) then
          do i = 1, nblock
            read(unit,*) elem_tag, in1, in2
            nbc = nbc + 1
            if (nodemap(in1) == 0 .or. nodemap(in2) == 0) then
              write(error_unit,*) 'Error: element node tag not found'
              error stop
            end if
            edgb(1,nbc) = nodemap(in1)
            edgb(2,nbc) = nodemap(in2)
            edgb(5,nbc) = curve_phys(ent_tag)
            edgb(6,nbc) = curve_phys(ent_tag)
            ndum1 = edgb(5,nbc)
            do j = 1, ntags
              if (ndum1 == ntags_num(j)) then
                edgb(5,nbc) = ntags_type(j)
              end if
            end do
          end do
        elseif (nnodes == 3) then
          do i = 1, nblock
            read(unit,*) elem_tag, in1, in2, in3
            nel = nel + 1
            elml(1,nel) = nodemap(in1)
            elml(2,nel) = nodemap(in2)
            elml(3,nel) = nodemap(in3)
            elml(9,nel) = 3
          end do
        elseif (nnodes == 4) then
          do i = 1, nblock
            read(unit,*) elem_tag, in1, in2, in3, in4
            nel = nel + 1
            elml(1,nel) = nodemap(in1)
            elml(2,nel) = nodemap(in2)
            elml(3,nel) = nodemap(in3)
            elml(4,nel) = nodemap(in4)
            elml(9,nel) = 4
          end do
        end if
      end do
      read(unit,*) dum_str
      close(unit)

      deallocate(block_tags, nodemap, curve_phys)

      nmm = neq * nel
      return
    end if

    if (dum_str /= '$Nodes') then
      write(error_unit,*) 'Error: expected $Nodes, found: ', trim(dum_str)
      error stop
    end if
    read(unit,*) nn
    do i = 1, nn
      read(unit,*) dum_int, rx, ry, dum_real
      nodel(1,i) = rx * coord_scale
      nodel(2,i) = ry * coord_scale
    end do
    read(unit,*) dum_str

    nbc = 0
    nel = 0
    if (dum_str /= '$EndNodes') then
      write(error_unit,*) 'Error: expected $EndNodes, found: ', trim(dum_str)
      error stop
    end if
    read(unit,*) dum_str
    if (dum_str /= '$Elements') then
      write(error_unit,*) 'Error: expected $Elements, found: ', trim(dum_str)
      error stop
    end if
    read(unit,*) nelt
    do i = 1, nelt
      read(unit,*) dum_int, dum_int, dum_tag

      if (dum_int == 1) then
        nbc = nbc + 1
        backspace unit
        read(unit,*) dum_int, dum_int, dum_int, edgb(5,nbc), dum_int, &
          edgb(1,nbc), edgb(2,nbc)
        edgb(6,nbc) = edgb(5,nbc)
        ndum1 = edgb(5,nbc)
        do j = 1, ntags
          if (ndum1 == ntags_num(j)) then
            edgb(5,nbc) = ntags_type(j)
          end if
        end do
      elseif (dum_int == 2) then
        nel = nel + 1
        backspace unit
        read(unit,*) dum_int, dum_int, dum_int, dum_int, dum_int, &
          elml(1,nel), elml(2,nel), elml(3,nel)
        elml(9,nel) = 3
      elseif (dum_int == 3) then
        nel = nel + 1
        backspace unit
        read(unit,*) dum_int, dum_int, dum_int, dum_int, dum_int, &
          elml(1,nel), elml(2,nel), elml(3,nel), elml(4,nel)
        elml(9,nel) = 4
      end if
    end do
    read(unit,*) dum_str
    close(unit)

    nmm = neq * nel
  end subroutine read_grid

  integer function gmsh_elem_nnodes(elem_type)
    ! Map Gmsh element type IDs to node counts.
    integer, intent(in) :: elem_type

    select case (elem_type)
    case (1)
      gmsh_elem_nnodes = 2
    case (2)
      gmsh_elem_nnodes = 3
    case (3)
      gmsh_elem_nnodes = 4
    case (15)
      gmsh_elem_nnodes = 1
    case default
      gmsh_elem_nnodes = 0
    end select
  end function gmsh_elem_nnodes

  subroutine parse_entity_line(line, dim, ent_tag, phys_tag)
    ! Parse $Entities lines to capture physical tags (Gmsh 4.x).
    character(len=*), intent(in) :: line
    integer, intent(in) :: dim
    integer, intent(out) :: ent_tag, phys_tag
    integer :: num_phys, num_bnd
    integer :: ntok, p
    character(len=32) :: tokens(512)

    call tokenize(line, tokens, ntok)
    p = 1
    if (ntok <= 0) then
      write(error_unit,*) 'Error: failed to parse $Entities line:'
      write(error_unit,*) trim(line)
      error stop
    end if

    read(tokens(p),*) ent_tag
    p = p + 1

    if (dim == 0) then
      p = p + 3
    else
      p = p + 6
    end if
    if (p > ntok) then
      write(error_unit,*) 'Error: failed to parse $Entities line:'
      write(error_unit,*) trim(line)
      error stop
    end if

    read(tokens(p),*) num_phys
    p = p + 1

    phys_tag = 0
    if (num_phys > 0) then
      read(tokens(p),*) phys_tag
    end if
    p = p + num_phys
    if (p > ntok) return

    read(tokens(p),*) num_bnd
    p = p + 1 + num_bnd
  end subroutine parse_entity_line

  subroutine tokenize(line, tokens, ntok)
    ! Split a line into whitespace-delimited tokens (fixed buffer size).
    character(len=*), intent(in) :: line
    character(len=32), intent(out) :: tokens(*)
    integer, intent(out) :: ntok
    integer :: i, linelen, start, tlen
    integer, parameter :: max_tok = 512

    ntok = 0
    linelen = len_trim(line)

    i = 1
    do while (i <= linelen)
      do while (i <= linelen .and. &
                (line(i:i) == ' ' .or. line(i:i) == char(9)))
        i = i + 1
      end do
      if (i > linelen) exit

      start = i
      do while (i <= linelen .and. &
                (line(i:i) /= ' ' .and. line(i:i) /= char(9)))
        i = i + 1
      end do

      tlen = i - start
      ntok = ntok + 1
      if (ntok > max_tok) then
        write(error_unit,*) 'Error: too many tokens in entity line'
        error stop
      end if
      tokens(ntok) = ' '
      if (tlen > len(tokens(ntok))) tlen = len(tokens(ntok))
      tokens(ntok)(1:tlen) = line(start:start+tlen-1)
    end do
  end subroutine tokenize

  subroutine read_int_list(unit, arr, n)
    ! Read an integer list that may span multiple lines.
    integer, intent(in) :: unit, n
    integer, intent(out) :: arr(n)
    integer :: idx, i, ntok
    character(len=1024) :: line
    character(len=32) :: tokens(512)

    idx = 0
    do while (idx < n)
      read(unit,'(A)') line
      call tokenize(line, tokens, ntok)
      do i = 1, ntok
        idx = idx + 1
        read(tokens(i),*) arr(idx)
        if (idx >= n) exit
      end do
    end do
  end subroutine read_int_list

  subroutine write_grid_vtk()
    ! Write a debug VTK mesh file (optional).
    character(len=15) :: fnamec
    integer, save :: icalld = -1
    integer :: unit
    integer :: i

    icalld = icalld + 1

    write(fnamec,'(A6,I5.5,A4)') 'grid_c', icalld, '.vtk'
    open(newunit=unit, file=fnamec, action='write', status='replace')

    write(unit,'(A26)') '# vtk DataFile Version 3.0'
    write(unit,*) '2D scalar data'
    write(unit,*) 'ASCII'
    write(unit,*) ''

    write(unit,*) 'DATASET UNSTRUCTURED_GRID'
    write(unit,*) 'POINTS', ' ', nn, ' ', ' float'
    do i = 1, nn
      write(unit,*) nodel(1,i), nodel(2,i), 0.0_real64
    end do
    write(unit,*) ''

    write(unit,*) 'CELLS ', nel, nel * 5
    do i = 1, nel
      if (elml(9,i) == 3) then
        write(unit,*) elml(9,i), elml(1,i)-1, elml(2,i)-1, elml(3,i)-1
      elseif (elml(9,i) == 4) then
        write(unit,*) elml(9,i), elml(1,i)-1, elml(2,i)-1, elml(3,i)-1, &
          elml(4,i)-1
      end if
    end do
    write(unit,*) ''

    write(unit,*) 'CELL_TYPES', nel
    do i = 1, nel
      if (elml(9,i) == 3) write(unit,*) 5
      if (elml(9,i) == 4) write(unit,*) 9
    end do
    write(unit,*) ''

    write(unit,*) 'CELL_DATA', nel
    write(unit,*) 'SCALARS  dumdata  float 1'
    write(unit,*) 'LOOKUP_TABLE default'
    do i = 1, nel
      write(unit,*) 0.0_real64
    end do
    write(unit,*) ''

    close(unit)
  end subroutine write_grid_vtk

  subroutine read_user_input_file()
    ! Legacy input path (not used by CLI mode).
    integer :: unit
    integer :: ios

    open(newunit=unit, file='prj.ini', form='formatted', action='read')
    read(unit,*) gmsh_filename
    write_vtk_mesh = 1
    verbosity = 1
    read(unit,*,iostat=ios) write_vtk_mesh
    if (ios == 0) then
      read(unit,*,iostat=ios) verbosity
      if (ios /= 0) verbosity = 1
    end if
    close(unit)
  end subroutine read_user_input_file
end module msh2fvc_ops

program msh2fvc_setup
  use msh2fvc_data
  use msh2fvc_ops
  implicit none

  call reset_arrays()

  call read_cmdline()

  call read_grid()
  call sort_grid()
  call write_new_grid()
end program msh2fvc_setup
