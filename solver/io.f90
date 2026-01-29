!------------------------------------------------------------------------------
! Case I/O (case.fvc / case.vars) and output writing.
!------------------------------------------------------------------------------
module io
  use iso_fortran_env, only: real64, output_unit, error_unit
  use data
  use utils
  implicit none
contains
  subroutine write_soln()
    character(len=32) :: fnamec
    character(len=256) :: outdir, outfile
    integer, save :: icalld = -1
    integer :: cells_size
    integer :: nbad_cfl
    real(real64) :: cfl_max, cfl_local, rho, uvel, vvel, p, a, speed, area, dx
    real(real64) :: ru, rv, re, etotal
    integer :: i
    integer :: unit

    icalld = icalld + 1
    cfl_max = 0.0_real64
    nbad_cfl = 0
    do i = 1, num_elements
      rho = state(i,1)
      if (rho <= 0.0_real64) then
        nbad_cfl = nbad_cfl + 1
        cycle
      end if
      uvel = state(i,2) / rho
      vvel = state(i,3) / rho
      p = (gamma_gas - 1.0_real64) * &
          (state(i,4) - 0.5_real64*(state(i,2)*state(i,2) + state(i,3)*state(i,3)) / rho)
      if (p <= 0.0_real64) then
        nbad_cfl = nbad_cfl + 1
        cycle
      end if
      area = element_geom(3,i)
      if (area <= 0.0_real64) then
        nbad_cfl = nbad_cfl + 1
        cycle
      end if
      dx = sqrt(area)
      a = sqrt(gamma_gas*p/rho)
      speed = sqrt(uvel*uvel + vvel*vvel)
      cfl_local = time_step * (speed + a) / dx
      if (cfl_local > cfl_max) cfl_max = cfl_local
    end do

    write (output_unit,'(A,I6,2X,ES12.6,2X,ES12.6,2X,A,ES12.6)') &
      'output_interval: ', icalld, simulation_time, time_step, 'CFLmax=', cfl_max
    if (nbad_cfl > 0) then
      write (output_unit,'(A,I6)') '  CFL skipped cells: ', nbad_cfl
    end if

    call ensure_output_dir()
    outdir = 'output/'

    cells_size = 0
    do i = 1, num_elements
      cells_size = cells_size + element_connectivity(9,i) + 1
    end do

    if (output_format_mode == 1 .or. output_format_mode == 2) then
      if (animate_output) then
        call write_soln_vtu(icalld, outdir)
        call write_pvd_series(icalld, outdir, '.vtu')
      else
        write(fnamec,'(A6,I5.5,A4)') 'soln_c', icalld, '.vtk'
        outfile = trim(outdir) // trim(fnamec)
        open(newunit=unit, file=outfile, action='write', status='replace')

        write(unit,'(A26)') '# vtk DataFile Version 3.0'
        write(unit,*) '2D scalar data'
        write(unit,*) 'ASCII'
        write(unit,*) ''

        write(unit,*) 'DATASET UNSTRUCTURED_GRID'
        write(unit,*) 'POINTS', ' ', num_nodes, ' ', ' float'
        do i = 1, num_nodes
          write(unit,*) node_coords(1,i), node_coords(2,i), 0.0_real64
        end do
        write(unit,*) ''

        write(unit,*) 'CELLS ', num_elements, cells_size
        do i = 1, num_elements
          if (element_connectivity(9,i) == 3) then
            write(unit,*) 3, element_connectivity(1,i)-1, element_connectivity(2,i)-1, element_connectivity(3,i)-1
          elseif (element_connectivity(9,i) == 4) then
            write(unit,*) 4, element_connectivity(1,i)-1, element_connectivity(2,i)-1, element_connectivity(3,i)-1, element_connectivity(4,i)-1
          end if
        end do
        write(unit,*) ''

        write(unit,*) 'CELL_TYPES', num_elements
        do i = 1, num_elements
          if (element_connectivity(9,i) == 3) write(unit,*) 5
          if (element_connectivity(9,i) == 4) write(unit,*) 9
        end do
        write(unit,*) ''

        write(unit,*) 'CELL_DATA', num_elements
        write(unit,*) 'SCALARS density float 1'
        write(unit,*) 'LOOKUP_TABLE default'
        do i = 1, num_elements
          write(unit,*) state(i,1)
        end do
        write(unit,*) 'SCALARS ru float 1'
        write(unit,*) 'LOOKUP_TABLE default'
        do i = 1, num_elements
          write(unit,*) state(i,2)
        end do
        write(unit,*) 'SCALARS rv float 1'
        write(unit,*) 'LOOKUP_TABLE default'
        do i = 1, num_elements
          write(unit,*) state(i,3)
        end do
        write(unit,*) 'SCALARS re float 1'
        write(unit,*) 'LOOKUP_TABLE default'
        do i = 1, num_elements
          write(unit,*) state(i,4)
        end do
        write(unit,*) 'VECTORS velocity float'
        do i = 1, num_elements
          rho = state(i,1)
          if (rho > 0.0_real64) then
            uvel = state(i,2) / rho
            vvel = state(i,3) / rho
          else
            uvel = 0.0_real64
            vvel = 0.0_real64
          end if
          write(unit,*) uvel, vvel, 0.0_real64
        end do
        write(unit,*) 'SCALARS pressure float 1'
        write(unit,*) 'LOOKUP_TABLE default'
        do i = 1, num_elements
          rho = state(i,1)
          ru = state(i,2)
          rv = state(i,3)
          re = state(i,4)
          if (rho > 0.0_real64) then
            p = (gamma_gas - 1.0_real64) * (re - 0.5_real64*(ru*ru + rv*rv) / rho)
          else
            p = 0.0_real64
          end if
          write(unit,*) p
        end do
        write(unit,*) 'SCALARS total_energy float 1'
        write(unit,*) 'LOOKUP_TABLE default'
        do i = 1, num_elements
          rho = state(i,1)
          re = state(i,4)
          if (rho > 0.0_real64) then
            etotal = re / rho
          else
            etotal = 0.0_real64
          end if
          write(unit,*) etotal
        end do
        write(unit,*) ''

        close(unit)
      end if
    end if

    if (output_format_mode == 0 .or. output_format_mode == 2) then
      write(fnamec,'(A6,I5.5,A4)') 'soln_c', icalld, '.dat'
      outfile = trim(outdir) // trim(fnamec)
      open(newunit=unit, file=outfile, action='write', status='replace')

      write(unit,*) 'TITLE = "Your Title Here"'
      write(unit,*) 'VARIABLES = "X", "Y", "Z", "Density", "RU", "RV", "RE"'
      write(unit,*) 'ZONE N=', num_nodes, ', E=', num_elements, ', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'

      do i = 1, num_nodes
        write(unit,*) node_coords(1,i), node_coords(2,i), 0.0_real64, state(i,1), state(i,2), state(i,3), state(i,4)
      end do

      do i = 1, num_elements
        if (element_connectivity(9,i) == 3) then
          write(unit,*) element_connectivity(1,i), element_connectivity(2,i), element_connectivity(3,i)
        elseif (element_connectivity(9,i) == 4) then
          write(unit,*) element_connectivity(1,i), element_connectivity(2,i), element_connectivity(3,i), element_connectivity(4,i)
        end if
      end do

      close(unit)
    end if
  end subroutine write_soln

  subroutine write_pvd_series(output_index, outdir, file_ext)
    integer, intent(in) :: output_index
    character(len=*), intent(in) :: outdir
    character(len=*), intent(in) :: file_ext
    integer :: unit
    integer :: idx
    real(real64) :: tval
    character(len=32) :: fnamec
    character(len=256) :: outfile

    outfile = trim(outdir) // 'soln.pvd'
    open(newunit=unit, file=outfile, action='write', status='replace')

    write(unit,'(A)') '<?xml version="1.0"?>'
    write(unit,'(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
    write(unit,'(A)') '  <Collection>'

    do idx = 0, output_index
      tval = real(idx * output_interval, real64) * time_step
      write(fnamec,'(A6,I5.5,A)') 'soln_c', idx, trim(file_ext)
      write(unit,'(A,F12.6,A,A,A)') '    <DataSet timestep="', tval, &
        '" group="" part="0" file="', trim(fnamec), '"/>'
    end do

    write(unit,'(A)') '  </Collection>'
    write(unit,'(A)') '</VTKFile>'
    close(unit)
  end subroutine write_pvd_series

  subroutine write_soln_vtu(output_index, outdir)
    integer, intent(in) :: output_index
    character(len=*), intent(in) :: outdir
    integer :: unit
    integer :: i
    integer :: offset
    real(real64) :: rho, uvel, vvel, p, ru, rv, re, etotal
    character(len=32) :: fnamec
    character(len=256) :: outfile

    write(fnamec,'(A6,I5.5,A4)') 'soln_c', output_index, '.vtu'
    outfile = trim(outdir) // trim(fnamec)
    open(newunit=unit, file=outfile, action='write', status='replace')

    write(unit,'(A)') '<?xml version="1.0"?>'
    write(unit,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(unit,'(A,I0,A,I0,A)') '  <UnstructuredGrid><Piece NumberOfPoints="', num_nodes, &
      '" NumberOfCells="', num_elements, '">'

    write(unit,'(A)') '    <Points>'
    write(unit,'(A)') '      <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    do i = 1, num_nodes
      write(unit,*) node_coords(1,i), node_coords(2,i), 0.0_real64
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '    </Points>'

    write(unit,'(A)') '    <Cells>'
    write(unit,'(A)') '      <DataArray type="Int32" Name="connectivity" format="ascii">'
    do i = 1, num_elements
      if (element_connectivity(9,i) == 3) then
        write(unit,*) element_connectivity(1,i)-1, element_connectivity(2,i)-1, element_connectivity(3,i)-1
      elseif (element_connectivity(9,i) == 4) then
        write(unit,*) element_connectivity(1,i)-1, element_connectivity(2,i)-1, &
          element_connectivity(3,i)-1, element_connectivity(4,i)-1
      end if
    end do
    write(unit,'(A)') '      </DataArray>'

    write(unit,'(A)') '      <DataArray type="Int32" Name="offsets" format="ascii">'
    offset = 0
    do i = 1, num_elements
      offset = offset + element_connectivity(9,i)
      write(unit,*) offset
    end do
    write(unit,'(A)') '      </DataArray>'

    write(unit,'(A)') '      <DataArray type="UInt8" Name="types" format="ascii">'
    do i = 1, num_elements
      if (element_connectivity(9,i) == 3) write(unit,*) 5
      if (element_connectivity(9,i) == 4) write(unit,*) 9
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '    </Cells>'

    write(unit,'(A)') '    <CellData Scalars="density">'
    write(unit,'(A)') '      <DataArray type="Float64" Name="density" format="ascii">'
    do i = 1, num_elements
      write(unit,*) state(i,1)
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '      <DataArray type="Float64" Name="ru" format="ascii">'
    do i = 1, num_elements
      write(unit,*) state(i,2)
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '      <DataArray type="Float64" Name="rv" format="ascii">'
    do i = 1, num_elements
      write(unit,*) state(i,3)
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '      <DataArray type="Float64" Name="re" format="ascii">'
    do i = 1, num_elements
      write(unit,*) state(i,4)
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '      <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii">'
    do i = 1, num_elements
      rho = state(i,1)
      if (rho > 0.0_real64) then
        uvel = state(i,2) / rho
        vvel = state(i,3) / rho
      else
        uvel = 0.0_real64
        vvel = 0.0_real64
      end if
      write(unit,*) uvel, vvel, 0.0_real64
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '      <DataArray type="Float64" Name="pressure" format="ascii">'
    do i = 1, num_elements
      rho = state(i,1)
      ru = state(i,2)
      rv = state(i,3)
      re = state(i,4)
      if (rho > 0.0_real64) then
        p = (gamma_gas - 1.0_real64) * (re - 0.5_real64*(ru*ru + rv*rv) / rho)
      else
        p = 0.0_real64
      end if
      write(unit,*) p
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '      <DataArray type="Float64" Name="total_energy" format="ascii">'
    do i = 1, num_elements
      rho = state(i,1)
      re = state(i,4)
      if (rho > 0.0_real64) then
        etotal = re / rho
      else
        etotal = 0.0_real64
      end if
      write(unit,*) etotal
    end do
    write(unit,'(A)') '      </DataArray>'
    write(unit,'(A)') '    </CellData>'

    write(unit,'(A)') '  </Piece></UnstructuredGrid>'
    write(unit,'(A)') '</VTKFile>'
    close(unit)
  end subroutine write_soln_vtu

  subroutine read_sorted_grid()
    character(len=256) :: mesh_dir, path
    logical :: use_grid
    integer :: ios
    integer :: i
    integer :: unit

    mesh_dir = ''
    inquire(file='grid/node_coords', exist=use_grid)
    if (use_grid) then
      mesh_dir = 'grid/'
    end if

    write(output_unit,*) 'Reading node_coords file'
    path = trim(mesh_dir) // 'node_coords'
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot open mesh file: ', trim(path)
      error stop
    end if
    read(unit,*) num_nodes
    do i = 1, num_nodes
      read(unit,*) node_coords(1,i)
      read(unit,*) node_coords(2,i)
    end do
    close(unit)
    write(output_unit,*) 'Finished reading node_coords file'

    write(output_unit,*) 'Reading element_geom file'
    path = trim(mesh_dir) // 'element_geom'
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot open mesh file: ', trim(path)
      error stop
    end if
    read(unit,*) num_elements
    do i = 1, num_elements
      read(unit,*) element_geom(1,i)
      read(unit,*) element_geom(2,i)
      read(unit,*) element_geom(3,i)
      read(unit,*) element_geom(4,i)
      read(unit,*) element_geom(5,i)
      read(unit,*) element_geom(6,i)
      read(unit,*) element_geom(7,i)
      read(unit,*) element_geom(8,i)
      read(unit,*) element_geom(9,i)
      read(unit,*) element_geom(10,i)
      read(unit,*) element_geom(11,i)
    end do
    close(unit)
    write(output_unit,*) 'Finished reading element_geom file'

    write(output_unit,*) 'Reading element_connectivity file'
    path = trim(mesh_dir) // 'element_connectivity'
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot open mesh file: ', trim(path)
      error stop
    end if
    read(unit,*) num_elements
    do i = 1, num_elements
      read(unit,*) element_connectivity(1,i)
      read(unit,*) element_connectivity(2,i)
      read(unit,*) element_connectivity(3,i)
      read(unit,*) element_connectivity(4,i)
      read(unit,*) element_connectivity(5,i)
      read(unit,*) element_connectivity(6,i)
      read(unit,*) element_connectivity(7,i)
      read(unit,*) element_connectivity(8,i)
      read(unit,*) element_connectivity(9,i)
    end do
    close(unit)
    write(output_unit,*) 'Finished reading element_connectivity file'

    write(output_unit,*) 'Reading edge_connectivity file'
    path = trim(mesh_dir) // 'edge_connectivity'
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot open mesh file: ', trim(path)
      error stop
    end if
    read(unit,*) num_edges
    do i = 1, num_edges
      read(unit,*) edge_connectivity(1,i)
      read(unit,*) edge_connectivity(2,i)
      read(unit,*) edge_connectivity(3,i)
      read(unit,*) edge_connectivity(4,i)
      read(unit,*) edge_connectivity(5,i)
      read(unit,*) edge_connectivity(6,i)
    end do
    close(unit)
    write(output_unit,*) 'Finished reading edge_connectivity file'

    call map_boundary_edges()

    num_state_entries = num_equations * num_elements
  end subroutine read_sorted_grid

  subroutine read_case_file()
    character(len=256) :: line
    logical :: have_phys, have_mesh
    integer :: ios
    integer :: unit

    have_phys = .false.
    have_mesh = .false.

    open(newunit=unit, file='case.fvc', form='formatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot open case.fvc'
      error stop
    end if
    call read_next_line(unit, line, ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot read case.fvc'
      error stop
    end if
    call lower_string(line)
    if (index(line, 'fvc2_case') /= 1) then
      write(error_unit,*) 'Error: invalid case.fvc header'
      error stop
    end if

    do
      call read_next_line(unit, line, ios)
      if (ios /= 0) exit
      call lower_string(line)
      if (line == '[phys_names]') then
        call read_phys_names_section(unit)
        have_phys = .true.
      elseif (line == '[vars]') then
        write(error_unit,*) 'Error: case.fvc should not contain vars'
        write(error_unit,*) '       Use a separate case.vars file.'
        error stop
      elseif (line == '[mesh]') then
        call read_mesh_section(unit)
        have_mesh = .true.
      end if
    end do
    close(unit)

    if (.not. have_phys) then
      write(error_unit,*) 'Error: missing phys_names section in case.fvc'
      error stop
    end if
    if (.not. have_mesh) then
      write(error_unit,*) 'Error: missing mesh section in case.fvc'
      error stop
    end if

    call read_vars_file()
    call map_boundary_edges()
  end subroutine read_case_file

  subroutine read_phys_names_section(unit)
    integer, intent(in) :: unit
    integer :: ios, n_all, dim, tag
    character(len=64) :: pname
    character(len=256) :: line
    integer :: i

    num_bc_defs = 0
    read(unit,*,iostat=ios) n_all
    if (ios /= 0) then
      write(error_unit,*) 'Error: malformed phys_names header in case.fvc'
      error stop
    end if

    do i = 1, n_all
      read(unit,*,iostat=ios) dim, tag, pname
      if (ios /= 0) then
        write(error_unit,*) 'Error: malformed phys_names entry in case.fvc'
        error stop
      end if
      call lower_string(pname)
      if (dim == 1) then
        num_bc_defs = num_bc_defs + 1
        if (num_bc_defs > max_boundary_conditions) then
          write(error_unit,*) 'Error: too many BC names in case.fvc'
          error stop
        end if
        bc_label(num_bc_defs) = trim(pname)
        bc_physical_tag(num_bc_defs) = tag
      end if
    end do

    call read_next_line(unit, line, ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: missing end_phys_names in case.fvc'
      error stop
    end if
    call lower_string(line)
    if (line /= '[end_phys_names]') then
      write(error_unit,*) 'Error: expected [end_phys_names] in case.fvc'
      error stop
    end if

    if (num_bc_defs == 0) then
      write(error_unit,*) 'Error: phys_names has no 1D entries'
      error stop
    end if
  end subroutine read_phys_names_section

  subroutine read_vars_file()
    real(real64) :: mmass, rideal
    integer :: unit, ios
    character(len=256) :: line
    character(len=64) :: key
    character(len=192) :: vals
    character(len=32) :: strval
    logical :: have_dt, have_iostep, have_nstepmax, have_gamma
    logical :: have_mmass, have_k, have_time_scheme, have_conv_scheme
    logical :: have_ic_type, have_ic_uniform, have_ic_split_x
    logical :: have_ic_left, have_ic_right, have_output_mode
    logical :: have_bc_block, in_bc_block
    logical :: have_animate_output
    integer :: ic_type_tmp
    integer :: i

    do i = 1, max_boundary_conditions
      bc_values(:,i) = 0.0_real64
      bc_kind(i) = 0
    end do
    have_dt = .false.
    have_iostep = .false.
    have_nstepmax = .false.
    have_gamma = .false.
    have_mmass = .false.
    have_k = .false.
    have_time_scheme = .false.
    have_conv_scheme = .false.
    have_ic_type = .false.
    have_ic_uniform = .false.
    have_ic_split_x = .false.
    have_ic_left = .false.
    have_ic_right = .false.
    have_output_mode = .false.
    have_animate_output = .false.
    have_bc_block = .false.
    in_bc_block = .false.

    initial_condition_type = 1
    ic_uniform_state = (/1.4_real64, 1.0_real64, 3.0_real64, 0.0_real64/)

    shock_split_x = 0.5_real64
    shock_left_state = (/1.0_real64, 1.0_real64, 0.0_real64, 0.0_real64/)
    shock_right_state = (/0.125_real64, 0.1_real64, 0.0_real64, 0.0_real64/)

    output_format_mode = 2
    animate_output = .false.

    open(newunit=unit, file='case.vars', form='formatted', status='old', iostat=ios)
    if (ios /= 0) then
      write(error_unit,*) 'Error: cannot open case.vars'
      error stop
    end if

    do
      call read_next_line(unit, line, ios)
      if (ios /= 0) exit
      if (in_bc_block) then
        if (trim(line) == '>') then
          in_bc_block = .false.
          have_bc_block = .true.
          cycle
        end if
        call parse_bc_entry(line, ios)
        if (ios /= 0) goto 982
        cycle
      end if

      call parse_kv_line(line, key, vals, ios)
      if (ios /= 0) goto 980

      if (key == 'time_step') then
        read(vals,*) time_step
        have_dt = .true.
      elseif (key == 'output_interval') then
        read(vals,*) output_interval
        have_iostep = .true.
      elseif (key == 'max_steps') then
        read(vals,*) max_steps
        have_nstepmax = .true.
      elseif (key == 'gamma') then
        read(vals,*) gamma_gas
        have_gamma = .true.
      elseif (key == 'molar_mass') then
        read(vals,*) mmass
        have_mmass = .true.
      elseif (key == 'diffusion_coeff') then
        read(vals,*) diffusion_coeff
        have_k = .true.
      elseif (key == 'time_discretization') then
        read(vals,*) time_integration_scheme
        have_time_scheme = .true.
      elseif (key == 'convective_scheme') then
        read(vals,*) convection_scheme
        have_conv_scheme = .true.
      elseif (key == 'boundary_conditions') then
        if (trim(vals) /= '<') goto 981
        in_bc_block = .true.
      elseif (key == 'initial_condition_type') then
        strval = ''
        read(vals,*,iostat=ios) strval
        if (ios /= 0) goto 950
        call lower_string(strval)
        if (strval == 'uniform') then
          initial_condition_type = 1
        elseif (strval == 'fwd_step' .or. strval == 'forward_step') then
          initial_condition_type = 2
        elseif (strval == 'shock_tube' .or. strval == 'shocktube') then
          initial_condition_type = 3
        else
          read(strval,*,iostat=ios) ic_type_tmp
          if (ios /= 0) goto 955
          initial_condition_type = ic_type_tmp
        end if
        have_ic_type = .true.
      elseif (key == 'ic_uniform_state') then
        read(vals,*) ic_uniform_state(1), ic_uniform_state(2), ic_uniform_state(3), ic_uniform_state(4)
        have_ic_uniform = .true.
      elseif (key == 'shock_split_x') then
        read(vals,*) shock_split_x
        have_ic_split_x = .true.
      elseif (key == 'shock_left_state') then
        read(vals,*) shock_left_state(1), shock_left_state(2), shock_left_state(3), shock_left_state(4)
        have_ic_left = .true.
      elseif (key == 'shock_right_state') then
        read(vals,*) shock_right_state(1), shock_right_state(2), shock_right_state(3), shock_right_state(4)
        have_ic_right = .true.
      elseif (key == 'output_format_mode') then
        strval = ''
        read(vals,*,iostat=ios) strval
        if (ios /= 0) goto 970
        call lower_string(strval)
        if (strval == 'tecplot') then
          output_format_mode = 0
        elseif (strval == 'vtk') then
          output_format_mode = 1
        elseif (strval == 'both') then
          output_format_mode = 2
        else
          read(strval,*,iostat=ios) output_format_mode
          if (ios /= 0) goto 970
        end if
        have_output_mode = .true.
      elseif (key == 'animate_output') then
        strval = ''
        read(vals,*,iostat=ios) strval
        if (ios /= 0) goto 970
        call lower_string(strval)
        if (strval == 'on' .or. strval == 'true' .or. strval == 'yes' .or. strval == '1') then
          animate_output = .true.
        elseif (strval == 'off' .or. strval == 'false' .or. strval == 'no' .or. strval == '0') then
          animate_output = .false.
        else
          goto 970
        end if
        have_animate_output = .true.
      else
        goto 9400
      end if
    end do

    if (in_bc_block) goto 981
    if (.not. have_dt) goto 910
    if (.not. have_iostep) goto 911
    if (.not. have_nstepmax) goto 912
    if (.not. have_gamma) goto 913
    if (.not. have_mmass) goto 914
    if (.not. have_k) goto 915
    if (.not. have_time_scheme) goto 916
    if (.not. have_conv_scheme) goto 917
    if (.not. have_bc_block) goto 918
    do i = 1, num_bc_defs
      if (bc_kind(i) == 0) goto 919
    end do

    if (have_ic_type) then
      if (initial_condition_type == 1 .and. .not. have_ic_uniform) goto 920
      if (initial_condition_type == 3) then
        if (.not. have_ic_split_x) goto 930
        if (.not. have_ic_left) goto 940
        if (.not. have_ic_right) goto 950
      end if
    end if

    rideal = 8.3144598_real64
    gas_constant  = rideal / mmass
    cp_gas    = gas_constant / (1.0_real64 - 1.0_real64/gamma_gas)
    cv_gas    = cp_gas / gamma_gas

    close(unit)
    return

910 continue
    write(error_unit,*) 'Error: missing key time_step in case.vars'
    error stop
911 continue
    write(error_unit,*) 'Error: missing key output_interval in case.vars'
    error stop
912 continue
    write(error_unit,*) 'Error: missing key max_steps in case.vars'
    error stop
913 continue
    write(error_unit,*) 'Error: missing key gamma in case.vars'
    error stop
914 continue
    write(error_unit,*) 'Error: missing key molar_mass in case.vars'
    error stop
915 continue
    write(error_unit,*) 'Error: missing key diffusion_coeff in case.vars'
    error stop
916 continue
    write(error_unit,*) 'Error: missing key time_discretization in case.vars'
    error stop
917 continue
    write(error_unit,*) 'Error: missing key convective_scheme in case.vars'
    error stop
918 continue
    write(error_unit,*) 'Error: missing boundary_conditions block in case.vars'
    error stop
919 continue
    write(error_unit,*) 'Error: boundary_conditions missing entry for: ', bc_label(i)
    error stop
920 continue
    write(error_unit,*) 'Error: initial_condition_type=uniform requires ic_uniform_state'
    error stop
930 continue
    write(error_unit,*) 'Error: initial_condition_type=shock_tube requires shock_split_x'
    error stop
940 continue
    write(error_unit,*) 'Error: initial_condition_type=shock_tube requires shock_left_state'
    error stop
950 continue
    write(error_unit,*) 'Error: initial_condition_type=shock_tube requires shock_right_state'
    error stop
955 continue
    write(error_unit,*) 'Error: invalid initial_condition_type in case.vars'
    error stop
970 continue
    write(error_unit,*) 'Error: invalid output_format_mode in case.vars'
    error stop
980 continue
    write(error_unit,*) 'Error: expected key: value format in case.vars'
    error stop
981 continue
    write(error_unit,*) 'Error: boundary_conditions must start with "<"'
    error stop
982 continue
    write(error_unit,*) 'Error: invalid boundary_conditions entry: ', line
    error stop
9400 continue
    write(error_unit,*) 'Error: unknown key in case.vars: ', key
    error stop
  end subroutine read_vars_file

  subroutine read_mesh_section(unit)
    integer, intent(in) :: unit
    integer :: ios
    character(len=256) :: line
    integer :: i, j

    read(unit,*,iostat=ios) num_nodes, num_elements, num_edges
    if (ios /= 0) then
      write(error_unit,*) 'Error: malformed mesh header in case.fvc'
      error stop
    end if

    call read_next_line(unit, line, ios)
    if (ios /= 0) goto 700
    call lower_string(line)
    if (line /= 'node_coords') goto 700
    do i = 1, num_nodes
      read(unit,*,iostat=ios) node_coords(1,i), node_coords(2,i)
      if (ios /= 0) goto 701
    end do

    call read_next_line(unit, line, ios)
    if (ios /= 0) goto 700
    call lower_string(line)
    if (line /= 'element_geom') goto 700
    do i = 1, num_elements
      read(unit,*,iostat=ios) (element_geom(j,i), j = 1, 11)
      if (ios /= 0) goto 701
    end do

    call read_next_line(unit, line, ios)
    if (ios /= 0) goto 700
    call lower_string(line)
    if (line /= 'element_connectivity') goto 700
    do i = 1, num_elements
      read(unit,*,iostat=ios) (element_connectivity(j,i), j = 1, 9)
      if (ios /= 0) goto 701
    end do

    call read_next_line(unit, line, ios)
    if (ios /= 0) goto 700
    call lower_string(line)
    if (line /= 'edge_connectivity') goto 700
    do i = 1, num_edges
      read(unit,*,iostat=ios) (edge_connectivity(j,i), j = 1, 6)
      if (ios /= 0) goto 701
    end do

    call read_next_line(unit, line, ios)
    if (ios /= 0) goto 702
    call lower_string(line)
    if (line /= '[end_mesh]') goto 702

    num_state_entries = num_equations * num_elements
    return

700 continue
    write(error_unit,*) 'Error: malformed mesh section in case.fvc'
    error stop
701 continue
    write(error_unit,*) 'Error: malformed mesh data in case.fvc'
    error stop
702 continue
    write(error_unit,*) 'Error: missing [end_mesh] in case.fvc'
    error stop
  end subroutine read_mesh_section

  subroutine map_boundary_edges()
    integer :: i

    do i = 1, num_edges
      if (edge_connectivity(4,i) == 0) then
        if (edge_connectivity(6,i) == 0) then
          write(error_unit,*) 'Error: boundary edge missing phys tag'
          error stop
        end if
        call find_bc_index(edge_connectivity(6,i), edge_connectivity(5,i))
        if (edge_connectivity(5,i) == 0) then
          write(error_unit,*) 'Error: no BC definition for phys tag:', edge_connectivity(6,i)
          error stop
        end if
      else
        edge_connectivity(5,i) = 0
      end if
    end do
  end subroutine map_boundary_edges

  subroutine parse_kv_line(line, key, vals, ios)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key, vals
    integer, intent(out) :: ios
    integer :: p

    p = index(line, ':')
    if (p <= 0) then
      ios = 1
      return
    end if

    key = adjustl(line(1:p-1))
    vals = adjustl(line(p+1:))
    call lower_string(key)
    call lower_string(vals)
    key = trim(key)
    vals = trim(vals)
    ios = 0
  end subroutine parse_kv_line

  subroutine parse_bc_entry(line, ios)
    character(len=*), intent(in) :: line
    integer, intent(out) :: ios
    character(len=128) :: lhs, rhs, bctype, args
    integer :: p, q, idx
    real(real64) :: v1, v2, v3, v4
    integer :: i

    ios = 0
    p = index(line, '=')
    if (p <= 0) then
      ios = 1
      return
    end if
    lhs = adjustl(line(1:p-1))
    rhs = adjustl(line(p+1:))
    call lower_string(lhs)
    call lower_string(rhs)
    lhs = trim(lhs)
    rhs = trim(rhs)

    idx = 0
    do i = 1, num_bc_defs
      if (trim(lhs) == trim(bc_label(i))) then
        idx = i
        exit
      end if
    end do
    if (idx == 0) then
      ios = 2
      return
    end if
    if (bc_kind(idx) /= 0) then
      ios = 3
      return
    end if

    q = index(rhs, '(')
    if (q > 0) then
      bctype = trim(rhs(1:q-1))
      args = rhs(q+1:)
      p = index(args, ')')
      if (p > 0) args = args(1:p-1)
    else
      bctype = trim(rhs)
      args = ''
    end if

    if (bctype == 'wall') then
      bc_kind(idx) = 1
    elseif (bctype == 'inlet') then
      bc_kind(idx) = 2
    elseif (bctype == 'extrapolated') then
      bc_kind(idx) = 3
    elseif (bctype == 'farfield') then
      bc_kind(idx) = 4
    else
      ios = 4
      return
    end if

    if (bc_kind(idx) == 2 .or. bc_kind(idx) == 4) then
      if (len_trim(args) <= 0) then
        ios = 5
        return
      end if
      read(args,*,iostat=ios) v1, v2, v3, v4
      if (ios /= 0) then
        ios = 5
        return
      end if
      bc_values(1,idx) = v1
      bc_values(2,idx) = v2
      bc_values(3,idx) = v3
      bc_values(4,idx) = v4
    end if
  end subroutine parse_bc_entry

  subroutine find_bc_index(tag, idx)
    integer, intent(in) :: tag
    integer, intent(out) :: idx
    integer :: i

    idx = 0
    do i = 1, num_bc_defs
      if (bc_physical_tag(i) == tag) then
        idx = i
        return
      end if
    end do
  end subroutine find_bc_index
end module io
