program fvc2_solver_main
  use iso_fortran_env, only: output_unit
  use data
  use io
  use user
  use solver_ops
  use convection
  implicit none

  integer :: istep

  call reset_arrays()

  ! Read mesh and runtime parameters, then initialize the flow field.
  call read_case_file()
  call precompute_edge_geometry()
  call user_ic()

  simulation_time = 0.0_real64
  do istep = 0, max_steps
    if (modulo(istep, output_interval) == 0) then
      call write_soln()
    end if

    simulation_time = simulation_time + time_step

    ! Forward Euler update (default). RK4 placeholder kept for teaching.
    if (time_integration_scheme == 1) then
      call update_primitive_fields()
      call zero_matricies(rhs, istep)

      do edge_id = 1, num_edges
        call compute_edge_values()

        if (edge_connectivity(5, edge_id) == 0) then
          call convective_flux_interior(rhs)
        else
          call convective_flux_boundary(rhs)
        end if
      end do

      call solve_equations_fwdE(rhs, 1)
    elseif (time_integration_scheme == 2) then
      ! RK4 - code here
    end if
  end do
end program fvc2_solver_main
