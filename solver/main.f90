program fvc2_solver_main
  use iso_fortran_env, only: error_unit, real64
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
    call advance_one_timestep(istep)
  end do

contains
  subroutine advance_one_timestep(istep)
    integer, intent(in) :: istep

    select case (time_integration_scheme)
    case (1)
      call advance_forward_euler(istep)
    case (2)
      call advance_rk4(istep)
    case default
      write(error_unit,*) 'Error: unsupported time_discretization value: ', time_integration_scheme
      error stop
    end select
  end subroutine advance_one_timestep

  subroutine advance_forward_euler(istep)
    integer, intent(in) :: istep

    ! Build the semi-discrete residual for the current state and apply
    ! a single forward-Euler update: U^(n+1) = U^n + dt * dU/dt(U^n).
    call assemble_convective_rhs(istep, rhs)
    call solve_equations_fwdE(rhs, 1)
  end subroutine advance_forward_euler

  subroutine assemble_convective_rhs(istep, rhs_vec)
    integer, intent(in) :: istep            ! Current time-step index (passed to zero_matricies)
    real(real64), intent(inout) :: rhs_vec(:) ! Residual vector to fill from edge fluxes

    ! Important: this routine is "state-local", not "time-aware".
    ! It does not know whether we are at U^n, U^(2), U^(3), or U^(4).
    ! It simply evaluates the residual for whatever conservative state is
    ! currently stored in the global array `state`.
    !
    ! Convert conservative state -> primitive fields (rho,u,v,p,a,h) for all cells.
    call update_primitive_fields()
    ! Reset residual accumulation before looping over edges.
    call zero_matricies(rhs_vec, istep)

    ! Accumulate net convective flux on every edge into rhs_vec.
    do edge_id = 1, num_edges
      ! Load left/right element indices, geometric normal, edge length, and offsets.
      call compute_edge_values()

      if (edge_connectivity(5, edge_id) == 0) then
        ! Interior edge: two-sided numerical flux contribution.
        call convective_flux_interior(rhs_vec)
      else
        ! Boundary edge: one-sided flux using BC ghost-state logic.
        call convective_flux_boundary(rhs_vec)
      end if
    end do
  end subroutine assemble_convective_rhs

  !------------------------------------------------------------------------------
  ! Build one time step using RK4.
  !
  ! Short RK4 summary:
  ! Evaluate four slopes (k1, k2, k3, k4) using the same residual operator
  ! at four states (current + three predictor states), then form U^(n+1) from
  ! a weighted combination of those slopes.
  !
  !------------------------------------------------------------------------------
  subroutine advance_rk4(istep)
    integer, intent(in) :: istep
    real(real64) :: stage2_factor   ! alpha2 in U^(2) = U^n + alpha2*dt*k1
    real(real64) :: stage3_factor   ! alpha3 in U^(3) = U^n + alpha3*dt*k2
    real(real64) :: stage4_factor   ! alpha4 in U^(4) = U^n + alpha4*dt*k3
    real(real64) :: w1, w2, w3, w4  ! Final RK weights in U^(n+1) = U^n + dt*sum(wi*ki)

    !--------------------------------------------------------------------------
    ! Stage A: save U^n so each RK stage is built from a common baseline.
    !--------------------------------------------------------------------------
    ! state_backup holds U^n for all cells/equations.
    state_backup = state

    !--------------------------------------------------------------------------
    ! Stage B: k1 = dU/dt(U^n)
    !--------------------------------------------------------------------------
    ! assemble_convective_rhs computes the finite-volume residual R(U^n)
    ! from all interior/boundary edge fluxes.
    call assemble_convective_rhs(istep, rhs)

    ! residual_to_time_derivative_inplace converts R(U) to dU/dt by dividing
    ! each equation entry by its control-volume area.
    call residual_to_time_derivative_inplace(rhs)

    ! rk4_storage(:,1) stores k1.
    rk4_storage(:,1) = rhs

    ! TODO(task 1): choose predictor factor for stage-2 state.
    ! (Build U^(2) = U^n + stage2_factor * dt * k1.)
    ! (Classical RK4 uses a half step here.)
    stage2_factor = 0.0_real64  ! REPLACE with correct value

    ! apply_stage_predictor writes state = U^n + stage2_factor*dt*k1.
    call apply_stage_predictor(state_backup, rk4_storage(:,1), stage2_factor)

    ! Stage C: k2 = dU/dt(U^(2))
    call assemble_convective_rhs(istep, rhs)
    call residual_to_time_derivative_inplace(rhs)

    ! rk4_storage(:,2) stores k2.
    rk4_storage(:,2) = rhs

    ! TODO(task 2): choose predictor factor for stage-3 state.
    ! (Build U^(3) = U^n + stage3_factor * dt * k2.)
    ! (Classical RK4 also uses a half step here.)
    stage3_factor = 0.0_real64  ! REPLACE with correct value

    ! apply_stage_predictor writes state = U^n + stage3_factor*dt*k2.
    call apply_stage_predictor(state_backup, rk4_storage(:,2), stage3_factor)

    ! Stage D: k3 = dU/dt(U^(3))
    call assemble_convective_rhs(istep, rhs)
    call residual_to_time_derivative_inplace(rhs)

    ! rk4_storage(:,3) stores k3.
    rk4_storage(:,3) = rhs

    ! TODO(task 3): choose predictor factor for stage-4 state.
    ! (Build U^(4) = U^n + stage4_factor * dt * k3.)
    stage4_factor = 0.0_real64  ! REPLACE with correct value

    ! apply_stage_predictor writes state = U^n + stage4_factor*dt*k3.
    call apply_stage_predictor(state_backup, rk4_storage(:,3), stage4_factor)

    ! Stage E: k4 = dU/dt(U^(4)).
    ! Reuse rhs as temporary k4 storage to avoid another full-size array.
    call assemble_convective_rhs(istep, rhs)
    call residual_to_time_derivative_inplace(rhs)

    ! TODO(task 4): choose final RK4 combination weights.
    ! (Combine U^(n+1) = U^n + dt * (w1*k1 + w2*k2 + w3*k3 + w4*k4).)
    w1 = 0.0_real64  ! REPLACE with correct value
    w2 = 0.0_real64  ! REPLACE with correct value
    w3 = 0.0_real64  ! REPLACE with correct value
    w4 = 0.0_real64  ! REPLACE with correct value

    ! combine_weighted_stage_slopes applies the final RK weighted sum.
    call combine_weighted_stage_slopes(state_backup, rk4_storage(:,1), rk4_storage(:,2), &
                                       rk4_storage(:,3), rhs, w1, w2, w3, w4)
  end subroutine advance_rk4

  subroutine residual_to_time_derivative_inplace(vec)
    real(real64), intent(inout) :: vec(:) ! Input: residual R(U); output: slope dU/dt
    integer :: j                          ! Flat equation index over all cells
    integer :: jel                        ! Cell index that owns vec(j)

    ! Finite-volume relation: dU/dt = R(U) / cell_area
    do j = 1, num_state_entries
      jel = (j - 1) / num_equations + 1
      vec(j) = vec(j) / element_geom(3, jel)
    end do
  end subroutine residual_to_time_derivative_inplace

  subroutine apply_stage_predictor(base_state, stage_slope, stage_factor)
    real(real64), intent(in) :: base_state(:,:) ! U^n, shape (num_elements, num_equations)
    real(real64), intent(in) :: stage_slope(:)  ! k_i slope used to form the next predictor state
    real(real64), intent(in) :: stage_factor    ! alpha_i predictor coefficient
    integer :: i                                ! Cell index
    integer :: ii                               ! Equation index within cell
    integer :: ic                               ! Flat vector index matching stage_slope layout

    ! Predictor form: U_stage = U^n + alpha_i * dt * k_i
    ic = 0
    do i = 1, num_elements
      do ii = 1, num_equations
        ic = ic + 1
        state(i, ii) = base_state(i, ii) + stage_factor * time_step * stage_slope(ic)
      end do
    end do
  end subroutine apply_stage_predictor

  subroutine combine_weighted_stage_slopes(base_state, k1, k2, k3, k4, w1, w2, w3, w4)
    real(real64), intent(in) :: base_state(:,:) ! U^n, shape (num_elements, num_equations)
    real(real64), intent(in) :: k1(:), k2(:), k3(:), k4(:) ! RK slopes dU/dt at stages 1..4
    real(real64), intent(in) :: w1, w2, w3, w4             ! RK combination weights
    integer :: i                                            ! Cell index
    integer :: ii                                           ! Equation index within cell
    integer :: ic                                           ! Flat vector index for k-arrays

    ! Final RK update: U^(n+1) = U^n + dt * (w1*k1 + w2*k2 + w3*k3 + w4*k4)
    ic = 0
    do i = 1, num_elements
      do ii = 1, num_equations
        ic = ic + 1
        state(i, ii) = base_state(i, ii) + time_step * (w1*k1(ic) + w2*k2(ic) + w3*k3(ic) + w4*k4(ic))
      end do
    end do
  end subroutine combine_weighted_stage_slopes
end program fvc2_solver_main
