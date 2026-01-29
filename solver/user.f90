!------------------------------------------------------------------------------
! User-controlled initial conditions and source terms.
!------------------------------------------------------------------------------
module user
  use iso_fortran_env, only: real64
  use data
  use thermo, only: MixtPerf_et
  implicit none
contains
  subroutine user_ic()
    if (initial_condition_type == 2) then
      call user_ic_fwd_step()
    elseif (initial_condition_type == 3) then
      call user_ic_shocktube()
    else
      call user_ic_uniform()
    end if
  end subroutine user_ic

  subroutine user_ic_uniform()
    ! Uniform freestream initial condition from case.vars (ic_uniform).
    real(real64) :: conservative_state(4)
    real(real64) :: density, pressure, velocity_x, velocity_y
    real(real64) :: total_specific_energy

    density = ic_uniform_state(1)
    pressure = ic_uniform_state(2)
    velocity_x = ic_uniform_state(3)
    velocity_y = ic_uniform_state(4)

    conservative_state(1) = density
    conservative_state(2) = density * velocity_x
    conservative_state(3) = density * velocity_y
    conservative_state(4) = -1.0_real64
    total_specific_energy = MixtPerf_et(conservative_state(1), conservative_state(2), &
                                        conservative_state(3), conservative_state(4), pressure)

    state(:,1) = density
    state(:,2) = density * velocity_x
    state(:,3) = density * velocity_y
    state(:,4) = density * total_specific_energy
  end subroutine user_ic_uniform

  subroutine user_ic_fwd_step()
    ! Legacy forward-step initial condition (hard-coded).
    real(real64) :: conservative_state(4)
    real(real64) :: density, pressure, velocity_x, velocity_y
    real(real64) :: total_specific_energy

    density = 1.4_real64
    pressure = 1.0_real64
    velocity_x = 3.0_real64
    velocity_y = 0.0_real64

    conservative_state(1) = density
    conservative_state(2) = density * velocity_x
    conservative_state(3) = density * velocity_y
    conservative_state(4) = -1.0_real64
    total_specific_energy = MixtPerf_et(conservative_state(1), conservative_state(2), &
                                        conservative_state(3), conservative_state(4), pressure)

    state(:,1) = density
    state(:,2) = density * velocity_x
    state(:,3) = density * velocity_y
    state(:,4) = density * total_specific_energy
  end subroutine user_ic_fwd_step

  subroutine user_ic_shocktube()
    ! Shock-tube initial condition (left/right states split by x location).
    real(real64) :: conservative_state(4)
    real(real64) :: density, pressure, velocity_x, velocity_y
    real(real64) :: total_specific_energy
    real(real64) :: centroid_x
    integer :: elem_idx

    do elem_idx = 1, num_elements
      centroid_x = element_geom(1, elem_idx)

      if (centroid_x <= shock_split_x) then
        density = shock_left_state(1)
        pressure = shock_left_state(2)
        velocity_x = shock_left_state(3)
        velocity_y = shock_left_state(4)
      else
        density = shock_right_state(1)
        pressure = shock_right_state(2)
        velocity_x = shock_right_state(3)
        velocity_y = shock_right_state(4)
      end if

      conservative_state(1) = density
      conservative_state(2) = density * velocity_x
      conservative_state(3) = density * velocity_y
      conservative_state(4) = -1.0_real64
      total_specific_energy = MixtPerf_et(conservative_state(1), conservative_state(2), &
                                          conservative_state(3), conservative_state(4), pressure)

      state(elem_idx,1) = density
      state(elem_idx,2) = density * velocity_x
      state(elem_idx,3) = density * velocity_y
      state(elem_idx,4) = density * total_specific_energy
    end do
  end subroutine user_ic_shocktube

  subroutine user_source(x_coord, y_coord, source_term)
    ! User-defined source term (placeholder, currently zero).
    real(real64), intent(in) :: x_coord, y_coord
    real(real64), intent(out) :: source_term(2)

    source_term(1) = 0.0_real64
    source_term(2) = 0.0_real64
  end subroutine user_source
end module user
