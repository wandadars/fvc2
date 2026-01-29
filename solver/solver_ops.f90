!------------------------------------------------------------------------------
! Core solver operations (residual build and simulation_time stepping).
!------------------------------------------------------------------------------
module solver_ops
  use iso_fortran_env, only: real64
  use data
  use user, only: user_source
  implicit none
contains
  subroutine precompute_edge_geometry()
    ! Precompute per-edge geometry and side indices (static mesh).
    integer :: edge_idx, elem_idx, side_idx, nsides
    integer :: node_a, node_b
    real(real64) :: x1, y1, x2, y2
    real(real64) :: xc_left, yc_left, xc_right, yc_right
    real(real64) :: dx, dy

    do edge_idx = 1, num_edges
      elem_left = edge_connectivity(3, edge_idx)
      elem_right = edge_connectivity(4, edge_idx)

      edge_side_left(edge_idx) = 0
      edge_side_right(edge_idx) = 0

      if (elem_left > 0) then
        nsides = element_connectivity(9, elem_left)
        do side_idx = 0, nsides - 1
          if (edge_idx == element_connectivity(5 + side_idx, elem_left)) then
            edge_side_left(edge_idx) = side_idx
            exit
          end if
        end do
      end if

      if (elem_right > 0) then
        nsides = element_connectivity(9, elem_right)
        do side_idx = 0, nsides - 1
          if (edge_idx == element_connectivity(5 + side_idx, elem_right)) then
            edge_side_right(edge_idx) = side_idx
            exit
          end if
        end do
      end if

      node_a = edge_connectivity(1, edge_idx)
      node_b = edge_connectivity(2, edge_idx)
      x1 = node_coords(1, node_a)
      y1 = node_coords(2, node_a)
      x2 = node_coords(1, node_b)
      y2 = node_coords(2, node_b)
      dx = x2 - x1
      dy = y2 - y1
      edge_length_arr(edge_idx) = sqrt(dx*dx + dy*dy)
      if (edge_length_arr(edge_idx) > 0.0_real64) then
        edge_tangent_arr(1, edge_idx) = dx / edge_length_arr(edge_idx)
        edge_tangent_arr(2, edge_idx) = dy / edge_length_arr(edge_idx)
      else
        edge_tangent_arr(:, edge_idx) = 0.0_real64
      end if

      if (elem_right > 0) then
        xc_left = element_geom(1, elem_left)
        yc_left = element_geom(2, elem_left)
        xc_right = element_geom(1, elem_right)
        yc_right = element_geom(2, elem_right)
        dx = xc_right - xc_left
        dy = yc_right - yc_left
        centroid_distance_arr(edge_idx) = sqrt(dx*dx + dy*dy)
        if (centroid_distance_arr(edge_idx) > 0.0_real64) then
          dir_to_neighbor_left(1, edge_idx) = dx / centroid_distance_arr(edge_idx)
          dir_to_neighbor_left(2, edge_idx) = dy / centroid_distance_arr(edge_idx)
          dir_to_neighbor_right(1, edge_idx) = -dir_to_neighbor_left(1, edge_idx)
          dir_to_neighbor_right(2, edge_idx) = -dir_to_neighbor_left(2, edge_idx)
        else
          dir_to_neighbor_left(:, edge_idx) = 0.0_real64
          dir_to_neighbor_right(:, edge_idx) = 0.0_real64
        end if
      else
        xc_left = element_geom(1, elem_left)
        yc_left = element_geom(2, elem_left)
        dx = 0.5_real64*(x1 + x2) - xc_left
        dy = 0.5_real64*(y1 + y2) - yc_left
        centroid_distance_arr(edge_idx) = sqrt(dx*dx + dy*dy)
        if (centroid_distance_arr(edge_idx) > 0.0_real64) then
          dir_to_neighbor_left(1, edge_idx) = dx / centroid_distance_arr(edge_idx)
          dir_to_neighbor_left(2, edge_idx) = dy / centroid_distance_arr(edge_idx)
        else
          dir_to_neighbor_left(:, edge_idx) = 0.0_real64
        end if
        dir_to_neighbor_right(:, edge_idx) = 0.0_real64
      end if

      if (elem_left > 0) then
        side_idx = edge_side_left(edge_idx)
        normal_dot_dir_left(edge_idx) = element_geom(4 + side_idx*2, elem_left) * dir_to_neighbor_left(1, edge_idx) + &
                                        element_geom(5 + side_idx*2, elem_left) * dir_to_neighbor_left(2, edge_idx)
      else
        normal_dot_dir_left(edge_idx) = 0.0_real64
      end if

      if (elem_right > 0) then
        side_idx = edge_side_right(edge_idx)
        normal_dot_dir_right(edge_idx) = element_geom(4 + side_idx*2, elem_right) * dir_to_neighbor_right(1, edge_idx) + &
                                         element_geom(5 + side_idx*2, elem_right) * dir_to_neighbor_right(2, edge_idx)
      else
        normal_dot_dir_right(edge_idx) = 0.0_real64
      end if
    end do
  end subroutine precompute_edge_geometry

  subroutine zero_matricies(bdum, istep)
    real(real64), intent(inout) :: bdum(num_state_entries)
    integer, intent(in) :: istep
    integer :: j

    do j = 1, num_state_entries
      bdum(j) = 0.0_real64
    end do
  end subroutine zero_matricies

  subroutine solve_equations_fwdE(bdum, istage)
    real(real64), intent(inout) :: bdum(num_state_entries)
    integer, intent(in) :: istage
    integer :: j, jel
    integer :: i, ii, ic

    do j = 1, num_state_entries
      jel = (j - 1) / num_equations + 1
      bdum(j) = bdum(j) / element_geom(3, jel) * time_step
    end do

    ic = 0
    do i = 1, num_elements
      do ii = 1, num_equations
        ic = ic + 1
        state(i, ii) = state(i, ii) + bdum(ic)
      end do
    end do
  end subroutine solve_equations_fwdE

  subroutine compute_sources(bdum)
    real(real64), intent(inout) :: bdum(num_state_entries)
    real(real64) :: rsource1(num_equations), rsource2(num_equations)
    integer :: j, i, iems1, iems2

    do j = 1, num_equations
      rsource1(j) = 0.0_real64
      rsource2(j) = 0.0_real64
    end do

    call user_source(element_geom(1, elem_left), element_geom(2, elem_left), rsource1)
    do i = 0, num_equations - 1
      iems1 = elem_left_offset + i
      rsource1(i+1) = rsource1(i+1) / 4.0_real64
      bdum(iems1) = bdum(iems1) + rsource1(i+1) * element_geom(3, elem_left)
    end do

    if (edge_connectivity(5, edge_id) == 0) then
      call user_source(element_geom(1, elem_right), element_geom(2, elem_right), rsource2)
      do i = 0, num_equations - 1
        iems2 = elem_right_offset + i
        rsource2(i+1) = rsource2(i+1) / 4.0_real64
        bdum(iems2) = bdum(iems2) + rsource2(i+1) * element_geom(3, elem_right)
      end do
    end if
  end subroutine compute_sources

  subroutine compute_edge_values()

    elem_left = edge_connectivity(3, edge_id)
    elem_left_offset = 1 + num_equations * (elem_left - 1)

    elem_left_rho_idx  = elem_left_offset + 0
    elem_left_rhou_idx = elem_left_offset + 1
    elem_left_rhov_idx = elem_left_offset + 2
    elem_left_rhoe_idx = elem_left_offset + 3

    side_index_left = edge_side_left(edge_id)

    if (edge_connectivity(5, edge_id) == 0) then
      elem_right = edge_connectivity(4, edge_id)
      elem_right_offset = 1 + num_equations * (elem_right - 1)

      elem_right_rho_idx  = elem_right_offset + 0
      elem_right_rhou_idx = elem_right_offset + 1
      elem_right_rhov_idx = elem_right_offset + 2
      elem_right_rhoe_idx = elem_right_offset + 3

      side_index_right = edge_side_right(edge_id)

      centroid_distance = centroid_distance_arr(edge_id)
      dir_to_neighbor_el1(1) = dir_to_neighbor_left(1, edge_id)
      dir_to_neighbor_el1(2) = dir_to_neighbor_left(2, edge_id)
      dir_to_neighbor_el2(1) = dir_to_neighbor_right(1, edge_id)
      dir_to_neighbor_el2(2) = dir_to_neighbor_right(2, edge_id)

      edge_length = edge_length_arr(edge_id)
      edge_tangent(1) = edge_tangent_arr(1, edge_id)
      edge_tangent(2) = edge_tangent_arr(2, edge_id)

      normal_dot_dir_el1 = normal_dot_dir_left(edge_id)
      normal_dot_dir_el2 = normal_dot_dir_right(edge_id)

    else
      centroid_distance = centroid_distance_arr(edge_id)
      dir_to_neighbor_el1(1) = dir_to_neighbor_left(1, edge_id)
      dir_to_neighbor_el1(2) = dir_to_neighbor_left(2, edge_id)
      dir_to_neighbor_el2(1) = 0.0_real64
      dir_to_neighbor_el2(2) = 0.0_real64

      edge_length = edge_length_arr(edge_id)
      edge_tangent(1) = edge_tangent_arr(1, edge_id)
      edge_tangent(2) = edge_tangent_arr(2, edge_id)

      normal_dot_dir_el1 = normal_dot_dir_left(edge_id)
      normal_dot_dir_el2 = 0.0_real64
    end if
  end subroutine compute_edge_values
end module solver_ops
