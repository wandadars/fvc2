!------------------------------------------------------------------------------
! Core data for the finite-volume solver (mesh, solution, and controls).
!------------------------------------------------------------------------------
module data
  use iso_fortran_env, only: real64
  implicit none

  ! Mesh sizing limits (static arrays for simplicity/teaching).
  integer, parameter :: max_mesh_items = 200000   ! max nodes/elements supported
  integer, parameter :: num_equations = 4         ! number of conserved equations
  integer, parameter :: max_state_entries = max_mesh_items * num_equations
  integer, parameter :: max_boundary_conditions = 20      ! total number of boundary conditions

  ! Grid sizes
  integer :: num_elements = 0, num_nodes = 0, num_edges = 0
  integer :: num_state_entries = 0

  ! Conservative state per element (rho, rho*u, rho*v, rho*E)
  real(real64) :: state(max_mesh_items, num_equations) = 0.0_real64
  real(real64) :: state_next(max_mesh_items, num_equations) = 0.0_real64
  real(real64) :: state_backup(max_mesh_items, num_equations) = 0.0_real64

  ! RHS vector (finite-volume residual accumulation)
  real(real64) :: rhs(max_state_entries) = 0.0_real64

  ! Mesh: node coordinates and per-element geometric data
  ! element_geom(:,:) stores centroid, area, and outward edge normals.
  real(real64) :: node_coords(2, max_mesh_items) = 0.0_real64
  real(real64) :: element_geom(11, max_mesh_items) = 0.0_real64

  ! Mesh connectivity (elements and edges)
  integer :: element_connectivity(9, max_mesh_items) = 0
  integer :: edge_connectivity(6, max_mesh_items) = 0

  ! Boundary-condition definitions (from case.vars)
  real(real64) :: bc_values(4, max_boundary_conditions) = 0.0_real64
  integer :: bc_kind(max_boundary_conditions) = 0
  integer :: bc_physical_tag(max_boundary_conditions) = 0
  integer :: num_bc_defs = 0
  character(len=64) :: bc_label(max_boundary_conditions) = ''

  ! Per-edge geometric vectors used during flux evaluation
  real(real64) :: dir_to_neighbor_el1(2) = 0.0_real64
  real(real64) :: dir_to_neighbor_el2(2) = 0.0_real64
  real(real64) :: edge_tangent(2) = 0.0_real64
  real(real64) :: centroid_distance = 0.0_real64
  real(real64) :: edge_length = 0.0_real64
  real(real64) :: normal_dot_dir_el1 = 0.0_real64
  real(real64) :: normal_dot_dir_el2 = 0.0_real64

  ! Precomputed edge geometry (arrays indexed by edge_id)
  real(real64) :: edge_length_arr(max_mesh_items) = 0.0_real64
  real(real64) :: centroid_distance_arr(max_mesh_items) = 0.0_real64
  real(real64) :: edge_tangent_arr(2, max_mesh_items) = 0.0_real64
  real(real64) :: dir_to_neighbor_left(2, max_mesh_items) = 0.0_real64
  real(real64) :: dir_to_neighbor_right(2, max_mesh_items) = 0.0_real64
  real(real64) :: normal_dot_dir_left(max_mesh_items) = 0.0_real64
  real(real64) :: normal_dot_dir_right(max_mesh_items) = 0.0_real64
  integer :: edge_side_left(max_mesh_items) = 0
  integer :: edge_side_right(max_mesh_items) = 0

  ! Scratch indices for the current edge/element pair
  integer :: edge_id = 0, elem_left = 0, elem_right = 0, side_index_left = 0, side_index_right = 0
  integer :: elem_left_offset = 0, elem_right_offset = 0
  integer :: elem_left_rho_idx = 0, elem_left_rhou_idx = 0, elem_left_rhov_idx = 0, elem_left_rhoe_idx = 0
  integer :: elem_right_rho_idx = 0, elem_right_rhou_idx = 0, elem_right_rhov_idx = 0, elem_right_rhoe_idx = 0

  ! Physical-name info (from case.fvc)
  integer :: ntags = 0
  integer :: ntags_dim(max_boundary_conditions) = 0
  integer :: ntags_num(max_boundary_conditions) = 0
  integer :: ntags_type(max_boundary_conditions) = 0

  ! Boundary edges count
  integer :: nbc = 0

  ! Scheme controls
  integer :: convection_scheme = 0, time_integration_scheme = 0, diffusion_scheme = 0
  integer :: output_interval = 0, max_steps = 0

  ! Gas properties
  real(real64) :: gamma_gas = 0.0_real64
  real(real64) :: gas_constant = 0.0_real64
  real(real64) :: cp_gas = 0.0_real64
  real(real64) :: cv_gas = 0.0_real64
  real(real64) :: diffusion_coeff = 0.0_real64

  ! Time integration state
  real(real64) :: time_step = 0.0_real64
  real(real64) :: simulation_time = 0.0_real64

  ! RK4 storage (placeholder for future extension)
  real(real64) :: rk4_storage(max_state_entries, 3) = 0.0_real64

  ! Initial-condition controls
  integer :: initial_condition_type = 0
  real(real64) :: ic_uniform_state(4) = 0.0_real64
  real(real64) :: shock_split_x = 0.0_real64
  real(real64) :: shock_left_state(4) = 0.0_real64
  real(real64) :: shock_right_state(4) = 0.0_real64

  ! Output format selector
  integer :: output_format_mode = 0
  logical :: animate_output = .false.
contains
  subroutine reset_arrays()
    ! Clear mesh and solution buffers before reading a case.
    node_coords = 0.0_real64
    element_geom = 0.0_real64
    element_connectivity = 0
    edge_connectivity = 0
    state = 0.0_real64
    state_next = 0.0_real64
    state_backup = 0.0_real64
    rhs = 0.0_real64
    edge_length_arr = 0.0_real64
    centroid_distance_arr = 0.0_real64
    edge_tangent_arr = 0.0_real64
    dir_to_neighbor_left = 0.0_real64
    dir_to_neighbor_right = 0.0_real64
    normal_dot_dir_left = 0.0_real64
    normal_dot_dir_right = 0.0_real64
    edge_side_left = 0
    edge_side_right = 0
  end subroutine reset_arrays
end module data
