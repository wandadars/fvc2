!------------------------------------------------------------------------------
! Convective fluxes (AUSM+ and Roe) for the Euler equations.
!------------------------------------------------------------------------------
module convection
  use iso_fortran_env, only: real64
  use data
  implicit none
contains
  subroutine convective_flux_interior(bdum)
    real(real64), intent(inout) :: bdum(num_state_entries)
    real(real64) :: flux(4), nhatx, nhaty

    nhatx = element_geom(4+side_index_left*2, elem_left)
    nhaty = element_geom(5+side_index_left*2, elem_left)

    call convective_flux_cached(elem_left, elem_right, nhatx, nhaty, flux)

    bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
    bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
    bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
    bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    bdum(elem_right_rho_idx)  = bdum(elem_right_rho_idx)  + flux(1) * edge_length
    bdum(elem_right_rhou_idx) = bdum(elem_right_rhou_idx) + flux(2) * edge_length
    bdum(elem_right_rhov_idx) = bdum(elem_right_rhov_idx) + flux(3) * edge_length
    bdum(elem_right_rhoe_idx) = bdum(elem_right_rhoe_idx) + flux(4) * edge_length
  end subroutine convective_flux_interior

  subroutine convective_flux_boundary(bdum)
    real(real64), intent(inout) :: bdum(num_state_entries)
    real(real64) :: flux(4), rnhatx, rnhaty
    real(real64) :: rdens, rpres, ru, rv, ret
    real(real64) :: uvel, vvel, vn
    real(real64) :: rho_l, u_l, v_l, p_l, a_l, h_l
    real(real64) :: rho_r, u_r, v_r, p_r, a_r, h_r
    integer :: idum, ibc

    idum = edge_connectivity(5, edge_id)
    if (idum <= 0) return
    ibc = bc_kind(idum)

    if (ibc == 3) then
      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      rho_l = density_arr(elem_left)
      u_l = velocity_x_arr(elem_left)
      v_l = velocity_y_arr(elem_left)
      p_l = pressure_arr(elem_left)
      a_l = sound_speed_arr(elem_left)
      h_l = enthalpy_arr(elem_left)

      call flux_from_primitives(rho_l, u_l, v_l, p_l, a_l, h_l, &
                                rho_l, u_l, v_l, p_l, a_l, h_l, &
                                rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    elseif (ibc == 2) then
      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      rdens = bc_values(bc_value_rho_idx, idum)
      rpres = bc_values(bc_value_p_idx, idum)
      ru    = bc_values(bc_value_u_idx, idum)
      rv    = bc_values(bc_value_v_idx, idum)
      if (rdens <= 0.0_real64 .or. rpres <= 0.0_real64) then
        rdens = ic_uniform_state(1)
        rpres = ic_uniform_state(2)
        ru    = ic_uniform_state(3)
        rv    = ic_uniform_state(4)
      end if

      rho_l = density_arr(elem_left)
      u_l = velocity_x_arr(elem_left)
      v_l = velocity_y_arr(elem_left)
      p_l = pressure_arr(elem_left)
      a_l = sound_speed_arr(elem_left)
      h_l = enthalpy_arr(elem_left)

      rho_r = rdens
      u_r = ru
      v_r = rv
      p_r = rpres
      a_r = sqrt(gamma_gas * p_r / rho_r)
      ret = p_r/(gamma_gas-1.0_real64)/rho_r + 0.5_real64*(u_r*u_r + v_r*v_r)
      h_r = ret + p_r / rho_r

      call flux_from_primitives(rho_l, u_l, v_l, p_l, a_l, h_l, &
                                rho_r, u_r, v_r, p_r, a_r, h_r, &
                                rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    elseif (ibc == 4) then
      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      rho_l = density_arr(elem_left)
      u_l = velocity_x_arr(elem_left)
      v_l = velocity_y_arr(elem_left)
      p_l = pressure_arr(elem_left)
      a_l = sound_speed_arr(elem_left)
      h_l = enthalpy_arr(elem_left)

      uvel = u_l
      vvel = v_l
      vn = uvel * rnhatx + vvel * rnhaty

      if (vn >= 0.0_real64) then
        rho_r = rho_l
        u_r = u_l
        v_r = v_l
        p_r = p_l
        a_r = a_l
        h_r = h_l
      else
        rdens = bc_values(bc_value_rho_idx, idum)
        rpres = bc_values(bc_value_p_idx, idum)
        ru    = bc_values(bc_value_u_idx, idum)
        rv    = bc_values(bc_value_v_idx, idum)
        if (rdens <= 0.0_real64 .or. rpres <= 0.0_real64) then
          rdens = ic_uniform_state(1)
          rpres = ic_uniform_state(2)
          ru    = ic_uniform_state(3)
          rv    = ic_uniform_state(4)
        end if
        rho_r = rdens
        u_r = ru
        v_r = rv
        p_r = rpres
        a_r = sqrt(gamma_gas * p_r / rho_r)
        ret = p_r/(gamma_gas-1.0_real64)/rho_r + 0.5_real64*(u_r*u_r + v_r*v_r)
        h_r = ret + p_r / rho_r
      end if

      call flux_from_primitives(rho_l, u_l, v_l, p_l, a_l, h_l, &
                                rho_r, u_r, v_r, p_r, a_r, h_r, &
                                rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    elseif (ibc == 5) then
      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      rho_l = density_arr(elem_left)
      u_l = velocity_x_arr(elem_left)
      v_l = velocity_y_arr(elem_left)
      p_l = pressure_arr(elem_left)
      a_l = sound_speed_arr(elem_left)
      h_l = enthalpy_arr(elem_left)

      call build_subsonic_outlet_ghost(idum, rnhatx, rnhaty, rho_l, u_l, v_l, p_l, a_l, h_l, &
                                       rho_r, u_r, v_r, p_r, a_r, h_r)

      call flux_from_primitives(rho_l, u_l, v_l, p_l, a_l, h_l, &
                                rho_r, u_r, v_r, p_r, a_r, h_r, &
                                rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    elseif (ibc == 1) then
      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      rho_l = density_arr(elem_left)
      u_l = velocity_x_arr(elem_left)
      v_l = velocity_y_arr(elem_left)
      p_l = pressure_arr(elem_left)
      a_l = sound_speed_arr(elem_left)
      h_l = enthalpy_arr(elem_left)

      rho_r = rho_l
      p_r = p_l
      a_r = a_l

      if (abs(rnhaty) < 1.0e-9_real64) then
        u_r = 0.0_real64
        v_r = v_l
      elseif (abs(rnhatx) < 1.0e-9_real64) then
        u_r = u_l
        v_r = 0.0_real64
      else
        u_r = 0.0_real64
        v_r = -u_r * rnhatx / rnhaty
      end if
      ret = p_r/(gamma_gas-1.0_real64)/rho_r + 0.5_real64*(u_r*u_r + v_r*v_r)
      h_r = ret + p_r / rho_r

      call flux_from_primitives(rho_l, u_l, v_l, p_l, a_l, h_l, &
                                rho_r, u_r, v_r, p_r, a_r, h_r, &
                                rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length
    end if
  end subroutine convective_flux_boundary

  !------------------------------------------------------------------------------
  ! Build the ghost (right) state for a subsonic outlet boundary face.
  !
  ! Responsibility of this routine:
  ! 1) Read the outlet target pressure p_out from already-parsed BC storage.
  ! 2) Decide whether the face is a valid subsonic outflow face.
  ! 3) If not valid, provide a safe fallback (extrapolated ghost = interior state).
  ! 4) If valid, implement the subsonic outlet ghost-state model
  !    through TODO tasks.
  !
  ! This routine is called before flux_from_primitives(...), so the outputs
  ! rho_r/u_r/v_r/p_r/a_r/h_r must be fully defined on every path.
  !------------------------------------------------------------------------------
  subroutine build_subsonic_outlet_ghost(idum, nx, ny, rho_l, u_l, v_l, p_l, a_l, h_l, &
                                         rho_r, u_r, v_r, p_r, a_r, h_r)
    integer, intent(in) :: idum                    ! BC definition index for this boundary edge
    real(real64), intent(in) :: nx, ny             ! Outward unit normal components from interior cell
    real(real64), intent(in) :: rho_l, u_l, v_l    ! Interior density and velocity components
    real(real64), intent(in) :: p_l, a_l, h_l      ! Interior pressure(p), sound speed(a), and total enthalpy(h)
    real(real64), intent(out) :: rho_r, u_r, v_r   ! Ghost density and velocity components (to be built)
    real(real64), intent(out) :: p_r, a_r, h_r     ! Ghost pressure, sound speed, and total enthalpy
    real(real64) :: vn                             ! Interior normal velocity (u dot n)
    real(real64) :: p_out                          ! User-prescribed outlet pressure from case.vars

    !--------------------------------------------------------------------------
    ! Stage A: gather quantities needed for the BC logic.
    !--------------------------------------------------------------------------
    vn = u_l * nx + v_l * ny ! Velocity normal(perpendicular) to the boundary face, positive for outflow.

    ! Parsed from case.vars as: outlet=subsonic_outlet(p_out)
    ! Parser storage rule:
    ! bc_values(bc_value_subsonic_outlet_p_idx, idum) = p_out
    p_out = bc_values(bc_value_subsonic_outlet_p_idx, idum)

    !--------------------------------------------------------------------------
    ! Stage B: guard conditions and robust fallback.
    !--------------------------------------------------------------------------
    ! For this tutorial model, only valid subsonic outflow uses the student BC:
    ! - vn > 0       : flow exits domain through this face
    ! - |vn| < a_l   : subsonic in normal direction
    ! - positive thermodynamic quantities
    ! Any other situation (backflow, supersonic, or invalid state) falls back to
    ! extrapolated ghost values so the solver remains stable.
    if (rho_l <= 0.0_real64 .or. p_l <= 0.0_real64 .or. a_l <= 0.0_real64 .or. &
        p_out <= 0.0_real64 .or. vn <= 0.0_real64 .or. abs(vn) >= a_l) then
      rho_r = rho_l
      u_r = u_l
      v_r = v_l
      p_r = p_l
      a_r = a_l
      h_r = h_l
      return
    end if

    !--------------------------------------------------------------------------
    ! Stage C (student implementation): fill in subsonic outlet ghost state.
    !--------------------------------------------------------------------------
    ! TODO(task 1): impose outlet static pressure from BC input (p_out).
    ! (The outlet static pressure is stored in the p_out variable)
    p_r = 0.0_real64  ! REPLACE with correct value

    ! TODO(task 2): set ghost velocity from the interior-adjacent cell.
    ! (Interior velocity components are available as u_l and v_l)
    u_r = 0.0_real64  ! REPLACE with correct value
    v_r = 0.0_real64  ! REPLACE with correct value

    ! TODO(task 3): compute ghost density from an isentropic relation.
    ! (For an isentropic ideal-gas update: p / rho^gamma = constant.)
    ! (Apply that between interior state "l" and ghost state "r", then solve
    ! algebraically for rho_r using rho_l, p_l, p_r, and gamma_gas.)
    rho_r = 0.0_real64  ! REPLACE with correct value

    ! TODO(task 4): recompute ghost sound speed and total enthalpy.
    ! (Sound speed for ideal gas: a_r^2 = gamma_gas * p_r / rho_r.)
    ! (Total enthalpy here is h = e_t + p/rho, where)
    ! (e_t = p/((gamma_gas-1)*rho) + 0.5*(u^2 + v^2).)
    ! Use the ghost state values in the expressions to obtain expressions for
    ! a_r and h_r.
    a_r = 0.0_real64  ! REPLACE with correct value
    h_r = 0.0_real64  ! REPLACE with correct value

    !--------------------------------------------------------------------------
    ! Stage D (after TODOs):
    ! At return, these outputs are consumed by flux_from_primitives(...) in
    ! convective_flux_boundary to compute the boundary flux contribution.
    !--------------------------------------------------------------------------
  end subroutine build_subsonic_outlet_ghost

  subroutine convective_flux_cached(elem_l, elem_r, nhatx, nhaty, flux)
    integer, intent(in) :: elem_l, elem_r
    real(real64), intent(in) :: nhatx, nhaty
    real(real64), intent(out) :: flux(4)

    call flux_from_primitives(density_arr(elem_l), velocity_x_arr(elem_l), velocity_y_arr(elem_l), &
                              pressure_arr(elem_l), sound_speed_arr(elem_l), enthalpy_arr(elem_l), &
                              density_arr(elem_r), velocity_x_arr(elem_r), velocity_y_arr(elem_r), &
                              pressure_arr(elem_r), sound_speed_arr(elem_r), enthalpy_arr(elem_r), &
                              nhatx, nhaty, flux)
  end subroutine convective_flux_cached

  subroutine flux_from_primitives(rhoL, uL, vL, pL, aL, HL, rhoR, uR, vR, pR, aR, HR, nx, ny, flux)
    real(real64), intent(in) :: rhoL, uL, vL, pL, aL, HL
    real(real64), intent(in) :: rhoR, uR, vR, pR, aR, HR
    real(real64), intent(in) :: nx, ny
    real(real64), intent(out) :: flux(4)

    select case (convection_scheme)
    case (0)
      flux = 0.0_real64
    case (1)
      call ausm_plus_prim(rhoL, uL, vL, pL, aL, HL, rhoR, uR, vR, pR, aR, HR, nx, ny, flux)
    case (2)
      call Roe_prim(rhoL, uL, vL, pL, aL, HL, rhoR, uR, vR, pR, aR, HR, nx, ny, flux)
    end select
  end subroutine flux_from_primitives

  subroutine Roe_prim(rhoL, uL, vL, pL, aL, HL, rhoR, uR, vR, pR, aR, HR, nx, ny, flux)
    real(real64), intent(in) :: rhoL, uL, vL, pL, aL, HL
    real(real64), intent(in) :: rhoR, uR, vR, pR, aR, HR
    real(real64), intent(in) :: nx, ny
    real(real64), intent(out) :: flux(4)

    real(real64) :: gamma
    real(real64) :: zero, fifth, half, one, two
    real(real64) :: tx, ty
    real(real64) :: vxL, vxR, vyL, vyR
    real(real64) :: vnL, vnR, vtL, vtR
    real(real64) :: RT, rho, vx, vy, H, a, vn, vt
    real(real64) :: drho, dvn, dvt, dp, dV(4)
    real(real64) :: ws(4), dws(4), Rv(4,4)
    real(real64) :: fL(4), fR(4), diss(4)
    integer :: i, j

    gamma = gamma_gas
    zero = 0.0_real64
    fifth = 0.2_real64
    half = 0.5_real64
    one = 1.0_real64
    two = 2.0_real64

    tx = -ny
    ty = nx

    vxL = uL
    vyL = vL
    vnL = vxL*nx + vyL*ny
    vtL = vxL*tx + vyL*ty

    vxR = uR
    vyR = vR
    vnR = vxR*nx + vyR*ny
    vtR = vxR*tx + vyR*ty

    RT = sqrt(rhoR / rhoL)
    rho = RT * rhoL
    vx = (vxL + RT * vxR) / (one + RT)
    vy = (vyL + RT * vyR) / (one + RT)
    H = (HL + RT * HR) / (one + RT)
    a = sqrt((gamma - one) * (H - half*(vx*vx + vy*vy)))
    vn = vx*nx + vy*ny
    vt = vx*tx + vy*ty

    drho = rhoR - rhoL
    dp = pR - pL
    dvn = vnR - vnL
    dvt = vtR - vtL

    dV(1) = (dp - rho*a*dvn) / (two*a*a)
    dV(2) = rho * dvt / a
    dV(3) = drho - dp / (a*a)
    dV(4) = (dp + rho*a*dvn) / (two*a*a)

    ws(1) = abs(vn - a)
    ws(2) = abs(vn)
    ws(3) = abs(vn)
    ws(4) = abs(vn + a)

    dws(1) = fifth
    if (ws(1) < dws(1)) ws(1) = half * (ws(1)*ws(1)/dws(1) + dws(1))
    dws(4) = fifth
    if (ws(4) > dws(4)) ws(4) = half * (ws(4)*ws(4)/dws(4) + dws(4))

    Rv(1,1) = one
    Rv(2,1) = vx - a*nx
    Rv(3,1) = vy - a*ny
    Rv(4,1) = H - vn*a

    Rv(1,2) = zero
    Rv(2,2) = a*tx
    Rv(3,2) = a*ty
    Rv(4,2) = vt*a

    Rv(1,3) = one
    Rv(2,3) = vx
    Rv(3,3) = vy
    Rv(4,3) = half*(vx*vx + vy*vy)

    Rv(1,4) = one
    Rv(2,4) = vx + a*nx
    Rv(3,4) = vy + a*ny
    Rv(4,4) = H + vn*a

    diss = zero
    do i = 1, 4
      do j = 1, 4
        diss(i) = diss(i) + ws(j) * dV(j) * Rv(i,j)
      end do
    end do

    fL(1) = rhoL * vnL
    fL(2) = rhoL * vnL * vxL + pL * nx
    fL(3) = rhoL * vnL * vyL + pL * ny
    fL(4) = rhoL * vnL * HL

    fR(1) = rhoR * vnR
    fR(2) = rhoR * vnR * vxR + pR * nx
    fR(3) = rhoR * vnR * vyR + pR * ny
    fR(4) = rhoR * vnR * HR

    flux(1) = half * (fL(1) + fR(1) - diss(1))
    flux(2) = half * (fL(2) + fR(2) - diss(2))
    flux(3) = half * (fL(3) + fR(3) - diss(3))
    flux(4) = half * (fL(4) + fR(4) - diss(4))
  end subroutine Roe_prim

  subroutine ausm_plus_prim(rhoL, uL, vL, pL, aL, HL, rhoR, uR, vR, pR, aR, HR, nx, ny, flx)
    real(real64), intent(in) :: rhoL, uL, vL, pL, aL, HL
    real(real64), intent(in) :: rhoR, uR, vR, pR, aR, HR
    real(real64), intent(in) :: nx, ny
    real(real64), intent(out) :: flx(4)
    real(real64) :: af, mf, mfa, mfm, mfp, ml, mla, mlp, mr, mra, mrm
    real(real64) :: pf, ql, qr, wtl, wtr

    ql = uL*nx + vL*ny
    qr = uR*nx + vR*ny

    af = 0.5_real64 * (aL + aR)
    ml  = ql / af
    mla = abs(ml)

    mr  = qr / af
    mra = abs(mr)

    if (mla <= 1.0_real64) then
      mlp = 0.25_real64*(ml+1.0_real64)*(ml+1.0_real64) + &
            0.125_real64*(ml*ml-1.0_real64)*(ml*ml-1.0_real64)
      wtl = 0.25_real64*(ml+1.0_real64)*(ml+1.0_real64)*(2.0_real64-ml) + &
            0.1875_real64*ml*(ml*ml-1.0_real64)*(ml*ml-1.0_real64)
    else
      mlp = 0.5_real64*(ml+mla)
      wtl = 0.5_real64*(1.0_real64+ml/mla)
    end if

    if (mra <= 1.0_real64) then
      mrm = -0.25_real64*(mr-1.0_real64)*(mr-1.0_real64) - &
            0.125_real64*(mr*mr-1.0_real64)*(mr*mr-1.0_real64)
      wtr = 0.25_real64*(mr-1.0_real64)*(mr-1.0_real64)*(2.0_real64+mr) - &
            0.1875_real64*mr*(mr*mr-1.0_real64)*(mr*mr-1.0_real64)
    else
      mrm = 0.5_real64*(mr-mra)
      wtr = 0.5_real64*(1.0_real64-mr/mra)
    end if

    mf  = mlp + mrm
    mfa = abs(mf)
    mfp = 0.5_real64*(mf+mfa)
    mfm = 0.5_real64*(mf-mfa)

    pf = wtl*pL + wtr*pR

    flx(1) = af * (mfp*rhoL      + mfm*rhoR)
    flx(2) = af * (mfp*rhoL*uL   + mfm*rhoR*uR) + pf*nx
    flx(3) = af * (mfp*rhoL*vL   + mfm*rhoR*vR) + pf*ny
    flx(4) = af * (mfp*rhoL*HL   + mfm*rhoR*HR)
  end subroutine ausm_plus_prim
end module convection
