!------------------------------------------------------------------------------
! Convective fluxes (AUSM+ and Roe) for the Euler equations.
!------------------------------------------------------------------------------
module convection
  use iso_fortran_env, only: real64
  use data
  use thermo, only: MixtPerf_P, MixtPerf_T, MixtPerf_A, MixtPerf_Ht
  implicit none
contains
  subroutine convective_flux_interior(bdum)
    real(real64), intent(inout) :: bdum(num_state_entries)
    real(real64) :: flux(4), ul1(4), ur1(4), nhatx, nhaty

    ul1 = state(elem_left, :)
    ur1 = state(elem_right, :)

    nhatx = element_geom(4+side_index_left*2, elem_left)
    nhaty = element_geom(5+side_index_left*2, elem_left)

    call convective_flux(ul1, ur1, nhatx, nhaty, flux)

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
    real(real64) :: flux(4), ul1(4), ur1(4), rnhatx, rnhaty
    real(real64) :: rdens, rpres, ru, rv, ret
    real(real64) :: uvel, vvel, vn
    integer :: idum, ibc

    idum = edge_connectivity(5, edge_id)
    if (idum <= 0) return
    ibc = bc_kind(idum)

    if (ibc == 3) then
      ul1 = state(elem_left, :)

      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      ur1 = state(elem_left, :)

      call convective_flux(ul1, ur1, rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    elseif (ibc == 2) then
      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      rdens = bc_values(1, idum)
      rpres = bc_values(2, idum)
      ru    = bc_values(3, idum)
      rv    = bc_values(4, idum)
      if (rdens <= 0.0_real64 .or. rpres <= 0.0_real64) then
        rdens = ic_uniform_state(1)
        rpres = ic_uniform_state(2)
        ru    = ic_uniform_state(3)
        rv    = ic_uniform_state(4)
      end if
      ret = (rpres/(gamma_gas-1.0_real64) + 0.5_real64*rdens*(ru*ru + rv*rv)) / rdens

      flux(1) = rdens * (ru*rnhatx + rv*rnhaty)
      flux(2) = rdens * ru * (ru*rnhatx + rv*rnhaty) + rpres * rnhatx
      flux(3) = rdens * rv * (ru*rnhatx + rv*rnhaty) + rpres * rnhaty
      flux(4) = (ret + rpres/rdens) * rdens * (ru*rnhatx + rv*rnhaty)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    elseif (ibc == 4) then
      ! Subsonic farfield: inflow uses specified freestream, outflow extrapolates interior.
      ul1 = state(elem_left, :)

      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      rdens = ul1(1)
      if (rdens > 0.0_real64) then
        uvel = ul1(2) / rdens
        vvel = ul1(3) / rdens
      else
        uvel = 0.0_real64
        vvel = 0.0_real64
      end if
      vn = uvel * rnhatx + vvel * rnhaty

      if (vn >= 0.0_real64) then
        ur1 = ul1
      else
        rdens = bc_values(1, idum)
        rpres = bc_values(2, idum)
        ru    = bc_values(3, idum)
        rv    = bc_values(4, idum)
        if (rdens <= 0.0_real64 .or. rpres <= 0.0_real64) then
          rdens = ic_uniform_state(1)
          rpres = ic_uniform_state(2)
          ru    = ic_uniform_state(3)
          rv    = ic_uniform_state(4)
        end if
        ret = (rpres/(gamma_gas-1.0_real64) + 0.5_real64*rdens*(ru*ru + rv*rv)) / rdens
        ur1(1) = rdens
        ur1(2) = rdens * ru
        ur1(3) = rdens * rv
        ur1(4) = rdens * ret
      end if

      call convective_flux(ul1, ur1, rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length

    elseif (ibc == 1) then
      ul1 = state(elem_left, :)

      rnhatx = element_geom(4+side_index_left*2, elem_left)
      rnhaty = element_geom(5+side_index_left*2, elem_left)

      ur1(1) = state(elem_left, 1)
      if (abs(rnhaty) < 1.0e-9_real64) then
        ur1(2) = 0.0_real64
        ur1(3) = state(elem_left, 3)
      elseif (abs(rnhatx) < 1.0e-9_real64) then
        ur1(2) = state(elem_left, 2)
        ur1(3) = 0.0_real64
      else
        ur1(2) = state(elem_left, 2) * 0.0_real64
        ur1(3) = -ur1(2) * rnhatx / rnhaty
      end if
      ur1(4) = state(elem_left, 4)

      call convective_flux(ul1, ur1, rnhatx, rnhaty, flux)

      bdum(elem_left_rho_idx)  = bdum(elem_left_rho_idx)  - flux(1) * edge_length
      bdum(elem_left_rhou_idx) = bdum(elem_left_rhou_idx) - flux(2) * edge_length
      bdum(elem_left_rhov_idx) = bdum(elem_left_rhov_idx) - flux(3) * edge_length
      bdum(elem_left_rhoe_idx) = bdum(elem_left_rhoe_idx) - flux(4) * edge_length
    end if
  end subroutine convective_flux_boundary

  subroutine convective_flux(ul1, ur1, nhatx, nhaty, flux)
    real(real64), intent(in) :: ul1(4), ur1(4), nhatx, nhaty
    real(real64), intent(out) :: flux(4)

    select case (convection_scheme)
    case (0)
      flux = 0.0_real64
    case (1)
      call ausm_plus(ul1, ur1, nhatx, nhaty, flux)
    case (2)
      call Roe(ul1, ur1, nhatx, nhaty, flux)
    end select
  end subroutine convective_flux

  subroutine Roe(uL, uR, nx, ny, flux)
    real(real64), intent(in) :: uL(4), uR(4)
    real(real64), intent(in) :: nx, ny
    real(real64), intent(out) :: flux(4)

    real(real64) :: gamma
    real(real64) :: zero, fifth, half, one, two
    real(real64) :: tx, ty
    real(real64) :: vxL, vxR, vyL, vyR
    real(real64) :: rhoL, rhoR, pL, pR
    real(real64) :: vnL, vnR, vtL, vtR
    real(real64) :: aL, aR, HL, HR, tL, tR
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

    rhoL = uL(1)
    vxL = uL(2) / uL(1)
    vyL = uL(3) / uL(1)
    vnL = vxL*nx + vyL*ny
    vtL = vxL*tx + vyL*ty

    pL = MixtPerf_P(uL(1), uL(2), uL(3), uL(4))
    tL = MixtPerf_T(uL(1), uL(2), uL(3), uL(4), pL)
    aL = MixtPerf_A(uL(1), uL(2), uL(3), uL(4), tL)
    hL = MixtPerf_Ht(uL(1), uL(2), uL(3), uL(4), pL)

    rhoR = uR(1)
    vxR = uR(2) / uR(1)
    vyR = uR(3) / uR(1)
    vnR = vxR*nx + vyR*ny
    vtR = vxR*tx + vyR*ty

    pR = MixtPerf_P(uR(1), uR(2), uR(3), uR(4))
    tR = MixtPerf_T(uR(1), uR(2), uR(3), uR(4), pR)
    aR = MixtPerf_A(uR(1), uR(2), uR(3), uR(4), tR)
    hR = MixtPerf_Ht(uR(1), uR(2), uR(3), uR(4), pR)

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
  end subroutine Roe

  subroutine ausm_plus(rul, rur, nx, ny, flx)
    real(real64), intent(in) :: rul(4), rur(4)
    real(real64), intent(in) :: nx, ny
    real(real64), intent(out) :: flx(4)

    real(real64) :: al, ar, pl, pr, rl, rr, ul, ur, vl, vr
    real(real64) :: tl, tr
    real(real64) :: af, mf, mfa, mfm, mfp, ml, mla, mlp, mr, mra, mrm
    real(real64) :: pf, ql, qr, wtl, wtr, Hl, Hr
    real(real64) :: retl, retr

    rl   = rul(1)
    ul   = rul(2) / rul(1)
    vl   = rul(3) / rul(1)
    retl = rul(4) / rul(1)

    rr   = rur(1)
    ur   = rur(2) / rur(1)
    vr   = rur(3) / rur(1)
    retr = rur(4) / rur(1)

    pl = MixtPerf_P(rul(1), rul(2), rul(3), rul(4))
    pr = MixtPerf_P(rur(1), rur(2), rur(3), rur(4))

    tl = MixtPerf_T(rul(1), rul(2), rul(3), rul(4), pl)
    tr = MixtPerf_T(rur(1), rur(2), rur(3), rur(4), pr)

    al = MixtPerf_A(rul(1), rul(2), rul(3), rul(4), tl)
    ar = MixtPerf_A(rur(1), rur(2), rur(3), rur(4), tr)

    hl = MixtPerf_Ht(rul(1), rul(2), rul(3), rul(4), pl)
    hr = MixtPerf_Ht(rur(1), rur(2), rur(3), rur(4), pr)

    ql = ul*nx + vl*ny
    qr = ur*nx + vr*ny

    af = 0.5_real64 * (al + ar)
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

    pf = wtl*pl + wtr*pr

    flx(1) = af * (mfp*rl      + mfm*rr)
    flx(2) = af * (mfp*rl*ul   + mfm*rr*ur) + pf*nx
    flx(3) = af * (mfp*rl*vl   + mfm*rr*vr) + pf*ny
    flx(4) = af * (mfp*rl*Hl   + mfm*rr*Hr)
  end subroutine ausm_plus
end module convection
