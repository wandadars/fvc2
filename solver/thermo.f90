!------------------------------------------------------------------------------
! Ideal-gas thermodynamic helper functions.
!------------------------------------------------------------------------------
module thermo
  use iso_fortran_env, only: real64
  use data, only: gamma_gas, gas_constant
  implicit none
contains
  ! Total enthalpy H_t = (rho*E + p) / rho
  real(real64) function MixtPerf_Ht(density, rho_u, rho_v, rho_E, pressure)
    real(real64), intent(in) :: density, rho_u, rho_v, rho_E, pressure
    MixtPerf_Ht = (rho_E + pressure) / density
  end function MixtPerf_Ht

  ! Temperature from ideal gas: T = p / (R * rho)
  real(real64) function MixtPerf_T(density, rho_u, rho_v, rho_E, pressure)
    real(real64), intent(in) :: density, rho_u, rho_v, rho_E, pressure
    MixtPerf_T = pressure / gas_constant / density
  end function MixtPerf_T

  ! Pressure from conservative variables (2D Euler):
  ! p = (gamma - 1) * (rho*E - 0.5*(rho*state^2 + rho*v^2)/rho)
  real(real64) function MixtPerf_P(density, rho_u, rho_v, rho_E)
    real(real64), intent(in) :: density, rho_u, rho_v, rho_E
    MixtPerf_P = (gamma_gas - 1.0_real64) * (rho_E - 0.5_real64 * (rho_u*rho_u + rho_v*rho_v) / density)
  end function MixtPerf_P

  ! Total specific energy e_t = (p/(gamma-1) + 0.5*rho*(state^2+v^2)) / rho
  real(real64) function MixtPerf_et(density, rho_u, rho_v, rho_E, pressure)
    real(real64), intent(in) :: density, rho_u, rho_v, rho_E, pressure
    MixtPerf_et = (pressure/(gamma_gas-1.0_real64) + 0.5_real64*(rho_u*rho_u + rho_v*rho_v)/density) / density
  end function MixtPerf_et

  ! Speed of sound a = sqrt(gamma * R * T)
  real(real64) function MixtPerf_A(density, rho_u, rho_v, rho_E, temperature)
    real(real64), intent(in) :: density, rho_u, rho_v, rho_E, temperature
    MixtPerf_A = sqrt(gamma_gas * gas_constant * temperature)
  end function MixtPerf_A
end module thermo
