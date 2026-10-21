//! Physical constants, CODATA 2018, as MHDSolve has always used them.

/// Boltzmann constant, J/K
pub const K_BOLTZ: f64 = 1.380649e-23;
/// Electron mass, kg
pub const E_MASS: f64 = 9.1093837015e-31;
/// Elementary charge, C
pub const E_CHARGE: f64 = 1.602176634e-19;
/// Planck constant, J s
pub const H_PLANCK: f64 = 6.62607015e-34;
/// The same double the Fortran literal 3.14159265358979323846d0 rounds to
pub const PI: f64 = std::f64::consts::PI;
/// 11604.5 K per eV
pub const KELVIN_PER_EV: f64 = E_CHARGE / K_BOLTZ;
/// Loschmidt number, m^-3 (273.15 K, 101325 Pa): reference for reduced mobility
pub const LOSCHMIDT: f64 = 2.6867811e25;
/// NRL electron-ion collision frequency coefficient, m^3 eV^(3/2)/s
pub const NRL_COEFFICIENT: f64 = 2.91e-12;
