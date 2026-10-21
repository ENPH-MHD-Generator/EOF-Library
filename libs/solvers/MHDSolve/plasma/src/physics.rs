//! Seeded-plasma physics: argon with potassium seed.
//!
//! Only the seed ionizes (two-temperature Saha at the electron temperature);
//! heavy particles are at the gas temperature. Electron collision data come from
//! case preparation (`electron_collisions.dat`, a Maxwellian EEDF at Te):
//! electron-neutral momentum transfer, recoil energy losses and K(4p)
//! excitation. Electron-ion (Coulomb) collisions use the NRL Coulomb logarithm.
//! Ions carry current with a polarization mobility (ion slip).
//!
//! A line-for-line port of the Fortran module this replaced. The operations
//! are kept in the same order, so the results agree to the last bit wherever
//! the platform's libm is the same.

use crate::collision::CollisionData;
use crate::constants::{
    E_CHARGE, E_MASS, H_PLANCK, KELVIN_PER_EV, K_BOLTZ, LOSCHMIDT, NRL_COEFFICIENT, PI,
};

/// The seed and the model options, filled from the case's Material section.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct SeedPlasma {
    /// Seed atoms per heavy particle
    pub seed_frac: f64,
    /// Seed ionization energy, J
    pub chi_j: f64,
    /// g_ion / g_neutral of the seed
    pub weight_ratio: f64,
    /// Ion mobility, m^2/(V s), at the Loschmidt density
    pub ion_reduced_mobility: f64,
    /// Lorentz electron transport (else drifting)
    pub lorentz: bool,
    /// K(4p) population at the gas temperature
    pub excited_at_gas: bool,
}

/// Local transport state at one point.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PlasmaTransport {
    /// Electron density, m^-3
    pub ne: f64,
    /// Coulomb logarithm
    pub ln_lambda: f64,
    /// Electron mobilities, parallel (DC), Pedersen and Hall, m^2/(V s)
    pub mu0: f64,
    pub mu_p: f64,
    pub mu_h: f64,
    /// Ion mobilities (Hall opposite)
    pub mu_i: f64,
    pub mu_ip: f64,
    pub mu_ih: f64,
    /// e / (me mu0), 1/s
    pub collision_frequency: f64,
}

/// Right-hand side S(T) of the Saha equation ne^2 / (nS - ne) = S, m^-3
pub fn saha_factor(pl: &SeedPlasma, t: f64) -> f64 {
    let expo = pl.chi_j / (K_BOLTZ * t);
    if expo > 700.0 {
        0.0
    } else {
        2.0 * pl.weight_ratio
            * ((2.0 * PI * E_MASS * K_BOLTZ * t) / (H_PLANCK * H_PLANCK)).powf(1.5)
            * (-expo).exp()
    }
}

/// Saha electron density at electron temperature `t` and seed density `ns`
pub fn saha_density(pl: &SeedPlasma, t: f64, ns: f64) -> f64 {
    // Positive root of ne^2 + S ne - S nS = 0, in a form that avoids
    // cancellation when S >> nS
    let s = saha_factor(pl, t);
    if s > 0.0 {
        2.0 * s * ns / (s + (s * s + 4.0 * s * ns).sqrt())
    } else {
        0.0
    }
}

/// d(ne)/dT of the Saha electron density `ne` at electron temperature `t`
/// and fixed seed density `ns`. From ne^2 + S ne - S nS = 0,
///   d(ne)/dS = (nS - ne) / (2 ne + S),  dS/dT = S (3/(2T) + chi/(kB T^2))
pub fn seed_density_derivative(pl: &SeedPlasma, t: f64, ns: f64, ne: f64) -> f64 {
    let s = saha_factor(pl, t);
    if s <= 0.0 || 2.0 * ne + s <= 0.0 {
        return 0.0;
    }
    s * (1.5 / t + pl.chi_j / (K_BOLTZ * t * t)) * (ns - ne).max(0.0) / (2.0 * ne + s)
}

/// NRL electron-ion Coulomb logarithm (Z = 1, Te below 10 eV), floored at 1
pub fn coulomb_logarithm(ne: f64, te: f64) -> f64 {
    let theta = te / KELVIN_PER_EV;
    if ne <= 0.0 {
        return 0.0;
    }
    (23.0 - ((ne * 1.0e-6).sqrt() * theta.powf(-1.5)).ln()).max(1.0)
}

/// Electron and ion mobilities at electron temperature `te`, heavy-particle
/// density `n`, electron density `ne` and |B|. Parallel (DC), Pedersen and
/// Hall mobilities; for ions the Hall current runs the opposite way.
///
/// Drifting: nu = N nu_en/N(theta) + 2.91e-12 ne lnL theta^-3/2 (NRL), and the
///   tensor of a single collision frequency, beta = mu |B|.
/// Lorentz: collision frequencies add inside the velocity average, so the
///   Pedersen and Hall mobilities are tabulated on (theta, ne lnL / N,
///   omega_ce / N) and are not those of any single frequency.
/// Ions: mu_i = mu_i0 N0 / N (polarization collisions, independent of speed).
pub fn plasma_transport(
    cd: &CollisionData,
    pl: &SeedPlasma,
    te: f64,
    n: f64,
    ne: f64,
    bmag: f64,
) -> PlasmaTransport {
    let theta = te / KELVIN_PER_EV;
    let mut tr = PlasmaTransport { ne, ln_lambda: coulomb_logarithm(ne, te), ..Default::default() };
    if pl.lorentz {
        let a = ne * tr.ln_lambda / n;
        let w = E_CHARGE * bmag / (E_MASS * n);
        tr.mu0 = cd.table3(&cd.ln_mu, te, a, w) / n;
        tr.mu_p = cd.table3(&cd.ln_mu_p, te, a, w) / n;
        tr.mu_h = cd.table3(&cd.ln_mu_h, te, a, w) / n;
        if bmag <= 0.0 {
            tr.mu_p = tr.mu0;
            tr.mu_h = 0.0;
        }
    } else {
        let nu = n * cd.table1(&cd.ln_nu_n, te)
            + NRL_COEFFICIENT * ne * tr.ln_lambda * theta.powf(-1.5);
        tr.mu0 = E_CHARGE / (E_MASS * nu);
        let beta = tr.mu0 * bmag;
        tr.mu_p = tr.mu0 / (1.0 + beta * beta);
        tr.mu_h = tr.mu0 * beta / (1.0 + beta * beta);
    }
    tr.collision_frequency = E_CHARGE / (E_MASS * tr.mu0);

    tr.mu_i = pl.ion_reduced_mobility * LOSCHMIDT / n;
    let beta_i = tr.mu_i * bmag;
    tr.mu_ip = tr.mu_i / (1.0 + beta_i * beta_i);
    tr.mu_ih = tr.mu_i * beta_i / (1.0 + beta_i * beta_i);
    tr
}

/// Elastic (recoil) energy loss coefficients of the electrons, W/(m^3 K), so
/// that the loss is (Cn + Cei) (Te - Tg): to neutrals, e ne N G(theta) / Te,
/// and to ions, e ne^2 lnL Gei(theta) / Te (the tables hold the
/// electron-temperature integrals; (1 - Tg/Te) = (Te - Tg)/Te).
pub fn elastic_loss_coefficients(
    cd: &CollisionData,
    te: f64,
    n: f64,
    tr: &PlasmaTransport,
) -> (f64, f64) {
    let cn = E_CHARGE * tr.ne * n * cd.table1(&cd.ln_elastic, te) / te;
    let cei = E_CHARGE * tr.ne * tr.ne * tr.ln_lambda * cd.table1(&cd.ln_ei, te) / te;
    (cn, cei)
}

/// Net inelastic (K 4s-4p) energy loss of the electrons, W/m^3, for K(4p)
/// populations at the gas temperature `tg` (`excited_at_gas`). Zero when the
/// excited states follow the electrons: excitation and de-excitation then
/// balance by detailed balance.
pub fn inelastic_loss(
    cd: &CollisionData,
    pl: &SeedPlasma,
    te: f64,
    tg: f64,
    ne: f64,
    nk: f64,
) -> f64 {
    if !pl.excited_at_gas || ne <= 0.0 || nk <= 0.0 {
        return 0.0;
    }
    let r = cd.exc_weight * (-(cd.exc_energy * KELVIN_PER_EV / tg).min(700.0)).exp();
    let y_low = 1.0 / (1.0 + r);
    let y_up = r / (1.0 + r);
    E_CHARGE * ne * nk * cd.exc_energy
        * (y_low * cd.table1(&cd.ln_kexc, te) - y_up * cd.table1(&cd.ln_ksup, te))
}

/// Electron temperature from the local balance of Joule heating of the
/// electrons, e ne (mu_0 E'_par^2 + mu_P E'_perp^2), against elastic and
/// inelastic losses, bracketed and bisected on [Tg, TeMax]. ne follows Te
/// through the Saha equation. The flag is set when the balance lies above
/// TeMax and Te is capped there.
#[allow(clippy::too_many_arguments)]
pub fn electron_temperature(
    cd: &CollisionData,
    pl: &SeedPlasma,
    tg: f64,
    n: f64,
    ns: f64,
    e_par2: f64,
    e_perp2: f64,
    bmag: f64,
    te_max: f64,
) -> (f64, bool) {
    let mut te = tg;
    if e_par2 + e_perp2 <= 0.0 || te_max <= tg {
        return (te, false);
    }
    let residual = |t: f64| -> f64 {
        let ne = saha_density(pl, t, ns);
        // A tiny floor keeps the per-electron balance defined in cold gas
        let ne_eff = ne.max(1.0e-12 * ns);
        let tr = plasma_transport(cd, pl, t, n, ne_eff, bmag);
        let (cn, cei) = elastic_loss_coefficients(cd, t, n, &tr);
        E_CHARGE * (tr.mu0 * e_par2 + tr.mu_p * e_perp2)
            - ((cn + cei) * (t - tg) + inelastic_loss(cd, pl, t, tg, ne_eff, (ns - ne).max(0.0)))
                / ne_eff
    };

    // Residual (heating - loss per electron) is positive at Tg; grow the upper
    // bracket until it turns negative or reaches the cap
    let mut te_lo = tg;
    let mut te_hi = (1.1 * tg).min(te_max);
    for _ in 0..60 {
        if residual(te_hi) <= 0.0 {
            break;
        }
        if te_hi >= te_max {
            return (te_max, true);
        }
        te_lo = te_hi;
        te_hi = (tg + 2.0 * (te_hi - tg)).min(te_max);
    }

    for _ in 0..100 {
        te = 0.5 * (te_lo + te_hi);
        if te_hi - te_lo < 1.0e-3 {
            break;
        }
        if residual(te) > 0.0 {
            te_lo = te;
        } else {
            te_hi = te;
        }
    }
    (te, false)
}
