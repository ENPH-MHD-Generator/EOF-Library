//! Closed forms, consistency between functions, and the example case's tables.

use std::path::{Path, PathBuf};

use mhd_plasma::collision::{self, locate};
use mhd_plasma::constants::{E_CHARGE, KELVIN_PER_EV, K_BOLTZ};
use mhd_plasma::physics::{
    coulomb_logarithm, electron_temperature, elastic_loss_coefficients, inelastic_loss,
    plasma_transport, saha_density, saha_factor, seed_density_derivative, SeedPlasma,
};

/// Potassium as the example case configures it
fn potassium(lorentz: bool, excited_at_gas: bool) -> SeedPlasma {
    SeedPlasma {
        seed_frac: 0.01,
        chi_j: 4.341 * E_CHARGE,
        weight_ratio: 0.5,
        ion_reduced_mobility: 2.43e-4,
        lorentz,
        excited_at_gas,
    }
}

fn example_tables() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../../tests/linearHall/electron_collisions.dat")
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs().max(a.abs())
}

#[test]
fn saha_density_solves_its_quadratic() {
    let pl = potassium(false, false);
    let ns = 2.9e22;
    for t in [1500.0, 2500.0, 3500.0, 6000.0] {
        let s = saha_factor(&pl, t);
        let ne = saha_density(&pl, t, ns);
        assert!(ne > 0.0 && ne <= ns);
        assert!(close(ne * ne, s * (ns - ne), 1e-12), "T = {t}");
    }
}

#[test]
fn saha_factor_vanishes_far_below_the_ionization_energy() {
    assert_eq!(saha_factor(&potassium(false, false), 10.0), 0.0);
    assert_eq!(saha_density(&potassium(false, false), 10.0, 1e22), 0.0);
}

#[test]
fn seed_density_derivative_matches_finite_differences() {
    let pl = potassium(false, false);
    let ns = 2.9e22;
    for t in [2000.0, 2500.0, 3000.0] {
        let dt = 1e-5 * t;
        let numeric = (saha_density(&pl, t + dt, ns) - saha_density(&pl, t - dt, ns)) / (2.0 * dt);
        let analytic = seed_density_derivative(&pl, t, ns, saha_density(&pl, t, ns));
        assert!(close(analytic, numeric, 1e-7), "T = {t}: {analytic} vs {numeric}");
    }
}

#[test]
fn coulomb_logarithm_matches_the_specification_oracle() {
    // ICPSIM's Appendix A: ln Lambda = 5.62 at ne = 1e19 m^-3, theta = 0.2 eV
    let te = 0.2 * KELVIN_PER_EV;
    assert!((coulomb_logarithm(1e19, te) - 5.62).abs() < 0.005);
    assert_eq!(coulomb_logarithm(0.0, te), 0.0);
    assert_eq!(coulomb_logarithm(1e40, te), 1.0);
}

#[test]
fn locate_brackets_and_clamps() {
    let grid = [0.0, 1.0, 2.0, 4.0];
    assert_eq!(locate(&grid, -1.0), (0, 0.0));
    assert_eq!(locate(&grid, 5.0), (2, 1.0));
    assert_eq!(locate(&grid, 1.0), (1, 0.0));
    assert_eq!(locate(&grid, 3.0), (2, 0.5));
}

#[test]
fn the_example_tables_read_and_interpolate_through_their_nodes() {
    let path = example_tables();
    let cd = collision::read(path.to_str().unwrap(), false).expect("example tables");
    assert_eq!(cd.ln_t.len(), 90);
    // At a node the interpolation returns the tabulated value
    let te = cd.ln_t[40].exp() * KELVIN_PER_EV;
    assert!(close(cd.table1(&cd.ln_nu_n, te), cd.ln_nu_n[40].exp(), 1e-12));
    // The example is drifting: asking for Lorentz must fail, with the reason
    let error = collision::read(path.to_str().unwrap(), true).unwrap_err();
    assert!(error.contains("electron_transport_model: lorentz"), "{error}");
}

#[test]
fn drifting_transport_is_a_single_collision_frequency() {
    let path = example_tables();
    let cd = collision::read(path.to_str().unwrap(), false).unwrap();
    let pl = potassium(false, false);
    let (te, n, ne, b) = (2800.0, 101325.0 / (K_BOLTZ * 2500.0), 1e20, 1.5);
    let tr = plasma_transport(&cd, &pl, te, n, ne, b);
    let beta = tr.mu0 * b;
    assert!(close(tr.mu_p, tr.mu0 / (1.0 + beta * beta), 1e-14));
    assert!(close(tr.mu_h, beta * tr.mu_p, 1e-14));
    assert!(close(tr.collision_frequency * tr.mu0 * 9.1093837015e-31, E_CHARGE, 1e-14));
    let (cn, cei) = elastic_loss_coefficients(&cd, te, n, &tr);
    assert!(cn > 0.0 && cei > 0.0);
}

#[test]
fn inelastic_loss_needs_quenched_excited_states() {
    let path = example_tables();
    let cd = collision::read(path.to_str().unwrap(), false).unwrap();
    assert_eq!(inelastic_loss(&cd, &potassium(false, false), 3000.0, 2500.0, 1e20, 1e22), 0.0);
    assert!(inelastic_loss(&cd, &potassium(false, true), 3000.0, 2500.0, 1e20, 1e22) > 0.0);
}

#[test]
fn electron_temperature_balances_heating_against_loss() {
    let path = example_tables();
    let cd = collision::read(path.to_str().unwrap(), false).unwrap();
    let pl = potassium(false, false);
    let (tg, n) = (2500.0, 101325.0 / (K_BOLTZ * 2500.0));
    let ns = 0.01 * n;
    // No field: the electrons stay at the gas temperature
    assert_eq!(electron_temperature(&cd, &pl, tg, n, ns, 0.0, 0.0, 1.5, 2e4), (tg, false));
    // A field heats them, more field more
    let (te1, capped1) = electron_temperature(&cd, &pl, tg, n, ns, 1e6, 0.0, 1.5, 2e4);
    let (te2, _) = electron_temperature(&cd, &pl, tg, n, ns, 4e6, 0.0, 1.5, 2e4);
    assert!(!capped1 && tg < te1 && te1 < te2);
    // An absurd field hits the cap and says so
    assert_eq!(electron_temperature(&cd, &pl, tg, n, ns, 1e14, 0.0, 1.5, 2e4), (2e4, true));
}

#[test]
fn lorentz_cubes_are_read_in_fortran_order() {
    // theta fastest: value = 1 + it + 10 ia + 100 iw on a 2 x 2 x 2 cube
    let dir = std::env::temp_dir().join(format!("mhd_plasma_cube_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("cube.dat");
    let cube: Vec<String> = (0..8)
        .map(|k| format!("{}", 1 + (k % 2) + 10 * ((k / 2) % 2) + 100 * (k / 4)))
        .collect();
    let column = "1 1";
    let text = format!(
        "! test\nformat 1\n1\ntheta 2\n0.1 1.0\nmomentum_frequency_n 2\n{column}\n\
         elastic_recoil_n 2\n{column}\nei_recoil_per_coulomb 2\n{column}\n\
         k_excitation_K 2\n{column}\nk_superelastic_K 2\n{column}\n\
         coulomb_parameter 2\n1e-3 1e-1\nmagnetic_parameter 2\n1e-2 1e0\n\
         mobility_n 8\n{c}\npedersen_mobility_n 8\n{c}\nhall_mobility_n 8\n{c}\nend 0\n",
        c = cube.join(" ")
    );
    std::fs::write(&path, text).unwrap();
    let cd = collision::read(path.to_str().unwrap(), true).unwrap();
    let te_hi = 1.0 * KELVIN_PER_EV;
    // (it, ia, iw) = (1, 0, 1) at the upper theta, lowest A, highest W
    assert!(close(cd.table3(&cd.ln_mu, te_hi, 1e-3, 1.0), 102.0, 1e-12));
    std::fs::remove_dir_all(&dir).unwrap();
}
