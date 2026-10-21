//! Electron collision tables: the reader for `electron_collisions.dat` and the
//! log-log interpolation MHDSolve evaluates them with.
//!
//! The file is written by case preparation (channelsim, from plasmakit) and is
//! the contract between the two: blocks of `<name> <count>` followed by
//! `count` values, `!` lines are comments, `end 0` closes it, and the first
//! block is `format 1`. Unknown blocks are skipped, so the writer can add data
//! for newer readers without breaking this one.
//!
//! Why interpolate in log-log space: these rates and mobilities span many
//! decades across the tables, and what matters to the conductivity and the
//! electron energy balance is their relative accuracy. In the cold wall layer
//! the values are orders of magnitude below those in the core, and
//! interpolating the values directly would let the larger neighbour dominate,
//! overestimating them there by far more than the table spacing suggests.
//! Between grid points the physics is close to a power law (collision rates
//! ~ theta^a; magnetized Pedersen and Hall mobilities ~ B^-2 and B^-1), which
//! log-log interpolation reproduces exactly, so the error is small and the
//! same at every magnitude. It also keeps every result positive. The worst
//! case, the Arrhenius-like K(4p) excitation rate, is ~0.7% off at 0.1 eV.
//! Outside a grid the end value is held.

use std::fs;

use crate::constants::KELVIN_PER_EV;

/// The tables of one file, as natural logarithms of the values.
#[derive(Debug, Clone, Default)]
pub struct CollisionData {
    pub file_name: String,
    /// ln theta, eV
    pub ln_t: Vec<f64>,
    /// ln of the Coulomb parameter ne lnL / N (Lorentz only)
    pub ln_a: Vec<f64>,
    /// ln of the magnetic parameter omega_ce / N (Lorentz only)
    pub ln_w: Vec<f64>,
    pub ln_nu_n: Vec<f64>,
    pub ln_elastic: Vec<f64>,
    pub ln_ei: Vec<f64>,
    pub ln_kexc: Vec<f64>,
    pub ln_ksup: Vec<f64>,
    /// Cubes in Fortran order, theta fastest: index it + nT (ia + nA iw)
    pub ln_mu: Vec<f64>,
    pub ln_mu_p: Vec<f64>,
    pub ln_mu_h: Vec<f64>,
    pub exc_energy: f64,
    pub exc_weight: f64,
}

fn safe_log(values: &[f64]) -> Vec<f64> {
    values.iter().map(|v| v.max(1.0e-300).ln()).collect()
}

/// Read and check a collision data file.
///
/// `lorentz` asks for the mobility cubes too. The error messages are the ones
/// MHDSolve reported when the reader was Fortran.
pub fn read(file_name: &str, lorentz: bool) -> Result<CollisionData, String> {
    let text = fs::read_to_string(file_name)
        .map_err(|_| format!("Cannot open electron collision data {file_name}"))?;
    let mut cd = CollisionData { exc_weight: 1.0, ..Default::default() };
    let (mut n_a, mut n_w) = (0usize, 0usize);
    let mut cubes: Vec<(&str, Vec<f64>)> = Vec::new();

    let mut lines = text.lines();
    while let Some(raw) = lines.next() {
        let line = raw.trim_start();
        if line.trim_end().is_empty() || line.starts_with('!') {
            continue;
        }
        let mut header = line.split_whitespace();
        let (name, count) = match (header.next(), header.next().map(str::parse::<usize>)) {
            (Some(name), Some(Ok(count))) => (name, count),
            _ => {
                return Err(format!(
                    "Malformed block header in {file_name}: {}",
                    line.trim_end()
                ))
            }
        };
        if name == "end" {
            break;
        }
        // The values may span lines; a line's leftovers after the last value
        // are discarded, as a Fortran list-directed read does
        let mut values = Vec::with_capacity(count);
        while values.len() < count {
            let Some(next) = lines.next() else { break };
            for token in next.split_whitespace() {
                if values.len() == count {
                    break;
                }
                match token.parse::<f64>() {
                    Ok(v) => values.push(v),
                    Err(_) => return Err(format!("Short block {name} in {file_name}")),
                }
            }
        }
        if values.len() < count {
            return Err(format!("Short block {name} in {file_name}"));
        }
        let scalar = values.first().copied().unwrap_or(0.0);
        match name {
            "format" => {
                if scalar.round() != 1.0 {
                    return Err(format!("Unsupported collision data format in {file_name}"));
                }
            }
            "theta" => cd.ln_t = values.iter().map(|v| v.ln()).collect(),
            "momentum_frequency_n" => cd.ln_nu_n = safe_log(&values),
            "elastic_recoil_n" => cd.ln_elastic = safe_log(&values),
            "ei_recoil_per_coulomb" => cd.ln_ei = safe_log(&values),
            "k_excitation_K" => cd.ln_kexc = safe_log(&values),
            "k_superelastic_K" => cd.ln_ksup = safe_log(&values),
            "excitation_energy_K" => cd.exc_energy = scalar,
            "excitation_weight_ratio_K" => cd.exc_weight = scalar,
            "coulomb_parameter" => {
                n_a = count;
                cd.ln_a = values.iter().map(|v| v.ln()).collect();
            }
            "magnetic_parameter" => {
                n_w = count;
                cd.ln_w = values.iter().map(|v| v.ln()).collect();
            }
            "mobility_n" | "pedersen_mobility_n" | "hall_mobility_n" => {
                // Checked against the grids once every block has been read
                let key = match name {
                    "mobility_n" => "mobility_n",
                    "pedersen_mobility_n" => "pedersen_mobility_n",
                    _ => "hall_mobility_n",
                };
                cubes.push((key, safe_log(&values)));
            }
            _ => {}
        }
    }

    for (key, cube) in cubes {
        if cube.len() != cd.ln_t.len() * n_a * n_w {
            return Err(format!(
                "Mobility table size does not match its grids in {file_name}"
            ));
        }
        match key {
            "mobility_n" => cd.ln_mu = cube,
            "pedersen_mobility_n" => cd.ln_mu_p = cube,
            _ => cd.ln_mu_h = cube,
        }
    }

    let complete = cd.ln_t.len() >= 2
        && [&cd.ln_nu_n, &cd.ln_elastic, &cd.ln_ei, &cd.ln_kexc, &cd.ln_ksup]
            .iter()
            .all(|column| !column.is_empty());
    if !complete {
        return Err(format!("Incomplete electron collision data in {file_name}"));
    }
    if lorentz && (cd.ln_mu.is_empty() || cd.ln_mu_p.is_empty() || cd.ln_mu_h.is_empty()) {
        return Err(format!(
            "Lorentz transport needs the mobility tables; regenerate {file_name} \
             with electron_transport_model: lorentz"
        ));
    }
    cd.file_name = file_name.to_owned();
    Ok(cd)
}

/// Bracketing index (0-based) and weight of x on an increasing grid, clamped
/// to its ends.
pub fn locate(grid: &[f64], x: f64) -> (usize, f64) {
    let n = grid.len();
    if x <= grid[0] {
        return (0, 0.0);
    } else if x >= grid[n - 1] {
        return (n - 2, 1.0);
    }
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi + 2) / 2 - 1; // the Fortran midpoint, in 1-based terms
        if grid[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo, (x - grid[lo]) / (grid[lo + 1] - grid[lo]))
}

impl CollisionData {
    fn n_t(&self) -> usize {
        self.ln_t.len()
    }

    fn n_a(&self) -> usize {
        self.ln_a.len()
    }

    /// exp of a log-tabulated column at electron temperature `te` [K] (log-log)
    pub fn table1(&self, ln_values: &[f64], te: f64) -> f64 {
        let (i, w) = locate(&self.ln_t, (te / KELVIN_PER_EV).ln());
        ((1.0 - w) * ln_values[i] + w * ln_values[i + 1]).exp()
    }

    /// exp of a log-tabulated cube at (theta, Coulomb parameter, magnetic parameter)
    pub fn table3(&self, ln_values: &[f64], te: f64, a: f64, w: f64) -> f64 {
        let (it, wt) = locate(&self.ln_t, (te / KELVIN_PER_EV).ln());
        let (ia, wa) = locate(&self.ln_a, a.max(1.0e-300).ln());
        let (iw, ww) = locate(&self.ln_w, w.max(1.0e-300).ln());
        let (nt, na) = (self.n_t(), self.n_a());
        let v = |t: usize, a: usize, w: usize| ln_values[t + nt * (a + na * w)];
        ((1.0 - ww)
            * ((1.0 - wa) * ((1.0 - wt) * v(it, ia, iw) + wt * v(it + 1, ia, iw))
                + wa * ((1.0 - wt) * v(it, ia + 1, iw) + wt * v(it + 1, ia + 1, iw)))
            + ww * ((1.0 - wa) * ((1.0 - wt) * v(it, ia, iw + 1) + wt * v(it + 1, ia, iw + 1))
                + wa * ((1.0 - wt) * v(it, ia + 1, iw + 1) + wt * v(it + 1, ia + 1, iw + 1))))
            .exp()
    }
}
