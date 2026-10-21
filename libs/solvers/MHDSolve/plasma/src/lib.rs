//! Plasma physics for EOF-Library's MHD solvers, behind a C ABI.
//!
//! Elmer's MHDSolve calls these through the Fortran interface module
//! `MHDPlasma.F90`; a C or C++ caller (an OpenFOAM solver) uses
//! `include/mhd_plasma.h`. The library holds no Elmer or OpenFOAM types: it
//! evaluates local physics from the case's exported data
//! (`electron_collisions.dat`) and the state handed to it.
//!
//! The collision tables are read once per process and kept in a process-wide
//! cache, as the Fortran module kept them. Every MPI rank is its own process.

pub mod collision;
pub mod constants;
pub mod physics;

use std::ffi::{c_char, CStr};
use std::sync::RwLock;

use collision::CollisionData;
pub use physics::{PlasmaTransport, SeedPlasma};

static TABLES: RwLock<Option<CollisionData>> = RwLock::new(None);

fn with_tables<T>(f: impl FnOnce(&CollisionData) -> T) -> T {
    let guard = TABLES.read().unwrap_or_else(|e| e.into_inner());
    let cd = guard
        .as_ref()
        .expect("mhd_plasma: collision data used before mhd_plasma_read_collision_data");
    f(cd)
}

/// Read the collision tables once; a later call with the same file does nothing.
///
/// Returns 0 on success. On failure returns 1 and writes a NUL-terminated
/// message of at most `message_len` bytes (including the NUL) to `message`.
///
/// # Safety
/// `file_name` must be a NUL-terminated string, and `message` must point to
/// `message_len` writable bytes.
#[no_mangle]
pub unsafe extern "C" fn mhd_plasma_read_collision_data(
    file_name: *const c_char,
    lorentz: bool,
    message: *mut c_char,
    message_len: usize,
) -> i32 {
    let name = CStr::from_ptr(file_name).to_string_lossy().into_owned();
    {
        let current = TABLES.read().unwrap_or_else(|e| e.into_inner());
        if current.as_ref().is_some_and(|cd| cd.file_name == name) {
            return 0;
        }
    }
    match collision::read(&name, lorentz) {
        Ok(cd) => {
            *TABLES.write().unwrap_or_else(|e| e.into_inner()) = Some(cd);
            0
        }
        Err(text) => {
            *TABLES.write().unwrap_or_else(|e| e.into_inner()) = None;
            if message_len > 0 {
                let bytes = text.as_bytes();
                let n = bytes.len().min(message_len - 1);
                std::ptr::copy_nonoverlapping(bytes.as_ptr().cast::<c_char>(), message, n);
                *message.add(n) = 0;
            }
            1
        }
    }
}

/// Saha electron density, m^-3, at electron temperature `t` [K] and seed density `ns`.
#[no_mangle]
pub extern "C" fn mhd_plasma_saha_density(pl: &SeedPlasma, t: f64, ns: f64) -> f64 {
    physics::saha_density(pl, t, ns)
}

/// d(ne)/dT of the Saha electron density at fixed seed density, m^-3/K.
#[no_mangle]
pub extern "C" fn mhd_plasma_seed_density_derivative(
    pl: &SeedPlasma,
    t: f64,
    ns: f64,
    ne: f64,
) -> f64 {
    physics::seed_density_derivative(pl, t, ns, ne)
}

/// Electron and ion mobilities and the Coulomb logarithm at one point.
#[no_mangle]
pub extern "C" fn mhd_plasma_transport(
    pl: &SeedPlasma,
    te: f64,
    n: f64,
    ne: f64,
    bmag: f64,
    tr: &mut PlasmaTransport,
) {
    *tr = with_tables(|cd| physics::plasma_transport(cd, pl, te, n, ne, bmag));
}

/// Elastic loss coefficients Cn and Cei, W/(m^3 K), for the transport state `tr`.
#[no_mangle]
pub extern "C" fn mhd_plasma_elastic_loss_coefficients(
    te: f64,
    n: f64,
    tr: &PlasmaTransport,
    cn: &mut f64,
    cei: &mut f64,
) {
    (*cn, *cei) = with_tables(|cd| physics::elastic_loss_coefficients(cd, te, n, tr));
}

/// Net K 4s-4p inelastic energy loss of the electrons, W/m^3.
#[no_mangle]
pub extern "C" fn mhd_plasma_inelastic_loss(
    pl: &SeedPlasma,
    te: f64,
    tg: f64,
    ne: f64,
    nk: f64,
) -> f64 {
    with_tables(|cd| physics::inelastic_loss(cd, pl, te, tg, ne, nk))
}

/// Local electron temperature, K; `at_max` is set when it is capped at `te_max`.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub extern "C" fn mhd_plasma_electron_temperature(
    pl: &SeedPlasma,
    tg: f64,
    n: f64,
    ns: f64,
    e_par2: f64,
    e_perp2: f64,
    bmag: f64,
    te_max: f64,
    at_max: &mut bool,
) -> f64 {
    let (te, capped) = with_tables(|cd| {
        physics::electron_temperature(cd, pl, tg, n, ns, e_par2, e_perp2, bmag, te_max)
    });
    *at_max = capped;
    te
}
