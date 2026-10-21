#ifndef MHD_PLASMA_H
#define MHD_PLASMA_H

/* Generated from src/lib.rs by cbindgen (tests/header.rs). Do not edit;
 * run `UPDATE_HEADER=1 cargo test` after changing the C ABI. The Fortran
 * declarations of the same functions are in ../MHDPlasma.F90. Link
 * libmhd_plasma.a with -lpthread -ldl -lm. */

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/*
 The seed and the model options, filled from the case's Material section.
 */
typedef struct MhdSeedPlasma {
  /*
   Seed atoms per heavy particle
   */
  double seed_frac;
  /*
   Seed ionization energy, J
   */
  double chi_j;
  /*
   g_ion / g_neutral of the seed
   */
  double weight_ratio;
  /*
   Ion mobility, m^2/(V s), at the Loschmidt density
   */
  double ion_reduced_mobility;
  /*
   Lorentz electron transport (else drifting)
   */
  bool lorentz;
  /*
   K(4p) population at the gas temperature
   */
  bool excited_at_gas;
} MhdSeedPlasma;

/*
 Local transport state at one point.
 */
typedef struct MhdPlasmaTransport {
  /*
   Electron density, m^-3
   */
  double ne;
  /*
   Coulomb logarithm
   */
  double ln_lambda;
  /*
   Electron mobilities, parallel (DC), Pedersen and Hall, m^2/(V s)
   */
  double mu0;
  double mu_p;
  double mu_h;
  /*
   Ion mobilities (Hall opposite)
   */
  double mu_i;
  double mu_ip;
  double mu_ih;
  /*
   e / (me mu0), 1/s
   */
  double collision_frequency;
} MhdPlasmaTransport;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/*
 Read the collision tables once; a later call with the same file does nothing.

 Returns 0 on success. On failure returns 1 and writes a NUL-terminated
 message of at most `message_len` bytes (including the NUL) to `message`.

 # Safety
 `file_name` must be a NUL-terminated string, and `message` must point to
 `message_len` writable bytes.
 */
int32_t mhd_plasma_read_collision_data(const char *file_name,
                                       bool lorentz,
                                       char *message,
                                       size_t message_len);

/*
 Saha electron density, m^-3, at electron temperature `t` [K] and seed density `ns`.
 */
double mhd_plasma_saha_density(const struct MhdSeedPlasma *pl, double t, double ns);

/*
 d(ne)/dT of the Saha electron density at fixed seed density, m^-3/K.
 */
double mhd_plasma_seed_density_derivative(const struct MhdSeedPlasma *pl,
                                          double t,
                                          double ns,
                                          double ne);

/*
 Electron and ion mobilities and the Coulomb logarithm at one point.
 */
void mhd_plasma_transport(const struct MhdSeedPlasma *pl,
                          double te,
                          double n,
                          double ne,
                          double bmag,
                          struct MhdPlasmaTransport *tr);

/*
 Elastic loss coefficients Cn and Cei, W/(m^3 K), for the transport state `tr`.
 */
void mhd_plasma_elastic_loss_coefficients(double te,
                                          double n,
                                          const struct MhdPlasmaTransport *tr,
                                          double *cn,
                                          double *cei);

/*
 Net K 4s-4p inelastic energy loss of the electrons, W/m^3.
 */
double mhd_plasma_inelastic_loss(const struct MhdSeedPlasma *pl,
                                 double te,
                                 double tg,
                                 double ne,
                                 double nk);

/*
 Local electron temperature, K; `at_max` is set when it is capped at `te_max`.
 */
double mhd_plasma_electron_temperature(const struct MhdSeedPlasma *pl,
                                       double tg,
                                       double n,
                                       double ns,
                                       double e_par2,
                                       double e_perp2,
                                       double bmag,
                                       double te_max,
                                       bool *at_max);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus

#endif  /* MHD_PLASMA_H */
