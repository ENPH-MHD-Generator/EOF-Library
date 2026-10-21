/* C interface of the mhd_plasma library (libs/solvers/MHDSolve/plasma).
 *
 * Kept by hand in step with src/lib.rs; the Fortran interface module
 * MHDPlasma.F90 declares the same functions. Link libmhd_plasma.a together
 * with -lpthread -ldl -lm.
 */
#ifndef MHD_PLASMA_H
#define MHD_PLASMA_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The seed and the model options */
typedef struct {
    double seed_frac;            /* seed atoms per heavy particle */
    double chi_j;                /* seed ionization energy, J */
    double weight_ratio;         /* g_ion / g_neutral of the seed */
    double ion_reduced_mobility; /* m^2/(V s) at the Loschmidt density */
    bool lorentz;                /* Lorentz electron transport (else drifting) */
    bool excited_at_gas;         /* K(4p) population at the gas temperature */
} MhdSeedPlasma;

/* Local transport state at one point */
typedef struct {
    double ne;                   /* electron density, m^-3 */
    double ln_lambda;            /* Coulomb logarithm */
    double mu0, mu_p, mu_h;      /* electron mobilities, m^2/(V s) */
    double mu_i, mu_ip, mu_ih;   /* ion mobilities (Hall opposite) */
    double collision_frequency;  /* e / (me mu0), 1/s */
} MhdPlasmaTransport;

int32_t mhd_plasma_read_collision_data(const char *file_name, bool lorentz,
                                       char *message, size_t message_len);
double mhd_plasma_saha_density(const MhdSeedPlasma *pl, double t, double ns);
double mhd_plasma_seed_density_derivative(const MhdSeedPlasma *pl, double t,
                                          double ns, double ne);
void mhd_plasma_transport(const MhdSeedPlasma *pl, double te, double n, double ne,
                          double bmag, MhdPlasmaTransport *tr);
void mhd_plasma_elastic_loss_coefficients(double te, double n,
                                          const MhdPlasmaTransport *tr,
                                          double *cn, double *cei);
double mhd_plasma_inelastic_loss(const MhdSeedPlasma *pl, double te, double tg,
                                 double ne, double nk);
double mhd_plasma_electron_temperature(const MhdSeedPlasma *pl, double tg, double n,
                                       double ns, double e_par2, double e_perp2,
                                       double bmag, double te_max, bool *at_max);

#ifdef __cplusplus
}
#endif

#endif
