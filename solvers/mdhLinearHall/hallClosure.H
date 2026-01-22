//
// Created by Joshua Riefman on 2026-01-20.
//

#ifndef EOF_LIBRARY_HALLCLOSURE_H
#define EOF_LIBRARY_HALLCLOSURE_H

namespace
{
struct SahaCoeffs
{
    scalar molarMass;   // kg/mol
    scalar chi_eV;      // eV
    scalar gi_over_gn;  // lumped ratio
    scalar alphaMin;
    scalar alphaMax;
};

inline scalar eV_to_J(const scalar eV)
{
    return eV * 1.602176634e-19;
}

inline scalar saha_S(const scalar Te, const SahaCoeffs& c)
{
    const scalar me = 9.1093837015e-31;
    const scalar kB = 1.380649e-23;
    const scalar h  = 6.62607015e-34;

    const scalar chi = eV_to_J(c.chi_eV);

    const scalar pref = Foam::pow((2.0*constant::mathematical::pi * me * kB * Te)/(h*h), 1.5);

    return pref * c.gi_over_gn * Foam::exp(-chi/(kB*Te));
}

// Returns ionization fraction alpha in [0,1] for single-stage Saha
inline scalar saha_alpha_from_rho(const scalar Te, const scalar rho, const SahaCoeffs& c)
{
    const scalar NA = 6.02214076e23;
    const scalar ntot = (rho / c.molarMass) * NA;  // 1/m^3

    if (ntot <= SMALL) return 0.0;

    const scalar S = saha_S(Te, c);
    if (S <= SMALL) return 0.0;

    // ne^2/(ntot-ne)=S => ne^2 + S ne - S ntot = 0
    const scalar disc = S*S + 4.0*S*ntot;
    if (disc <= 0) return 0.0;

    const scalar ne = 0.5 * (-S + Foam::sqrt(disc));
    scalar alpha = ne/ntot;

    alpha = min(max(alpha, c.alphaMin), c.alphaMax);
    return alpha;
}

void computeAlphaNe
(
    const volScalarField& Te,
    volScalarField& alpha,
    volScalarField& ne,
    const scalar rhoConst,
    const SahaCoeffs& saha
)
{
    const scalar NA = 6.02214076e23;
    const scalar ntot = (rhoConst / saha.molarMass) * NA;

    scalarField& aI  = alpha.primitiveFieldRef();
    scalarField& neI = ne.primitiveFieldRef();
    const scalarField& TeI = Te.internalField();

    forAll(TeI, celli)
    {
        const scalar Tec = max(TeI[celli], SMALL);
        const scalar a = saha_alpha_from_rho(Tec, rhoConst, saha);
        aI[celli]  = a;
        neI[celli] = a * ntot;
    }

    alpha.correctBoundaryConditions();
    ne.correctBoundaryConditions();
}

// Placeholder conductivity model: sigma = sigmaRef * clamp(alpha/alphaRef)
// Replace later with collision frequency based sigma.
void computeSigmaFromAlpha
(
    const volScalarField& alpha,
    volScalarField& sigma,
    const scalar sigmaRef,
    const scalar alphaRef,
    const scalar sigmaMin,
    const scalar sigmaMax
)
{
    const dimensionedScalar sigmaRefDim("sigmaRef", sigma.dimensions(), sigmaRef);
    const dimensionedScalar sigmaMinDim("sigmaMin", sigma.dimensions(), sigmaMin);
    const dimensionedScalar sigmaMaxDim("sigmaMax", sigma.dimensions(), sigmaMax);

Info<< "alphaRef=" << alphaRef << " sigmaRef=" << sigmaRef << endl;
const scalarField& a = alpha.internalField();
scalar aMin = GREAT, aMax = -GREAT;
forAll(a, i) { aMin = min(aMin, a[i]); aMax = max(aMax, a[i]); }
Info<< "alpha local min/max=" << aMin << " " << aMax << endl;

    sigma = sigmaRefDim * (alpha / max(alphaRef, SMALL));  // now RHS has conductivity dims
    sigma.max(sigmaMinDim);
    sigma.min(sigmaMaxDim);
    //sigma.correctBoundaryConditions();
}

// Placeholder energy exchange coefficient: Cep = Cep0 * ne
// Replace later with proper electron-heavy collision energy transfer.
void computeCep
(
    const volScalarField& ne,
    volScalarField& Cep,
    const scalar Cep0,
    const scalar CepMin
)
{
    const dimensionedScalar Cep0Dim("Cep0", Cep.dimensions()*dimVolume, Cep0);
    // Explanation: Cep0 must have units so that Cep0 * ne -> Cep.
    // Since ne is 1/m^3, Cep0Dim should be (Cep * m^3).

    const dimensionedScalar CepMinDim("CepMin", Cep.dimensions(), CepMin);

    Cep = Cep0Dim * ne;
    Cep.max(CepMinDim);
    Cep.correctBoundaryConditions();
}
} // namespace

#endif //EOF_LIBRARY_HALLCLOSURE_H