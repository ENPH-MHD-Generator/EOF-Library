// -------------------------------------------------------------------------
// Compressible gas (rhoPimpleFoam): thermophysical model, density, velocity
// -------------------------------------------------------------------------

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

// Absolute pressure [Pa]
volScalarField& p = thermo.p();
// Gas temperature from the energy equation; sent to Elmer
volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

// -------------------------------------------------------------------------
// Coupling controls (constant/couplingProperties, written by caseprep).
// Elmer is re-solved when U, T or p has changed by more than these relative
// tolerances since the last update, or after maxStepsBetweenUpdates steps.
// Zero tolerances update every step.
// -------------------------------------------------------------------------

IOdictionary couplingProperties
(
    IOobject
    (
        "couplingProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const scalar couplingVelocityTolerance =
    couplingProperties.lookupOrDefault<scalar>("velocityTolerance", 0);
const scalar couplingTemperatureTolerance =
    couplingProperties.lookupOrDefault<scalar>("temperatureTolerance", 0);
const scalar couplingPressureTolerance =
    couplingProperties.lookupOrDefault<scalar>("pressureTolerance", 0);
const label couplingMaxSteps =
    couplingProperties.lookupOrDefault<label>("maxStepsBetweenUpdates", 0);

Info<< "Elmer coupling tolerances: U " << couplingVelocityTolerance
    << "  T " << couplingTemperatureTolerance
    << "  p " << couplingPressureTolerance
    << "  max steps between updates " << couplingMaxSteps << nl << endl;

label stepsSinceElmer = 0;
label nElmerUpdates = 0;

// Velocity, temperature and pressure at the last Elmer update
volVectorField U_old
(
    IOobject("U_old", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    U
);
volScalarField T_sent
(
    IOobject("T_sent", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    T
);
volScalarField p_sent
(
    IOobject("p_sent", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    p
);

// -------------------------------------------------------------------------
// MHD coupling fields
// -------------------------------------------------------------------------

Info<< "Reading magnetic field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading electric potential\n" << endl;
volScalarField potential
(
    IOobject
    (
        "Potential",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing electric field (E = -grad(potential))\n" << endl;
volVectorField electric_field
(
    IOobject
    (
        "electric_field",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Fields received from Elmer only have cell values; zeroGradient patches copy
// them to the boundary (correctBoundaryConditions), which keeps the boundary
// and the VTK point data from reading as zero.

Info<< "Initializing current density components (Jx, Jy, Jz)\n" << endl;
volScalarField Jx
(
    IOobject("Jx", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimCurrent/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField Jy
(
    IOobject("Jy", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimCurrent/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField Jz
(
    IOobject("Jz", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimCurrent/dimArea, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Initializing reconstructed current vector J (from components)\n" << endl;
volVectorField J_dens
(
    IOobject
    (
        "J",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimCurrent/dimArea, vector::zero)
);

Info<< "Initializing Lorentz force field fLorentz\n" << endl;
volVectorField fLorentz
(
    IOobject
    (
        "fLorentz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    J_dens ^ B
);

// Joule heating J^2/sigma [W/m^3]
Info<< "Initializing field JH\n" << endl;
volScalarField JH
(
    IOobject
    (
        "JH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField JH_recv
(
    IOobject("JH_recv", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    JH
);

// -------------------------------------------------------------------------
// Scalar velocity and magnetic field components for Elmer coupling
// -------------------------------------------------------------------------

Info<< "Initializing velocity and magnetic field components\n" << endl;
volScalarField Ux
(
    IOobject("Ux", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);
volScalarField Uy
(
    IOobject("Uy", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);
volScalarField Uz
(
    IOobject("Uz", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);
volScalarField Bx
(
    IOobject("Bx", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", B.dimensions(), 0.0)
);
volScalarField By
(
    IOobject("By", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", B.dimensions(), 0.0)
);
volScalarField Bz
(
    IOobject("Bz", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", B.dimensions(), 0.0)
);

// -------------------------------------------------------------------------
// Seeded-plasma state computed by Elmer (Saha equation), written for viewing
// -------------------------------------------------------------------------

Info<< "Initializing plasma fields (elcond_elmer, ionizationFraction, Te)\n" << endl;
volScalarField elcond_elmer
(
    IOobject("elcond_elmer", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// Electron mole fraction ne / n_heavy
volScalarField ionizationFraction
(
    IOobject("ionizationFraction", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Te
(
    IOobject("Te", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// Semi-implicit Lorentz drag coefficient K = sigma |B|^2 [kg/(m^3 s)]. The
// momentum equation adds -K (U - U_old) to the lagged J x B force: it damps
// the explicit coupling in transients and vanishes at convergence. It ignores
// the Hall reduction sigma B^2/(1 + beta^2), so it over-damps (safely).
volScalarField lorentzDamping
(
    IOobject("lorentzDamping", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    elcond_elmer*magSqr(B)
);

#include "createMRF.H"
#include "createFvOptions.H"
