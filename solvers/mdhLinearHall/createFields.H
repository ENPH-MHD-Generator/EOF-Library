// -------------------------------------------------------------------------
// Primary flow fields
// -------------------------------------------------------------------------

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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
const scalar couplingReferencePressure =
    couplingProperties.lookupOrDefault<scalar>("referencePressure", 101325);

Info<< "Elmer coupling tolerances: U " << couplingVelocityTolerance
    << "  T " << couplingTemperatureTolerance
    << "  p " << couplingPressureTolerance
    << "  max steps between updates " << couplingMaxSteps << nl << endl;

label stepsSinceElmer = 0;
label nElmerUpdates = 0;

Info<< "Initializing field U_old\n" << endl;
volVectorField U_old
(
    IOobject
    (
        "U_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);

// -------------------------------------------------------------------------
// MHD coupling fields
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// MHD coupling fields (Elmer -> OpenFOAM)
// Elmer exports Volume Current as components: Jx, Jy, Jz
// We reconstruct J as a volVectorField
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


Info<< "Initializing current density components (Jx, Jy, Jz)\n" << endl;

// Fields received from Elmer only have cell values; zeroGradient patches copy
// them to the boundary (correctBoundaryConditions), which keeps the boundary
// and the VTK point data from reading as zero.

volScalarField Jx
(
    IOobject
    (
        "Jx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrent/dimArea,   
        0.0
    ),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Jy
(
    IOobject
    (
        "Jy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrent/dimArea,
        0.0
    ),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Jz
(
    IOobject
    (
        "Jz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrent/dimArea,
        0.0
    ),
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
    dimensionedVector
    (
        "zero",
        dimCurrent/dimArea,
        vector::zero
    )
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
    J_dens ^ B // This is cross product for openFOAM
);



Info<< "Initializing field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

Info<< "Initializing field JH_recv\n" << endl;
volScalarField JH_recv
(
    IOobject
    (
        "JH_recv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    JH
);

// -------------------------------------------------------------------------
// Fluxes and transport
// -------------------------------------------------------------------------

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
singlePhaseTransportModel laminarTransport(U, phi);

// Constant fluid properties
dimensionedScalar rhoConst("rho", dimDensity, laminarTransport);
dimensionedScalar cp("cp", dimSpecificHeatCapacity, laminarTransport);
dimensionedScalar kappa("kappa", dimPower/dimLength/dimTemperature, laminarTransport);

// Density field
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    rhoConst
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// -------------------------------------------------------------------------
// Gravity and pressure
// -------------------------------------------------------------------------

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

// Temperature and pressure at the last Elmer update
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

label pRefCell = 0;
scalar pRefValue = 0.0;

setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());

// -------------------------------------------------------------------------
// MRF (after mesh & transport are ready)
// -------------------------------------------------------------------------

#include "createMRF.H"


// -------------------------------------------------------------------------
// Elmer update control
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Scalar velocity components for Elmer coupling
// -------------------------------------------------------------------------

Info<< "Initializing velocity components (Ux, Uy, Uz)\n" << endl;

volScalarField Ux
(
    IOobject
    (
        "Ux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimVelocity,
        0.0
    )
);

volScalarField Uy
(
    IOobject
    (
        "Uy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimVelocity,
        0.0
    )
);

volScalarField Uz
(
    IOobject
    (
        "Uz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimVelocity,
        0.0
    )
);

// Initialize once
Ux = U.component(vector::X);
Uy = U.component(vector::Y);
Uz = U.component(vector::Z);


// -------------------------------------------------------------------------
// Scalar magnetic field components for Elmer coupling
// -------------------------------------------------------------------------

Info<< "Initializing magnetic field components (Bx, By, Bz)\n" << endl;

volScalarField Bx
(
    IOobject
    (
        "Bx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        B.dimensions(),
        0.0
    )
);

volScalarField By
(
    IOobject
    (
        "By",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        B.dimensions(),
        0.0
    )
);

volScalarField Bz
(
    IOobject
    (
        "Bz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        B.dimensions(),
        0.0
    )
);

// Initialize once
Bx = B.component(vector::X);
By = B.component(vector::Y);
Bz = B.component(vector::Z);

// -------------------------------------------------------------------------
// Seeded-plasma state computed by Elmer (Saha equation), written for viewing
// -------------------------------------------------------------------------

Info<< "Initializing plasma fields (elcond_elmer, ionizationFraction, Te)\n" << endl;

volScalarField elcond_elmer
(
    IOobject
    (
        "elcond_elmer",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// Electron mole fraction ne / n_heavy
volScalarField ionizationFraction
(
    IOobject
    (
        "ionizationFraction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Te
(
    IOobject
    (
        "Te",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
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
    IOobject
    (
        "lorentzDamping",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    elcond_elmer*magSqr(B)
);
