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

Info<< "Initializing field JxB\n" << endl;
volVectorField JxB
(
    IOobject
    (
        "JxB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field JxB_recv\n" << endl;
volVectorField JxB_recv
(
    IOobject
    (
        "JxB_recv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    JxB
);

Info<< "Initializing field elcond\n" << endl;
volScalarField elcond
(
    IOobject
    (
        "elcond",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimensionSet(-1, -3, 3, 0, 0, 2, 0),
        0.0
    )
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
// Electrical conductivity (single-body enforced)
// -------------------------------------------------------------------------

dimensionedScalar elcond_melt
(
    "elcond_melt",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    laminarTransport
);

// Force single, nonzero conductivity everywhere
elcond = elcond_melt;

// -------------------------------------------------------------------------
// Elmer update control
// -------------------------------------------------------------------------

scalar maxRelDiff
(
    readScalar(laminarTransport.lookup("maxRelDiff"))
);
