//
// Created by Joshua Riefman on 2026-01-20.
//

#ifndef EOF_LIBRARY_SETHALLCLOSURECOEFFS_H

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Sub-dictionaries
const dictionary& sahaDict = transportProperties.subDict("saha");
const dictionary& sigDict  = transportProperties.subDict("sigmaModel");
const dictionary& exDict   = transportProperties.subDict("energyExchange");
const dictionary& teDict   = transportProperties.subDict("TeCoupling");

// Saha coefficients
SahaCoeffs saha;
saha.molarMass  = readScalar(sahaDict.lookup("molarMass"));
saha.chi_eV     = readScalar(sahaDict.lookup("chi_eV"));
saha.gi_over_gn = readScalar(sahaDict.lookup("gi_over_gn"));
saha.alphaMin   = readScalar(sahaDict.lookup("alphaMin"));
saha.alphaMax   = readScalar(sahaDict.lookup("alphaMax"));

// Sigma knobs
const scalar sigmaMin = readScalar(sigDict.lookup("sigmaMin"));
const scalar sigmaMax = readScalar(sigDict.lookup("sigmaMax"));
const scalar sigmaRef = readScalar(sigDict.lookup("sigmaRef"));
const scalar alphaRef = readScalar(sigDict.lookup("alphaRef"));

// Exchange knobs
const scalar Cep0   = readScalar(exDict.lookup("Cep0"));
const scalar CepMin = readScalar(exDict.lookup("CepMin"));


// Iteration controls
const label  TeMaxIter = readLabel(teDict.lookup("maxIter"));
const scalar TeTolRel  = readScalar(teDict.lookup("tolRel"));
const scalar TeOmega   = readScalar(teDict.lookup("omega"));
const scalar TeMin     = readScalar(teDict.lookup("TeMin"));
const scalar TeMax     = readScalar(teDict.lookup("TeMax"));

// Constant rho from your transportProperties
const scalar rhoConstVal = readScalar(transportProperties.lookup("rho"));

#define EOF_LIBRARY_SETHALLCLOSURECOEFFS_H

#endif //EOF_LIBRARY_SETHALLCLOSURECOEFFS_H