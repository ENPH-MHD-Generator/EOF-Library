/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    mdhLinearHall

Description
    EOF coupled OpenFOAM/Elmer solver (Hall MHD channel).

    Pattern intentionally mirrors the EOF reference solvers:
      - Construct Elmer sender/receiver normally (constructor does handshake).
      - Do ONE initial send/recv before the OF time loop.
      - During the time loop, couple every step (robust) using sendStatus(runTime.run()).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "Elmer.H"
#include "hallClosure.H"

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"
    #include "setHallClosureCoeffs.H"

    turbulence->validate();

    // Explicitly disable LTS in this solver (keep compilation simple for OF6 setups)
    const bool LTS = false;
    (void)LTS; // silence unused warning if your includes don’t reference it

    #include "readTimeControls.H" // reads time controls from the control dict
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<<"\n Initializing Te...\n" << endl;

    if (gMax(Te) <= SMALL) { Te = T; Te.correctBoundaryConditions(); }

    Info<< "\nStarting time loop\n" << endl;

    // ---------------------------------------------------------------------
    // Initial coupling (mirrors EOF test solver style)
    // ---------------------------------------------------------------------

    // Construct scalar component fields from vectors
    Ux = U.component(vector::X);
    Uy = U.component(vector::Y);
    Uz = U.component(vector::Z);

    Bx = B.component(vector::X);
    By = B.component(vector::Y);
    Bz = B.component(vector::Z);

    Info<< "\nStarting time loop\n" << endl;

    // Send fields to Elmer
    Info<< "Constructing Elmer receiving...\n" << endl;
    bool receiverInitialized = false;
    Elmer<fvMesh> receiving(mesh, -1, false); // Do not initialize yet

    Info<< "Constructing Elmer sending...\n" << endl;
    Elmer<fvMesh> sending(mesh, 1);

    Info<< "Elmer objects constructed.\n" << endl;

    Info<< "\nStarting time loop3\n" << endl;
    for (label k = 0; k < TeMaxIter; ++k)
    {
        Info<< "\nStarting time loop-1\n" << endl;
        sending.sendStatus(1); //  1 = ok / continue
        Info<< "\nStarting time loop0\n" << endl;

        // 2) Saha -> ne
        computeAlphaNe(Te, alpha, ne, rhoConstVal, saha);
        Info<< "\nStarting time loop0.5e\n" << endl;

        // 3) ne -> sigma (placeholder from alpha)
        computeSigmaFromAlpha(alpha, elcond, sigmaRef, alphaRef, sigmaMin, sigmaMax);
        Info<< "\nStarting time loop4\n" << endl;
        // Send sigma + U,B to Elmer
        Info<< "about to send elcond\n" << endl;
        sending.sendScalar(elcond);
        Info<< "sent elcond\n" << endl;
        sending.sendScalar(Ux);
        sending.sendScalar(Uy);
        sending.sendScalar(Uz);
        sending.sendScalar(Bx);
        sending.sendScalar(By);
        sending.sendScalar(Bz);
        Info<< "\nStarting time loop5\n" << endl;
        // Receive fields from Elmer

        if (!receiverInitialized) {
            receiving.initialize(); // Initialize now that we are ready
            receiverInitialized = true;
        }

        receiving.sendStatus(1);
        Info<< "\nStarting time loop6\n" << endl;
        // Receive J from Elmer
        receiving.recvScalar(Jx);
        receiving.recvScalar(Jy);
        receiving.recvScalar(Jz);
        receiving.recvScalar(JH_recv);
        JH = JH_recv;
        Info<< "\nStarting time loop7\n" << endl;
        // Reconstruct J_dens from component fields
        // Brackets define a local scope in OF6
        {
            vectorField& J_dens_in = J_dens.primitiveFieldRef();
            const scalarField & Jx_in = Jx.internalField();
            const scalarField & Jy_in = Jy.internalField();
            const scalarField & Jz_in = Jz.internalField();

            forAll(J_dens_in, celli)
            {
                J_dens_in[celli] = vector(Jx_in[celli], Jy_in[celli], Jz_in[celli]);
            }

            forAll(J_dens.boundaryField(), patchi)
            {
                vectorField& J_dens_p = J_dens.boundaryFieldRef()[patchi];
                const scalarField & Jx_p = Jx.boundaryField()[patchi];
                const scalarField & Jy_p = Jy.boundaryField()[patchi];
                const scalarField & Jz_p = Jz.boundaryField()[patchi];

                forAll(J_dens_p, facei)
                {
                    J_dens_p[facei] = vector(Jx_p[facei], Jy_p[facei], Jz_p[facei]);
                }
            }
            Info<< "\nStarting time loop8\n" << endl;
        }

        // 4) Energy balance: Cep*(Te - T) = eta*|J|^2 = |J|^2/sigma
        computeCep(ne, Cep, Cep0, CepMin);
        const dimensionedScalar sigmaMinDim("sigmaMin", elcond.dimensions(), sigmaMin);
        const dimensionedScalar CepMinDim("CepMin",  Cep.dimensions(),  CepMin);

        // TeNew = T + |J|^2/(sigma*Cep)
        // (sigma already clamped; also clamp Cep)
        TeNew = T + magSqr(J_dens) / (max(elcond, sigmaMinDim) * max(Cep, CepMinDim));
        Info<< "\nStarting time loop9\n" << endl;
        // 5) relax + clamp
        Te = (1.0 - TeOmega)*Te + TeOmega*TeNew;
        Te.max(TeMin);
        Te.min(TeMax);
        Te.correctBoundaryConditions();

        // residual check
        const scalarField dTe = TeNew.internalField() - Te.internalField();
        Info<< "\nStarting time loop10\n" << endl;
        scalar maxAbs = 0.0;
        scalar maxDen = 0.0;

        forAll(dTe, i)
        {
            maxAbs = max(maxAbs, mag(dTe[i]));
        }

        // denominator: max |Te|
        const scalarField& TeI = Te.internalField();
        forAll(TeI, i)
        {
            maxDen = max(maxDen, mag(TeI[i]));
        }
        Info<< "\nStarting time loop11\n" << endl;
        const scalar rel = maxAbs / (maxDen + SMALL);
        Info<< "Te iter " << k
            << " rel=" << rel
            << " Te[min/max]=" << gMin(Te) << " " << gMax(Te)
            << " alpha[min/max]=" << gMin(alpha) << " " << gMax(alpha)
            << " sigma[min/max]=" << gMin(elcond) << " " << gMax(elcond)
            << nl << endl;

        if (rel < TeTolRel) break;
    }

    Info<< "\nStarting time loop12\n" << endl;

    fLorentz = J_dens ^ B;
    {
        scalar fMax = 0.0;

        const vectorField& fI = fLorentz.internalField();

        forAll(fI, celli)
        {
            const scalar magFi = mag(fI[celli]);
            if (magFi > fMax)
            {
                fMax = magFi;
            }
        }

        Info<< "DEBUG fLorentz max |F| = " << fMax << nl << endl;
    }


    // ---------------------------------------------------------------------
    // OpenFOAM time loop
    // ---------------------------------------------------------------------

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "deltaT(fixed?) = " << runTime.deltaTValue() << nl << endl;

        const int status = runTime.run();

        Info<< "Coupling status = " << status
            << "  time=" << runTime.timeName()
            << "  runTime.run()=" << runTime.run()
            << nl << endl;

        receiving.sendStatus(status);
        sending.sendStatus(status);

        // -----------------------------------------------------------------
        // Coupling step EVERY time step (robust)
        // -----------------------------------------------------------------

        // Construct scalar component fields from vectors
        Ux = U.component(vector::X);
        Uy = U.component(vector::Y);
        Uz = U.component(vector::Z);

        Bx = B.component(vector::X);
        By = B.component(vector::Y);
        Bz = B.component(vector::Z);

        elcond = elcond_melt;
        Info<< "elcond min/max = " << gMin(elcond) << " " << gMax(elcond) << nl << endl;
        sending.sendScalar(elcond);
        sending.sendScalar(Ux);
        sending.sendScalar(Uy);
        sending.sendScalar(Uz);
        sending.sendScalar(Bx);
        sending.sendScalar(By);
        sending.sendScalar(Bz);

        receiving.recvScalar(Jx);
        receiving.recvScalar(Jy);
        receiving.recvScalar(Jz);
        receiving.recvScalar(JH_recv);
        JH  = JH_recv;

        // Reconstruct J_dens from component fields
        // Brackets define a local scope in OF6
        {
            vectorField& J_dens_in = J_dens.primitiveFieldRef();
            const scalarField & Jx_in = Jx.internalField();
            const scalarField & Jy_in = Jy.internalField();
            const scalarField & Jz_in = Jz.internalField();

            forAll(J_dens_in, celli)
            {
                J_dens_in[celli] = vector(Jx_in[celli], Jy_in[celli], Jz_in[celli]);
            }

            forAll(J_dens.boundaryField(), patchi)
            {
                vectorField& J_dens_p = J_dens.boundaryFieldRef()[patchi];
                const scalarField & Jx_p = Jx.boundaryField()[patchi];
                const scalarField & Jy_p = Jy.boundaryField()[patchi];
                const scalarField & Jz_p = Jz.boundaryField()[patchi];

                forAll(J_dens_p, facei)
                {
                    J_dens_p[facei] = vector(Jx_p[facei], Jy_p[facei], Jz_p[facei]);
                }
            }
        }
        fLorentz = J_dens ^ B;
        {
            scalar fMax = 0.0;

            const vectorField& fI = fLorentz.internalField();

            forAll(fI, celli)
            {
                const scalar magFi = mag(fI[celli]);
                if (magFi > fMax)
                {
                    fMax = magFi;
                }
            }

            Info<< "DEBUG fLorentz max |F| = " << fMax << nl << endl;
        }

        // -----------------------------------------------------------------
        // PIMPLE loop
        // -----------------------------------------------------------------
        while (pimple.loop())
        {
            laminarTransport.correct();

            #include "UEqn.H"
            #include "TEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}
