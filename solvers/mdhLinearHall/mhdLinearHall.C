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
      - No manual .initialize() calls.

    Notes:
      - This assumes your Elmer SIF exports:
          Target Variable 1 = (vector) JxB
          Target Variable 2 = (scalar) Joule Heating
        i.e., the receive order below matches the SIF target order.

      - If your Elmer side actually exports "Volume Current" (J) rather than JxB,
        then you should rename fields accordingly OR compute JxB = J ^ B in OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "Elmer.H"

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

    turbulence->validate();

    // Explicitly disable LTS in this solver (keep compilation simple for OF6 setups)
    const bool LTS = false;
    (void)LTS; // silence unused warning if your includes don’t reference it

    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    // ---------------------------------------------------------------------
    // Initial coupling (mirrors EOF test solver style)
    // ---------------------------------------------------------------------

    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh, 1);     //  1 = send
    sending.sendStatus(1);              //  1 = ok / continue
    elcond = elcond_melt;
    sending.sendScalar(elcond);

    // Receive fields from Elmer
    Elmer<fvMesh> receiving(mesh, -1);  // -1 = receive
    receiving.sendStatus(1);

    // Receive in the SAME order Elmer exports "Target Variable i"
    receiving.recvScalar(Jx);
    receiving.recvScalar(Jy);
    receiving.recvScalar(Jz);
    receiving.recvScalar(JH_recv);
    // Reconstruct JxB from component fields
    // Brackets define a local scope in OF6
    {
        vectorField& JxBif = JxB.primitiveFieldRef();
        const scalarField& Jxif = Jx.internalField();
        const scalarField& Jyif = Jy.internalField();
        const scalarField& Jzif = Jz.internalField();

        forAll(JxBif, celli)
        {
            JxBif[celli] = vector(Jxif[celli], Jyif[celli], Jzif[celli]);
        }

        forAll(JxB.boundaryField(), patchi)
        {
            vectorField& JxBp = JxB.boundaryFieldRef()[patchi];
            const scalarField& Jxp = Jx.boundaryField()[patchi];
            const scalarField& Jyp = Jy.boundaryField()[patchi];
            const scalarField& Jzp = Jz.boundaryField()[patchi];

            forAll(JxBp, facei)
            {
                JxBp[facei] = vector(Jxp[facei], Jyp[facei], Jzp[facei]);
            }
        }
    }

    receiving.recvScalar(JH_recv);      // only if Elmer exports Joule Heating
    JH  = JH_recv;

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

        // -----------------------------------------------------------------
        // Coupling step EVERY time step (robust)
        // IMPORTANT: use sendStatus(runTime.run()) to match EOF semantics.
        // -----------------------------------------------------------------
        sending.sendStatus(runTime.run());

        elcond = elcond_melt;
        Info<< "elcond min/max = " << gMin(elcond) << " " << gMax(elcond) << nl << endl;
        sending.sendScalar(elcond);

        receiving.sendStatus(runTime.run());

        receiving.recvScalar(Jx);
        receiving.recvScalar(Jy);
        receiving.recvScalar(Jz);
        receiving.recvScalar(JH_recv);

        // Reconstruct JxB from component fields
        // Brackets define a local scope in OF6
        {
            vectorField& JxBif = JxB.primitiveFieldRef();
            const scalarField& Jxif = Jx.internalField();
            const scalarField& Jyif = Jy.internalField();
            const scalarField& Jzif = Jz.internalField();

            forAll(JxBif, celli)
            {
                JxBif[celli] = vector(Jxif[celli], Jyif[celli], Jzif[celli]);
            }

            // If you care about boundaries too, set them similarly:
            forAll(JxB.boundaryField(), patchi)
            {
                vectorField& JxBp = JxB.boundaryFieldRef()[patchi];
                const scalarField& Jxp = Jx.boundaryField()[patchi];
                const scalarField& Jyp = Jy.boundaryField()[patchi];
                const scalarField& Jzp = Jz.boundaryField()[patchi];

                forAll(JxBp, facei)
                {
                    JxBp[facei] = vector(Jxp[facei], Jyp[facei], Jzp[facei]);
                }
            }
        }
        JH  = JH_recv;

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

    // ---------------------------------------------------------------------
    // Clean shutdown (match EOF convention)
    // ---------------------------------------------------------------------
    sending.sendStatus(0);
    receiving.sendStatus(0);

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
