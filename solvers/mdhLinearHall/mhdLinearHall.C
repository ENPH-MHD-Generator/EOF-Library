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

    Notes:
      - This assumes your Elmer SIF exports:
          Target Variable 1 = (vector) J_dens
          Target Variable 2 = (scalar) Joule Heating
        i.e., the receive order below matches the SIF target order.
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


    // For coupling stabilization.
    // TODO: Move this to a centralized config
    const scalar alphaLorentz = 0.2;   


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
    const volVectorField fNew(J_dens ^ B);
    fLorentz = alphaLorentz*fNew + (1.0 - alphaLorentz)*fLorentz;

    // ---------------------------------------------------------------------
    // OpenFOAM time loop
    // ---------------------------------------------------------------------

    const int numIter = 5000;
    for (int iter = 0; iter < numIter; ++iter)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "deltaT(fixed?) = " << runTime.deltaTValue() << nl << endl;

        const int status = (iter < numIter-1) ? 1 : 0;

        Info<< "Coupling status = " << status
            << "  time=" << runTime.timeName()
            << "  runTime.run()=" << runTime.run()
            << nl << endl;

        receiving.sendStatus(status);
        sending.sendStatus(status);

        // -----------------------------------------------------------------
        // Coupling step EVERY time step (robust)
        // -----------------------------------------------------------------

        elcond = elcond_melt;
        Info<< "elcond min/max = " << gMin(elcond) << " " << gMax(elcond) << nl << endl;
        sending.sendScalar(elcond);

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
        const volVectorField fNew(J_dens ^ B);
        fLorentz = alphaLorentz*fNew + (1.0 - alphaLorentz)*fLorentz;

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

// ************************************************************************* //
