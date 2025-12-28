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
    EOF coupled OpenFOAM/Elmer solver.

    IMPORTANT:
      - Disable LTS (local time stepping) to avoid createRDeltaT/localEulerDdt
        compile issues in your current OF6 build setup.
      - Fix EOF coupling deadlock by deferring receiver initialization
        (init=false) and NEVER doing blocking recv() before the time loop.
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

    // ------------------------------------------------------------
    // Disable local time stepping (LTS) explicitly.
    // This avoids createRDeltaT.H / localEulerDdt / fvc::smooth issues.
    // ------------------------------------------------------------
    const bool LTS = false;

    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    // ------------------------------------------------------------
    // EOF coupling objects
    //
    // sending: init immediately is OK (Elmer starts in OpenFOAM2Elmer stage)
    // receiving: MUST NOT init immediately (would deadlock in initialize()).
    // ------------------------------------------------------------
    Elmer<fvMesh> sending(mesh,  1, true);     // mode=+1, init now
    Elmer<fvMesh> receiving(mesh, -1, false);  // mode=-1, init later
    bool receivingInitialized = false;

    // If you have JxB_recv/JH_recv in createFields.H, keep them.
    // Otherwise, receive directly into JxB/JH.
    // Here I assume your createFields.H defines: JxB_recv, JH_recv, JxB, JH.
    // If not, adjust names accordingly.

    // Optional: if you want to always couple every time step, keep as-is.
    // If you want a trigger, base it on U change (as you started doing).

    while (runTime.run())
    {
        // No LTS branch at all (keeps compilation simple)
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ------------------------------------------------------------
        // Coupling step EVERY time step (robust)
        // ------------------------------------------------------------
        sending.sendStatus(1);

        // Your choice: either use the existing elcond field directly
        // or set it each step. Keep it simple:
        elcond = elcond_melt;
        Info<< "elcond min/max = " << gMin(elcond) << " " << gMax(elcond) << nl << endl;
        sending.sendScalar(elcond);

        // Init reverse mapping ONCE, but only after the first sendScalar()
        if (!receivingInitialized)
        {
            receiving.initialize();
            receivingInitialized = true;
        }

        receiving.sendStatus(1);
        Info<< "BEFORE recvVector" << endl;
        receiving.recvVector(JxB_recv);
        Info<< "AFTER recvVector" << endl;

        Info<< "BEFORE recvScalar" << endl;
        receiving.recvScalar(JH_recv);
        Info<< "AFTER recvScalar" << endl;


        // Apply received fields
        JxB = JxB_recv;
        JH  = JH_recv;

        // ------------------------------------------------------------
        // PIMPLE loop
        // ------------------------------------------------------------
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

    // Clean shutdown
    sending.sendStatus(0);
    if (receivingInitialized)
    {
        receiving.sendStatus(0);
    }

    Info<< "End\n" << endl;
    return 0;
}

// ************************************************************************* //
