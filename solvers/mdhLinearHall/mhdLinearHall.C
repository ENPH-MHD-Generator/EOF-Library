/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mdhLinearHall = Elmer (EM) + interFoam + TEqn + TableVisc

Description
    Solver for 2 incompressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

-------------------------------------------------------------------------------

Original interFoam solver is part of OpenFOAM

Thermal solver taken from Qingming Liu:
http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2011/QingmingLiu/Project_QingmingLIU-final.pdf

Modified by: Juris Vencels

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "Elmer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    Info<< "Debug: entering initial Elmer coupling before first time step" << nl << endl;

    // Send fields to Elmer
    Info<< "Debug: creating Elmer sender and sending initial elcond field" << nl << endl;
    Elmer<fvMesh> sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    elcond = elcond_melt;
    sending.sendScalar(elcond);

    if (Pstream::master())
    {
        Info<< "Initial Elmer send (elcond) completed" << nl << endl;
    }

    // Receive fields from Elmer
    Info<< "Debug: creating Elmer receiver and waiting for initial JxB/JH" << nl << endl;
    Elmer<fvMesh> receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(JxB_recv);
    receiving.recvScalar(JH_recv);

    if (Pstream::master())
    {
        Info<< "Initial Elmer receive (JxB, JH) completed" << nl << endl;
    }

    while (runTime.run())
    {
        Info<< "Debug: top of time-step while(runTime.run()), time = "
            << runTime.timeName() << nl << endl;

        JxB = JxB_recv;
        JH = JH_recv;

        if (Pstream::master())
        {
            Info<< "Entering time step, current time = " << runTime.timeName()
                << ", maxRelDiff trigger = " << maxRelDiff << nl << endl;
        }

        Info<< "Debug: about to read time controls" << nl << endl;
        #include "readTimeControls.H"
        Info<< "Debug: finished readTimeControls.H, current deltaT = "
            << runTime.deltaTValue() << nl << endl;

        if (LTS)
        {
            Info<< "Debug: using local time stepping, entering setRDeltaT.H" << nl << endl;
            #include "setRDeltaT.H"
            Info<< "Debug: finished setRDeltaT.H" << nl << endl;
        }
        else
        {
            Info<< "Debug: computing Courant number" << nl << endl;
            #include "CourantNo.H"
            Info<< "Debug: finished CourantNo.H, entering setDeltaT.H" << nl << endl;
            #include "setDeltaT.H"
            Info<< "Debug: finished setDeltaT.H, new deltaT = "
                << runTime.deltaTValue() << nl << endl;
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        Info<< "Debug: entering PIMPLE loop for this time step" << nl << endl;
        while (pimple.loop())
        {
            Info<< "Debug:  PIMPLE outer iteration, calling laminarTransport.correct()" << nl << endl;
            laminarTransport.correct();

            Info<< "Debug:  about to assemble/solve UEqn (UEqn.H)" << nl << endl;
            #include "UEqn.H"
            Info<< "Debug:  finished UEqn.H" << nl << endl;

            Info<< "Debug:  about to assemble/solve TEqn (TEqn.H)" << nl << endl;
            #include "TEqn.H"
            Info<< "Debug:  finished TEqn.H" << nl << endl;

            // --- Pressure corrector loop
            Info<< "Debug:  entering pressure corrector loop (pEqn.H)" << nl << endl;
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
            Info<< "Debug:  finished pressure corrector loop" << nl << endl;

            if (pimple.turbCorr())
            {
                Info<< "Debug:  turbulence correction enabled, calling turbulence->correct()" << nl << endl;
                turbulence->correct();
            }
        }
        Info<< "Debug: leaving PIMPLE loop for this time step" << nl << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        // Check whether we need to update electromagnetic solution with Elmer.
        // Use relative change in velocity magnitude as a trigger, similar to
        // the mhdVxBPimpleFoam solver. If maxRelDiff == 0, always update.
        dimensionedScalar smallU
        (
            "smallU",
            dimensionSet(0, 1, -1, 0, 0, 0 ,0),
            1e-6
        );

        scalar maxRelDiff_local =
            (max(mag(U_old - U)/(average(mag(U)) + smallU))).value();

        if (Pstream::master())
        {
            Info<< "Computed maxRelDiff_local = " << maxRelDiff_local << nl << endl;
        }

        bool doElmer = false;
        if
        (
            maxRelDiff_local > maxRelDiff
         && (maxRelDiff < SMALL || maxRelDiff + SMALL <= 1.0)
        )
        {
            doElmer = true;

            if (Pstream::master())
            {
                Info<< "Elmer update triggered at time = " << runTime.timeName()
                    << " (maxRelDiff_local = " << maxRelDiff_local << ")" << nl << endl;
            }
        }

        if (doElmer || !runTime.run())
        {
            if (Pstream::master())
            {
                Info<< "Performing Elmer coupling step at time = "
                    << runTime.timeName() << nl << endl;
            }

            U_old = U;

            // Send fields to Elmer
            sending.sendStatus(runTime.run());
            elcond = elcond_melt;
            sending.sendScalar(elcond);

            // Receive fields from Elmer
            receiving.sendStatus(runTime.run());
            receiving.recvVector(JxB_recv);
            receiving.recvScalar(JH_recv);

            if (Pstream::master())
            {
                Info<< "Elmer coupling step completed at time = "
                    << runTime.timeName() << nl << endl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //


