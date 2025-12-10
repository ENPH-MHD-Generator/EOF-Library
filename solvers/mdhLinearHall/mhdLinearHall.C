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

    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh,1); // 1=send, -1=receive
    sending.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    elcond = elcond_melt;
    sending.sendScalar(elcond);

    // Receive fields from Elmer
    Elmer<fvMesh> receiving(mesh,-1); // 1=send, -1=receive
    receiving.sendStatus(1); // 1=ok, 0=lastIter, -1=error
    receiving.recvVector(JxB_recv);
    receiving.recvScalar(JH_recv);

    while (runTime.run())
    {
        JxB = JxB_recv;
        JH = JH_recv;

        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            laminarTransport.correct();

            #include "UEqn.H"

            #include "TEqn.H"

            // --- Pressure corrector loop
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

        bool doElmer = false;
        if
        (
            maxRelDiff_local > maxRelDiff
         && (maxRelDiff < SMALL || maxRelDiff + SMALL <= 1.0)
        )
        {
            doElmer = true;
        }

        if (doElmer || !runTime.run())
        {
            U_old = U;

            // Send fields to Elmer
            sending.sendStatus(runTime.run());
            elcond = elcond_melt;
            sending.sendScalar(elcond);

            // Receive fields from Elmer
            receiving.sendStatus(runTime.run());
            receiving.recvVector(JxB_recv);
            receiving.recvScalar(JH_recv);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //


