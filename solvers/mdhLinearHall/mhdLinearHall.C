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

    Compressible laminar gas flow (rhoPimpleFoam, static mesh) with the
    Lorentz force J x B in the momentum equation and the electromagnetic power
    J.E in the total energy equation. Elmer solves the electrical problem and
    the seeded-plasma state from the velocity, gas temperature and absolute
    pressure sent to it.

    Pattern intentionally mirrors the EOF reference solvers:
      - Construct Elmer sender/receiver normally (constructor does handshake).
      - Do ONE initial send/recv before the OF time loop.
      - During the time loop, couple every step (robust) using sendStatus(runTime.run()).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "Elmer.H"
#include "zeroGradientFvPatchFields.H"

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

    turbulence->validate();

    #include "readTimeControls.H" // reads time controls from the control dict
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

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


    // Send fields to Elmer
    Elmer<fvMesh> sending(mesh, 1);     //  1 = send
    sending.sendStatus(1);              //  1 = ok / continue
    U_old = U;  // velocity Elmer computes J from; used by the Lorentz damping
    sending.sendScalar(Ux);
    sending.sendScalar(Uy);
    sending.sendScalar(Uz);
    sending.sendScalar(Bx);
    sending.sendScalar(By);
    sending.sendScalar(Bz);
	sending.sendScalar(p);
	sending.sendScalar(T);

    // Receive fields from Elmer
    Elmer<fvMesh> receiving(mesh, -1);  // -1 = receive
    receiving.sendStatus(1);

    receiving.recvScalar(Jx);
    receiving.recvScalar(Jy);
    receiving.recvScalar(Jz);
    receiving.recvScalar(JH_recv);
    receiving.recvScalar(potential);
    JH  = JH_recv;
    receiving.recvScalar(elcond_elmer);
    receiving.recvScalar(ionizationFraction);
    receiving.recvScalar(Te);

    // Without this, parallel decomposition fails when writing fields with non-zero values
    Jx.correctBoundaryConditions();
    Jy.correctBoundaryConditions();
    Jz.correctBoundaryConditions();
    JH.correctBoundaryConditions();
    potential.correctBoundaryConditions();
    elcond_elmer.correctBoundaryConditions();
    ionizationFraction.correctBoundaryConditions();
    Te.correctBoundaryConditions();

    // Compute electric field E = -grad(potential)
    electric_field = -fvc::grad(potential);
    electric_field.correctBoundaryConditions();

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
    lorentzDamping = elcond_elmer*magSqr(B);
    {
        // Mechanical power the flow loses to the field; matches Elmer's P_emf
        const scalar Pmech = -gSum
        (
            (U.primitiveField() & fLorentz.primitiveField())*mesh.V().field()
        );
        Info<< "Lorentz force: max |F| = " << gMax(mag(fLorentz.primitiveField()))
            << " N/m^3  P_mech = " << Pmech << " W" << nl << endl;
    }


    // ---------------------------------------------------------------------
    // OpenFOAM time loop
    // ---------------------------------------------------------------------

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "deltaT(fixed?) = " << runTime.deltaTValue() << nl << endl;

        const int status = runTime.run();

        // -----------------------------------------------------------------
        // Decide whether to update the electrical solution. Elmer's problem
        // is quasi-static (J depends only on the current U, T, p), so it only
        // needs re-solving once its inputs have changed by more than the
        // tolerances in constant/couplingProperties since the last update.
        // Always update on the last step so Elmer receives the final status.
        // -----------------------------------------------------------------
        ++stepsSinceElmer;

        const scalar dUrel =
            gMax(mag(U.primitiveField() - U_old.primitiveField()))
           /max(gMax(mag(U_old.primitiveField())), SMALL);
        // Temperature: RMS relative change weighted by Joule heating power.
        // Conductivity only matters where current flows; the cold, current-free
        // layer at the walls changes every step as it develops and would
        // otherwise trigger an update on every step. Volume-weighted until
        // Elmer has produced any current.
        scalar dTrel = 0;
        {
            const scalarField relT
            (
                (T.primitiveField() - T_sent.primitiveField())
               /max(T_sent.primitiveField(), SMALL)
            );
            scalarField weight(max(JH.primitiveField(), scalar(0))*mesh.V().field());
            scalar weightSum = gSum(weight);
            if (weightSum <= VSMALL)
            {
                weight = mesh.V().field();
                weightSum = gSum(weight);
            }
            dTrel = Foam::sqrt(gSum(weight*sqr(relT))/max(weightSum, VSMALL));
        }
        const scalar dprel =
            gMax(mag(p.primitiveField() - p_sent.primitiveField()))
           /max(gMax(p_sent.primitiveField()), SMALL);

        const bool doElmer =
            status != 1
         || dUrel > couplingVelocityTolerance
         || dTrel > couplingTemperatureTolerance
         || dprel > couplingPressureTolerance
         || (couplingMaxSteps > 0 && stepsSinceElmer >= couplingMaxSteps);

        if (doElmer)
        {
            Info<< "Elmer update after " << stepsSinceElmer << " step(s); input change"
                << " U " << dUrel << "  T " << dTrel << "  p " << dprel << nl << endl;
            stepsSinceElmer = 0;
            ++nElmerUpdates;

            receiving.sendStatus(status);
            sending.sendStatus(status);

            // Construct scalar component fields from vectors
            Ux = U.component(vector::X);
            Uy = U.component(vector::Y);
            Uz = U.component(vector::Z);

            Bx = B.component(vector::X);
            By = B.component(vector::Y);
            Bz = B.component(vector::Z);

            // Inputs Elmer computes J from; U_old also sets the Lorentz damping
            U_old = U;
            T_sent = T;
            p_sent = p;
            sending.sendScalar(Ux);
            sending.sendScalar(Uy);
            sending.sendScalar(Uz);
            sending.sendScalar(Bx);
            sending.sendScalar(By);
            sending.sendScalar(Bz);
            sending.sendScalar(p);
            sending.sendScalar(T);

            receiving.recvScalar(Jx);
            receiving.recvScalar(Jy);
            receiving.recvScalar(Jz);
            receiving.recvScalar(JH_recv);
            receiving.recvScalar(potential);
            JH  = JH_recv;
            receiving.recvScalar(elcond_elmer);
            receiving.recvScalar(ionizationFraction);
            receiving.recvScalar(Te);
            Info<< "elcond min/max = " << gMin(elcond_elmer) << " " << gMax(elcond_elmer)
                << "  ionizationFraction max = " << gMax(ionizationFraction) << nl << endl;

            // Without this, parallel decomposition fails when writing fields with non-zero values
            Jx.correctBoundaryConditions();
            Jy.correctBoundaryConditions();
            Jz.correctBoundaryConditions();
            JH.correctBoundaryConditions();
            potential.correctBoundaryConditions();
            elcond_elmer.correctBoundaryConditions();
            ionizationFraction.correctBoundaryConditions();
            Te.correctBoundaryConditions();

            // Compute electric field E = -grad(potential)
            electric_field = -fvc::grad(potential);
            electric_field.correctBoundaryConditions();
        
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
            // Staleness of the force used since the previous update: the relative
            // (volume-weighted L2) jump in J x B when it is refreshed. This is the
            // error carried by skipping Elmer updates; keep it to a few percent.
            const vectorField fLorentzPrev(fLorentz.primitiveField());
            fLorentz = J_dens ^ B;
            lorentzDamping = elcond_elmer*magSqr(B);
            {
                const scalarField& V = mesh.V().field();
                const scalar forceChange = Foam::sqrt
                (
                    gSum(magSqr(fLorentz.primitiveField() - fLorentzPrev)*V)
                   /max(gSum(magSqr(fLorentz.primitiveField())*V), VSMALL)
                );

                // Mechanical power the flow loses to the field; matches Elmer's P_emf
                const scalar Pmech = -gSum((U.primitiveField() & fLorentz.primitiveField())*V);

                Info<< "Lorentz force: max |F| = " << gMax(mag(fLorentz.primitiveField()))
                    << " N/m^3  P_mech = " << Pmech << " W"
                    << "  change since last update " << forceChange
                    << "  (update " << nElmerUpdates << ")" << nl << endl;

                // dt / tau_mag with tau_mag = rho/(sigma B^2): above ~1 the lagged
                // force alone would be unstable and the implicit damping carries it
                Info<< "Lorentz damping: max dt/tau_mag = "
                    << gMax(lorentzDamping.primitiveField()/rho.primitiveField())*runTime.deltaTValue()
                    << nl << endl;
            }
        }

        // -----------------------------------------------------------------
        // PIMPLE loop
        // -----------------------------------------------------------------
        while (pimple.loop())
        {
            if (pimple.firstIter() && !pimple.simpleRho())
            {
                #include "rhoEqn.H"
            }

            #include "UEqn.H"
            #include "EEqn.H"

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

        rho = thermo.rho();

        {
            // Extremes of the compressible state, for monitoring
            Info<< "rho min/max = " << gMin(rho.primitiveField()) << " " << gMax(rho.primitiveField())
                << "  p min/max = " << gMin(p.primitiveField()) << " " << gMax(p.primitiveField())
                << "  T min/max = " << gMin(T.primitiveField()) << " " << gMax(T.primitiveField())
                << "  max Mach = "
                << gMax(mag(U.primitiveField())*sqrt(psi.primitiveField()/thermo.gamma()().primitiveField()))
                << nl << endl;
        }

        runTime.write();


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "Elmer updates: " << nElmerUpdates << nl << endl;
    Info<< "End\n" << endl;
    return 0;
}
