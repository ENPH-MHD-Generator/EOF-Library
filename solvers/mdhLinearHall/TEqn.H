{
    // Effective thermal diffusivity: conductive + turbulent part
    volScalarField kappaf
    (
        IOobject
        (
            "kappaf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        kappa/cp + rho*turbulence->nut()/0.7
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvm::laplacian(kappaf, T)
     ==
        JH/cp
    );

    TEqn.relax();
    TEqn.solve();
}
