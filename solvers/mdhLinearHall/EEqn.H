{
    volScalarField& he = thermo.he();

    // Electromagnetic power delivered to the gas, J.E = J^2/sigma + U.(J x B).
    // The energy equation carries the kinetic energy K, so this is the total
    // work: Joule heating minus the mechanical power extracted by the generator
    // (U.(J x B) < 0). The force is the one the momentum equation applies,
    // including the Lorentz damping, so kinetic and total energy stay consistent.
    const volScalarField electromagneticPower
    (
        "electromagneticPower",
        JH + (U & (fLorentz - lorentzDamping*(U - U_old)))
    );

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        electromagneticPower
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
