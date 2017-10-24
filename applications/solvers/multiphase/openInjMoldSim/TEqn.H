{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - fvc::laplacian(mixture.kappa(), T)/mixture.Cp()
      - mixture.mu() * pow(strig,2)/mixture.Cp()//viscous heating
      - sqrt(T * mixture.psi() * (mixture.Cp() - mixture.Cv()) ) // <- T*alpha_v
        * (
            fvc::ddt(p)
          + fvc::div(phi, p)
          )// Dp/Dt
          /mixture.Cp()
    );

    TEqn.relax();
    TEqn.solve();

    mixture.correct();

    Info<< "min(T) " << min(T).value() << endl;
}
