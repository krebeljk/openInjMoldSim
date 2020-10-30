    if(FiberSuspension) //NEW - Kerstin
    {
        StressFiber = 2.0*visc*FiberVolfrac*symm(mu1*D + mu2*ClosureA4);
    }
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + turbulence->divDevRhoReff(U)
      - fvc::div(elSigDev)
      - fvc::div(StressFiber)
      + fvc::grad(pAux)//Kristjan: auxiliary pressure for pMin equilibrium
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        K = 0.5*magSqr(U);
    }
