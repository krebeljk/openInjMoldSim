if (sldDictIO.headerOk())
{
    fvSymmTensorMatrix elSigDevEqn(
      fvm::ddt(elSigDev)
    + fvm::div(phi,elSigDev)
    - fvc::div(phi)*elSigDev
      ==
      twoSymm(elSigDev & fvc::grad(U))
    + shrMod * dev(twoSymm(fvc::grad(U)))
      * pos(shrRateLimEl-shrRate)
      * pos(visc-viscLimEl)
    );
    elSigDevEqn.relax();
    elSigDevEqn.solve();
    elSigDev = dev(elSigDev);
}
