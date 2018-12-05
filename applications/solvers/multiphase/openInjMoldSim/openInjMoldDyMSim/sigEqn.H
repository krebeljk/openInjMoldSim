if (sldDictIO.headerOk())
{
    fvSymmTensorMatrix elSigDevEqn(
      fvm::ddt(elSigDev)
    + fvm::div(phi,elSigDev)
    + fvm::SuSp(-fvc::div(phi),elSigDev)
      ==
      dev(
          twoSymm(elSigDev & fvc::grad(U))
        + shrMod * twoSymm(fvc::grad(U))
         )
      * pos(shrRateLimEl-shrRate)
      * pos(visc-viscLimEl)
    );
    elSigDevEqn.relax();
    elSigDevEqn.solve();
    elSigDev = dev(elSigDev);
}
