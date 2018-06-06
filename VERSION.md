# openInjMoldSim v1.1 (2.11.2017)
Version v1.1 introduces:
- The proper energy equation.
- Tabulated specific heat constant/cpTable.
- Tabulated thermal conductivity constant/kappaTable.
- Elastic solidification defined by constant/solidificationProperties.

## Commentary
Input specific heat J/kgK (kappaTable) at constant pressure for a list of temperatures in K.
Input thermal conductivity in W/mK (kappaTable) at temperatures in K and pressures in Pa.
Values outside of the provided range will be taken from the table edge.

Elastic shear modulus will be effected when viscosity surpasses the limit specified in the solidificationProperties
dictionary. To control this according to temperature, change the TnoFlow temperature of the CrossWLF model. According to
experience, elastic behavior should be limited to slow flow. This is effected with the shearRate limit, also given in the
solidificationProperties. Not providing the solidificationProperties dictionary will make the solver skip the elastic
calculation and only use the CrossWLF viscosity. The fvSchemes for elasticity must still be provided even though they are not used.

## Version v1.1.1
- Important bug fix for thermal conductivity.
- Parameter deltaTempInterp in crossWLF.
- Renamed lambda to kappa.

The [bug fix for kappa](https://bugs.openfoam.org/view.php?id=2532) was incorporated. This was needed for
semi-crystalline simulation where Cv is significantly less than Cp (specific heat at constant pressure).

The parameter deltaTempInterp can now be chosen or left at 5K. It allows to make a gradual transition to the maximal
viscosity on solidification.

The naming of constant thermal conductivity as lambda was changed to kappa.

## Version v1.1.2
- Bug fix. Thermal conductivity in the energy equation calculated properly.
- test/CrossWLF tests multiple temperatures and numerical parameters are corrected (see jupyter notebook).
