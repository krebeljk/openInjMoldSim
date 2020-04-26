# openInjMoldDyMSimAmClr v1.0

This is a fork of [openInjMoldSim](https://github.com/krebeljk/openInjMoldSim) for amorphous polymer volumetric relaxation modeling.  

This fork and the cases are described in a [paper](https://rdcu.be/but3t).

This code is not part of OpenFOAM. It only uses OpenFOAM v3.0.1. To learn about OpenFOAM see [openfoam.org](https://openfoam.org/), unofficial wiki
[openfoamwiki.net](https://openfoamwiki.net/index.php/Main_Page) and the community forum
[cfd-online.com](https://www.cfd-online.com/Forums/openfoam/).

## Getting Started

An OpenFOAM 3.0.1 installation is required to run the code.

Run `wmake` in the [`openInjMoldDyMSim`](applications/solvers/multiphase/openInjMoldSim/openInjMoldDyMSim) directory.

Run `Allwmake` in [`externalWallHeatFluxTemperatureP`](/tutorials/externalWallHeatFluxTemperatureP) tutorial to include the pressure dependent thermal boundary condition.  

Run `AllRun` in the [`tutorials`](/tutorials/) directory to run the example simulation. Then run `paraFoam` in the case directory to view the results.

## Functionality

* Compressible, non-isothermal, laminar cavity flow.
* Tait equation of state with viscoelastic volume relaxation.
* Cross-WLF viscosity model.
* Elastic solidification.
* Cases model a compliant mold wall implemented through dynamic mesh and a coded boundary condition.

## Feedback

Any feedback is apreciated - krebeljk()gmail.com.

## Authors

* **Kristjan Krebelj** - Assembled and tested the code. Maintainer of the repository.
* **Janez Turk** - Introduced the key modifications to the original OpenFOAM library.

## License

This project is licensed under the GPU License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* The work was supported by the [Laboratory for Numerical Modelling and Simulation - LNMS](http://lab.fs.uni-lj.si/lnms/).
