# openInjMoldSim [`v1.1.1`](VERSION.md)

This is an [OpenFOAM 3.0.1](http://openfoamwiki.net/index.php/Installation) solver for simulation of injection molding filling, packing and cooling stages.
It is a modification of the compressibleInterFoam solver distributed with OpenFOAM.
The solver has so far been validated on the [2D demo](https://nbviewer.jupyter.org/github/krebeljk/openInjMoldSim/blob/master/tutorials/Tutorials.ipynb) case for amorphous polystyrene.
A simulation of (semi-crystalline) HDPE is under development.
Currently, the simulation is prohibitively inconvenient for typical industrial use.

See the [tutorials](https://nbviewer.jupyter.org/github/krebeljk/openInjMoldSim/blob/master/tutorials/Tutorials.ipynb) for more info.

See the official website [openfoam.org](https://openfoam.org/), unofficial wiki
[openfoamwiki.net](https://openfoamwiki.net/index.php/Main_Page) and the community forum
[cfd-online.com](https://www.cfd-online.com/Forums/openfoam/).

![openInjMoldSimAxiFill](openInjMoldSimAxiFill.gif)

## Getting Started

An [OpenFOAM 3.0.1 installation](http://openfoamwiki.net/index.php/Installation) is required to run the code. This code was run on Linux Ubuntu 14.04 LTS and Centos 6.6 (Rocks 6.2 Sidewinder).

Run `Allwmake` in the solver directory to compile the solver.

Run `AllRun` in the [`fill_pack`](/tutorials/demo/fill_pack) directory to run the example simulation. Then run `paraFoam` in the case directory to view the results.

## Functionality

* Compressible, non-isothermal, laminar cavity flow.
* Tait equation of state.
* Cross-WLF viscosity model.
* Specific heat and thermal conductivity may depend on temperature and pressure (tabular form).
* Elastic stress in the solid phase.

## Feedback

Any feedback is apreciated - krebeljk()gmail.com.

## Authors

* **Kristjan Krebelj** - Assembled and tested the code. Maintainer of the repository.
* **Janez Turk** - Introduced the key modifications to the original OpenFOAM library.

## License

This project is licensed under the GPU License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* The work was supported by the [Laboratory for Numerical Modelling and Simulation - LNMS](http://lab.fs.uni-lj.si/lnms/).

## Demo case
![Demo geometry](/tutorials/demo/demo_geom.png)
![Pressure evolution](/tutorials/demo/plot.png)
