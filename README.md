# openInjMoldSim [`v1.1.1`](VERSION.md)

This is an OpenFOAM-3.0.1 solver for simulation of injection molding filling, packing and cooling stages. It is a
modification of the compressibleInterFoam solver distributed with OpenFOAM. The solver has so far been used on the small
2D demo case for amorphous polystyrene. A simulation of (semi-crystalline) HDPE is under development. Currently,
the simulation is prohibitively inconvenient for typical industrial use.

See the official website [openfoam.org](https://openfoam.org/), unofficial wiki
[openfoamwiki.net](https://openfoamwiki.net/index.php/Main_Page) and the community forum
[cfd-online.com](https://www.cfd-online.com/Forums/openfoam/).

## Getting Started

An OpenFOAM 3.0.1 installation is required to run the code. This code was run on Linux Ubuntu 14.04 LTS and Centos 6.6 (Rocks 6.2 Sidewinder).

Run `Allwmake` in the solver directory to compile the solver.

Run `AllRun` in the [`fill_pack`](/tutorials/demo/fill_pack) directory to run the example simulation. Then run `paraFoam` in the case directory to view the results.

## Functionality

* Compressible, non-isothermal, laminar cavity flow.
* Tait equation of state.
* Cross-WLF viscosity model.

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
