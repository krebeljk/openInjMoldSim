# crossWLF

Compile the mojLibForces library before running this test. It is found in the solver directory - run Allwmake.

This is a virtual rheometer. The case simulates perfect shearing flow. The resulting wall shear stress is divided by the
shear rate to obtain the viscosity in the model. The result must match the prescribed Cross-WLF viscosity. This is
verified by running the analyze.py python script.
