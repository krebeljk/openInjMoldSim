# shearModulus

This case is similar to the crossWLF test but tests the shear mosulus. The Allrun script first runs a virtual rheometer
to establish shearing flow without elasticity.
Begining with the last frame of the first simulation, a second simulation starts which includes elasticity. The shearing
deformation causes elSigDevxy shear stress. The elSigDevxy is written in the second column in the results file or it can
be viewed with the paraFoam viewer. Shear modulus is defined in the constant/solidificationProperties dictionary. It is
first disabled with by setting the shearRate limit in the negative. The Allrun script than changes the limit to a very
high number to include the elasticity.

elSigDevxy = shearModulus * gamma_xy
           = shearModulus * deltaY/L_x
           = shearModulus * v_y/L_x * time

gamma_xy   ... engineering shear deformation (2*epsilon_xy)
deltaY     ... displacement of the far right side of the domain
L_x        ... wall thicness (see the constant/polyMesh/blockMeshDict file)
v_y        ... right side wall velocity (see the U file)
time       ... duration of the simulation with elasticiy

