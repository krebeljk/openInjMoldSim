/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    compressibleInterFoam

Description
    Solver for 2 compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e.  laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "mojRhoThermo.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "mojTwoPhaseMixtureThermo.H"
#include "mojTurbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info << "openInjMoldSim v7.1 (fiber)" << endl;
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldsFiber.H" //NEW - Kerstin
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "createTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaEqnsSubCycle.H"

            #include "FiberClosure.H" //NEW - Kerstin
            #include "A2Eqn.H" //NEW - Kerstin

            // correct interface on first PIMPLE corrector
           // if (pimple.corr() == 1)
            if (pimple.firstPimpleIter())  //cbpf on 24/09/2019 compatibility to of6
            {
                interface.correct();
            }

            solve(fvm::ddt(rho) + fvc::div(rhoPhi));

            //update state
            shrRate = sqrt(2.0*symm(fvc::grad(U))&&symm(fvc::grad(U)));
            mixture.correct();
            visc = mixture.mu();
            mojKappaOut = mixture.kappa();

            //Kristjan: Elastic deviatoric stress equation
            if (sldDictIO.typeHeaderOk<IOdictionary>()) //VER SE ESTÁ BEM
            {
                fvSymmTensorMatrix elSigDevEqn(
                  fvm::ddt(elSigDev)
                + fvm::div(phi,elSigDev)
                + fvm::SuSp(-fvc::div(phi),elSigDev)
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

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }


            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        Tc = T - twoSevenThree;
        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
