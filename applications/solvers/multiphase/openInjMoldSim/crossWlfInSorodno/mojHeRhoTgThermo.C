/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "mojHeRhoTgThermo.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::mojHeRhoTgThermo<BasicPsiThermo, MixtureType>::calculate()
{
    //const scalarField& hCells = this->he().internalField(); //Kristjan: governed by TEqn
    const scalarField& pCells = this->p_.internalField();
    const scalarField& strigCells = this->strig_.internalField();

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& rhoCells = this->rho_.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();
    scalarField& vfCells = this->vf_.internalField();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        //TCells[celli] = mixture_.THE //Kristjan: governed by TEqn
        //(
        //    hCells[celli],
        //    pCells[celli],
        //    TCells[celli]
        //);

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli], strigCells[celli]);

        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
        vfCells[celli] =   mixture_.veq(pCells[celli], TCells[celli])
                         - mixture_.vg(pCells[celli], TCells[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pstrig = this->strig_.boundaryField()[patchi];
	//fvPatchScalarField& pU = strig.boundaryField()[patchi]; // tukaj boundary

        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];

        //fvPatchScalarField& ph = this->he().boundaryField()[patchi];//Kristjan: governed by TEqn


        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];
        fvPatchScalarField& pvf = this->vf_.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                //ph[facei] = mixture_.HE(pp[facei], pT[facei]);//Kristjan: governed by TEqn


                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei], pstrig[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvf[facei] =    mixture_.veq(pp[facei], pT[facei])
                              - mixture_.vg(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                //pT[facei] = mixture_.THE(ph[facei], pp[facei], pT[facei]);//Kristjan: governed by TEqn

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei], pstrig[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                pvf[facei] =    mixture_.veq(pp[facei], pT[facei])
                              - mixture_.vg(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::mojHeRhoTgThermo<BasicPsiThermo, MixtureType>::mojHeRhoTgThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    mojHeThermo<BasicPsiThermo, MixtureType>(mesh, phaseName),
    vf_
    (
        IOobject
        (
            "vf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(-1, 3, 0, 0, 0)
    )
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::mojHeRhoTgThermo<BasicPsiThermo, MixtureType>::~mojHeRhoTgThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::mojHeRhoTgThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering mojHeRhoTgThermo<MixtureType>::correct()" << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "exiting mojHeRhoTgThermo<MixtureType>::correct()" << endl;
    }
}


// ************************************************************************* //
