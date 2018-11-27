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

    volScalarField vfeq = vf_; // allocate equil free vol
    volScalarField vg = vf_; // allocate glassy vol
    dimensionedScalar tauMin = scalar(10.0) * this->T_.mesh().time().deltaT();// min relaxation time
    volScalarField tauRlx
    (
        IOobject
        (
            "tauRlx",
            this->T_.mesh().time().timeName(),
            this->T_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->T_.mesh(),
        dimensionedScalar("tauRlxInit", dimensionSet(0,0,1,0,0,0,0), 1.0)
    );

    scalarField& TCells = this->T_.internalField();
    scalarField& psiCells = this->psi_.internalField();
    scalarField& vfeqCells = vfeq.internalField();
    scalarField& vgCells = vg.internalField();
    scalarField& muCells = this->mu_.internalField();
    scalarField& alphaCells = this->alpha_.internalField();
    scalarField& sfTauRlx = tauRlx.internalField();

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
        vfeqCells[celli] = mixture_.vfeq(pCells[celli], TCells[celli]);
        vgCells[celli] = mixture_.vg(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli], strigCells[celli]);

        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
        sfTauRlx[celli] = mixture_.tauRlx(pCells[celli], TCells[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pstrig = this->strig_.boundaryField()[patchi];
	//fvPatchScalarField& pU = strig.boundaryField()[patchi]; // tukaj boundary

        fvPatchScalarField& ppsi = this->psi_.boundaryField()[patchi];
        fvPatchScalarField& pvfeq = vfeq.boundaryField()[patchi];
        fvPatchScalarField& pvg = vg.boundaryField()[patchi];

        //fvPatchScalarField& ph = this->he().boundaryField()[patchi];//Kristjan: governed by TEqn


        fvPatchScalarField& pmu = this->mu_.boundaryField()[patchi];
        fvPatchScalarField& palpha = this->alpha_.boundaryField()[patchi];
        fvPatchScalarField& ptauRlx = tauRlx.boundaryField()[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                //ph[facei] = mixture_.HE(pp[facei], pT[facei]);//Kristjan: governed by TEqn


                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                pvfeq[facei] = mixture_.vfeq(pp[facei], pT[facei]);
                pvg[facei] = mixture_.vg(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei], pstrig[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                ptauRlx[facei] = mixture_.tauRlx(pp[facei], pT[facei]);
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
                pvfeq[facei] = mixture_.vfeq(pp[facei], pT[facei]);
                pvg[facei] = mixture_.vg(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei], pstrig[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
                ptauRlx[facei] = mixture_.tauRlx(pp[facei], pT[facei]);
            }
        }
    }
    // calculate lagging density

    // relaxation time

    const surfaceScalarField& phi = this->db().objectRegistry::lookupObject<surfaceScalarField>("phi");

    // volume relaxation
    fvScalarMatrix vfEqn
    (
        //(
        fvm::ddt(vf_)
      + fvm::div(phi, vf_)
      + fvm::SuSp(-fvc::div(phi), vf_)
        //) ... Dvf_/Dt
        ==
        - (vf_ - vfeq)/max(tauRlx, tauMin)
     );
    vfEqn.relax();
    vfEqn.solve();

    // get current density
    this->rho_ = scalar(1)/(vg + vf_);
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("vZero", dimensionSet(-1,3,0,0,0,0,0), 0.0),
        "zeroGradient"
    )
{
    if(max(vf_).value() < VSMALL)//not a restart -> init vf
    {
        const scalarField& pCells = this->p_.internalField();
        scalarField& TCells = this->T_.internalField();
        scalarField& vfCells = this->vf_.internalField();

        forAll(TCells, celli)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->cellMixture(celli);

            vfCells[celli] = mixture_.vfeq(pCells[celli], TCells[celli]);
        }

        forAll(this->T_.boundaryField(), patchi)
        {
            fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
            fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
            fvPatchScalarField& pvf = this->vf_.boundaryField()[patchi];

            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pvf[facei] = mixture_.vfeq(pp[facei], pT[facei]);
            }
        }
    }
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
