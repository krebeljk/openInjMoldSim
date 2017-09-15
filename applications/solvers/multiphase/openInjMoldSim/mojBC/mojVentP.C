/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "mojVentP.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mojVentP
::mojVentP
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    alphaLimit_(0),
    alphaName_("none")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::mojVentP
::mojVentP
(
    const mojVentP& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    alphaLimit_(ptf.alphaLimit_),
    alphaName_(ptf.alphaName_)
{}


Foam::mojVentP
::mojVentP
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    alphaLimit_(readScalar(dict.lookup("alphaLimit"))),
    alphaName_(dict.lookup("alpha"))
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(patchInternalField());
    }
}


Foam::mojVentP
::mojVentP
(
    const mojVentP& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    alphaLimit_(ptf.alphaLimit_),
    alphaName_(ptf.alphaName_)
{}


Foam::mojVentP
::mojVentP
(
    const mojVentP& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    alphaLimit_(ptf.alphaLimit_),
    alphaName_(ptf.alphaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mojVentP::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField alphap =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);

    alphap = max(alphap, scalar(0));
    alphap = min(alphap, scalar(1));

    valueFraction() = pos(alphaLimit_ - alphap);
    //alphLim>alpha => valFrac=1 (U=0)
    //else valFrac=0 (gradU=0)

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::mojVentP::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("alphaLimit") << alphaLimit_
        << token::END_STATEMENT << nl;
    os.writeKeyword("alpha") << alphaName_
        << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       mojVentP
   );
}


// ************************************************************************* //
