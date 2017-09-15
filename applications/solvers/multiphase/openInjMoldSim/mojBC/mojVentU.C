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

#include "mojVentU.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mojVentU
::mojVentU
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF),
    alphaLimit_(0),
    alphaName_("none")
{
    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


Foam::mojVentU
::mojVentU
(
    const mojVentU& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    alphaLimit_(ptf.alphaLimit_),
    alphaName_(ptf.alphaName_)
{}


Foam::mojVentU
::mojVentU
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    alphaLimit_(readScalar(dict.lookup("alphaLimit"))),
    alphaName_(dict.lookup("alpha"))
{
    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchVectorField::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}


Foam::mojVentU
::mojVentU
(
    const mojVentU& ptf
)
:
    mixedFvPatchField<vector>(ptf),
    alphaLimit_(ptf.alphaLimit_),
    alphaName_(ptf.alphaName_)
{}


Foam::mojVentU
::mojVentU
(
    const mojVentU& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    alphaLimit_(ptf.alphaLimit_),
    alphaName_(ptf.alphaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mojVentU::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField alphap =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);

    alphap = max(alphap, scalar(0));
    alphap = min(alphap, scalar(1));

    valueFraction() = pos(alphap - alphaLimit_);
    //packing: alpha > alphaLim => valFrac=1 (U=0)
    //else filling: valFrac=0 (gradU=0)

    mixedFvPatchField<vector>::updateCoeffs();
}


void Foam::mojVentU::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

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
       fvPatchVectorField,
       mojVentU
   );
}


// ************************************************************************* //
