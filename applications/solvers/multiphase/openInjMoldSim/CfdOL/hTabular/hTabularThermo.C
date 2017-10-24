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

#include "hTabularThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    Istream& is
)
:
    EquationOfState(is),
    Hf_(readScalar(is))
{
    Hf_ *= this->W();
    cpTable = interpolation2DTable<scalar>("constant/cpTable");
    //hTable = interpolation2DTable<scalar>("constant/hTable");
    cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    //hTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
/*
    hCoeffs_ = CpCoeffs_.integral();
    sCoeffs_ = CpCoeffs_.integralMinus1();

    // Offset h poly so that it is relative to the enthalpy at Tstd
    hCoeffs_[0] += Hf_ - hCoeffs_.value(specie::Tstd);

    // Offset s poly so that it is relative to the entropy at Tstd
    sCoeffs_[0] += Sf_ - sCoeffs_.value(specie::Tstd);
    */
}


template<class EquationOfState>
Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf")))
{
    Hf_ *= this->W();
    cpTable = interpolation2DTable<scalar>("constant/cpTable");
    //hTable = interpolation2DTable<scalar>("constant/hTable");
    cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    //hTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hTabularThermo<EquationOfState>::write
(
    Ostream& os
) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Hf", Hf_/this->W());
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hTabularThermo<EquationOfState>& pt
)
{
    os  << static_cast<const EquationOfState&>(pt) << tab
        << pt.Hf_/pt.W() << tab;

    os.check
    (
        "operator<<"
        "("
            "Ostream&, "
            "const hTabularThermo<EquationOfState>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
