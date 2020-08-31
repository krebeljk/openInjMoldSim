/*---------------------------------------------------------------------------*\
modified from: OpenFOAM-3.0.0/src/thermophysicalModels/specie/transport/sutherland/sutherlandTransport.C
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "polymerPVT.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::polymerPVT<Specie>::polymerPVT(Istream& is)
:
    Specie(is),
    b1m_(readScalar(is)),
    b2m_(readScalar(is)),
    b3m_(readScalar(is)),
    b4m_(readScalar(is)),
    b1s_(readScalar(is)),
    b2s_(readScalar(is)),
    b3s_(readScalar(is)),
    b4s_(readScalar(is)),
    b5_(readScalar(is)),
    b6_(readScalar(is)),
    b7_(readScalar(is)),
    b8_(readScalar(is)),
    b9_(readScalar(is))
{
    is.check("polymerPVT<Specie>::polymerPVT(Istream& is)");
}


template<class Specie>
Foam::polymerPVT<Specie>::polymerPVT(const dictionary& dict)
:
    Specie(dict),
    b1m_(readScalar(dict.subDict("equationOfState").lookup("b1m"))),
    b2m_(readScalar(dict.subDict("equationOfState").lookup("b2m"))),
    b3m_(readScalar(dict.subDict("equationOfState").lookup("b3m"))),
    b4m_(readScalar(dict.subDict("equationOfState").lookup("b4m"))),
    b1s_(readScalar(dict.subDict("equationOfState").lookup("b1s"))),
    b2s_(readScalar(dict.subDict("equationOfState").lookup("b2s"))),
    b3s_(readScalar(dict.subDict("equationOfState").lookup("b3s"))),
    b4s_(readScalar(dict.subDict("equationOfState").lookup("b4s"))),
    b5_(readScalar(dict.subDict("equationOfState").lookup("b5"))),
    b6_(readScalar(dict.subDict("equationOfState").lookup("b6"))),
    b7_(readScalar(dict.subDict("equationOfState").lookup("b7"))),
    b8_(readScalar(dict.subDict("equationOfState").lookup("b8"))),
    b9_(readScalar(dict.subDict("equationOfState").lookup("b9")))
{}


// * * * * * * * * * * * * * * * Podatki * * * * * * * * * * * * * //
template<class Specie>
Foam::scalar Foam::polymerPVT<Specie>::C_ = 0.0894; // to se pozneje posplosi

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::polymerPVT<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("b1m", b1m_);
    dict.add("b2m", b2m_);
    dict.add("b3m", b3m_);
    dict.add("b4m", b4m_);
    dict.add("b1s", b1s_);
    dict.add("b2s", b2s_);
    dict.add("b3s", b3s_);
    dict.add("b4s", b4s_);
    dict.add("b5", b5_);
    dict.add("b6", b6_);
    dict.add("b7", b7_);
    dict.add("b8", b8_);
    dict.add("b9", b9_);


    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const polymerPVT<Specie>& pf)
{
    os  << static_cast<const Specie&>(pf)
        << token::SPACE << pf.b1m_
        << token::SPACE << pf.b2m_
        << token::SPACE << pf.b3m_
        << token::SPACE << pf.b4m_
        << token::SPACE << pf.b1s_
        << token::SPACE << pf.b2s_
        << token::SPACE << pf.b3s_
        << token::SPACE << pf.b4s_
        << token::SPACE << pf.b5_
        << token::SPACE << pf.b6_
        << token::SPACE << pf.b7_
        << token::SPACE << pf.b8_
        << token::SPACE << pf.b9_;


    os.check("Ostream& operator<<(Ostream&, const polymerPVT<Specie>&)");
    return os;
}


// ************************************************************************* //
