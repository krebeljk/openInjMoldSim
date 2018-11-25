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

#include "taitTg.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::taitTg<Specie>::taitTg(Istream& is)
:
    Specie(is),
    D2_(readScalar(is)),
    D3_(readScalar(is)),
    A1_(readScalar(is)),
    A2_(readScalar(is)),
    tauRlx_(readScalar(is)),
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
    b10_(readScalar(is))
{
    is.check("taitTg<Specie>::taitTg(Istream& is)");
}


template<class Specie>
Foam::taitTg<Specie>::taitTg(const dictionary& dict)
:
    Specie(dict),
    D2_(readScalar(dict.subDict("equationOfState").lookup("D2"))),
    D3_(readScalar(dict.subDict("equationOfState").lookup("D3"))),
    A1_(readScalar(dict.subDict("equationOfState").lookup("A1"))),
    A2_(readScalar(dict.subDict("equationOfState").lookup("A2"))),
    tauRlx_(readScalar(dict.subDict("equationOfState").lookup("tauRlx"))),
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
    b10_(readScalar(dict.subDict("equationOfState").lookup("b10")))
{
    Info << "TaitTg:"    << endl;
    Info <<   "D2       " << D2_
         << "\nD3       " << D3_
         << "\nA1       " << A1_
         << "\nA2       " << A2_
         << "\ntauRlx   " << tauRlx_
         << "\nb1m      " << b1m_
         << "\nb2m      " << b2m_
         << "\nb3m      " << b3m_
         << "\nb4m      " << b4m_
         << "\nb1s      " << b1s_
         << "\nb2s      " << b2s_
         << "\nb3s      " << b3s_
         << "\nb4s      " << b4s_
         << "\nb5       " << b5_
         << "\nb6       " << b6_
         << "\nb10;     " << b10_
         << endl;
}


// * * * * * * * * * * * * * * * Podatki * * * * * * * * * * * * * //
template<class Specie>
Foam::scalar Foam::taitTg<Specie>::C_ = 0.0894; // to se pozneje posplosi

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::taitTg<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("D2", D2_);
    dict.add("D3", D3_);
    dict.add("A1", A1_);
    dict.add("A2", A2_);
    dict.add("tauRlx", tauRlx_);
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
    dict.add("b10", b10_);


    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const taitTg<Specie>& pf)
{
    os  << static_cast<const Specie&>(pf)
        << token::SPACE << pf.D2_
        << token::SPACE << pf.D3_
        << token::SPACE << pf.A1_
        << token::SPACE << pf.A2_
        << token::SPACE << pf.tauRlx_
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
        << token::SPACE << pf.b10_;


    os.check("Ostream& operator<<(Ostream&, const taitTg<Specie>&)");
    return os;
}


// ************************************************************************* //
