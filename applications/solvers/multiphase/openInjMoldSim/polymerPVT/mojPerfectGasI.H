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

\*---------------------------------------------------------------------------*/

#include "mojPerfectGas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::mojPerfectGas<Specie>::mojPerfectGas(const Specie& sp)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::mojPerfectGas<Specie>::mojPerfectGas
(
    const word& name,
    const mojPerfectGas<Specie>& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::mojPerfectGas<Specie> >
Foam::mojPerfectGas<Specie>::clone() const
{
    return autoPtr<mojPerfectGas<Specie> >(new mojPerfectGas<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::mojPerfectGas<Specie> >
Foam::mojPerfectGas<Specie>::New(Istream& is)
{
    return autoPtr<mojPerfectGas<Specie> >(new mojPerfectGas<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::mojPerfectGas<Specie> > Foam::mojPerfectGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<mojPerfectGas<Specie> >(new mojPerfectGas<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::mojPerfectGas<Specie>::rho(scalar p, scalar T) const
{
    return p/(this->R()*T);
}

template<class Specie>
inline Foam::scalar Foam::mojPerfectGas<Specie>::cTExp(scalar p, scalar T) const
{
    return 1/T;
}

template<class Specie>
inline Foam::scalar Foam::mojPerfectGas<Specie>::s(scalar p, scalar T) const
{
    return -RR*log(p/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::mojPerfectGas<Specie>::psi(scalar p, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::mojPerfectGas<Specie>::Z(scalar p, scalar T) const
{
    return 1;
}


template<class Specie>
inline Foam::scalar Foam::mojPerfectGas<Specie>::cpMcv(scalar p, scalar T) const
{
    return RR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::mojPerfectGas<Specie>::operator+=(const mojPerfectGas<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::mojPerfectGas<Specie>::operator-=(const mojPerfectGas<Specie>& pg)
{
    Specie::operator-=(pg);
}


template<class Specie>
inline void Foam::mojPerfectGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::mojPerfectGas<Specie> Foam::operator+
(
    const mojPerfectGas<Specie>& pg1,
    const mojPerfectGas<Specie>& pg2
)
{
    return mojPerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::mojPerfectGas<Specie> Foam::operator-
(
    const mojPerfectGas<Specie>& pg1,
    const mojPerfectGas<Specie>& pg2
)
{
    return mojPerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::mojPerfectGas<Specie> Foam::operator*
(
    const scalar s,
    const mojPerfectGas<Specie>& pg
)
{
    return mojPerfectGas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::mojPerfectGas<Specie> Foam::operator==
(
    const mojPerfectGas<Specie>& pg1,
    const mojPerfectGas<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
