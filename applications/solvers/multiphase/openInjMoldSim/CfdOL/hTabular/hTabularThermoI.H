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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const EquationOfState& pt,
    const scalar Hf
)
:
    EquationOfState(pt),
    Hf_(Hf)
{
    cpTable = interpolation2DTable<scalar>("constant/cpTable");
    //hTable = interpolation2DTable<scalar>("constant/hTable");
    cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    //hTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const hTabularThermo& pt
)
:
    EquationOfState(pt),
    Hf_(pt.Hf_)
{
    cpTable = interpolation2DTable<scalar>("constant/cpTable");
    //hTable = interpolation2DTable<scalar>("constant/hTable");
    cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    //hTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>::hTabularThermo
(
    const word& name,
    const hTabularThermo& pt
)
:
    EquationOfState(name, pt),
    Hf_(pt.Hf_)
{
    cpTable = interpolation2DTable<scalar>("constant/cpTable");
    //hTable = interpolation2DTable<scalar>("constant/hTable");
    cpTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
    //hTable.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::limit
(
    const scalar T
) const
{
  return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::cp
(
    const scalar p, const scalar T
) const
{
    //Info << "cp for T= " << T << ", p= " << p  << ": " << cpTable(T, p) * this->W() << endl;
    return cpTable(T, p) * this->W();
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::ha
(
    const scalar p, const scalar T
) const
{
    //Info << "ha for T= " << T << ": " << ((hTable(T, p) + hc()) * this->W()) << endl;
    //return (hTable(T, p) + hc()) * this->W();
    return hc() * this->W();
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::hs
(
    const scalar p, const scalar T
) const
{
    //Info << "hs for T= " << T << ", p= " << p << ": " << (hTable(T, p) * this->W()) << endl;
    //abort(FatalError);
    //return (hTable(T, p)) * this->W();
    return hc()* this->W();
}


template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::hc()
const
{
    return Hf_;
}

template<class EquationOfState>
inline Foam::scalar Foam::hTabularThermo<EquationOfState>::s
(
    const scalar p, const scalar T
) const
{
    notImplemented
    (
        "scalar hTabularThermo<EquationOfState>::"
        "s(const scalar p, const scalar T) const"
    );
    return T;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState>&
Foam::hTabularThermo<EquationOfState>::operator=
(
    const hTabularThermo<EquationOfState>& pt
)
{
    EquationOfState::operator=(pt);

    Hf_ = pt.Hf_;

    return *this;
}


template<class EquationOfState>
inline void Foam::hTabularThermo<EquationOfState>::operator+=
(
    const hTabularThermo<EquationOfState>& pt
)
{
}


template<class EquationOfState>
inline void Foam::hTabularThermo<EquationOfState>::operator-=
(
    const hTabularThermo<EquationOfState>& pt
)
{
}


template<class EquationOfState>
inline void Foam::hTabularThermo<EquationOfState>::operator*=
(
    const scalar s
)
{
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState> Foam::operator+
(
    const hTabularThermo<EquationOfState>& pt1,
    const hTabularThermo<EquationOfState>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs += pt2;

    scalar molr1 = pt1.nMoles()/eofs.nMoles();
    scalar molr2 = pt2.nMoles()/eofs.nMoles();

    return hTabularThermo<EquationOfState>
    (
        eofs,
        molr1*pt1.Hf_ + molr2*pt2.Hf_
    );
}


template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState> Foam::operator-
(
    const hTabularThermo<EquationOfState>& pt1,
    const hTabularThermo<EquationOfState>& pt2
)
{
    EquationOfState eofs = pt1;
    eofs -= pt2;

    scalar molr1 = pt1.nMoles()/eofs.nMoles();
    scalar molr2 = pt2.nMoles()/eofs.nMoles();

    return hTabularThermo<EquationOfState>
    (
        eofs,
        molr1*pt1.Hf_ - molr2*pt2.Hf_
    );
}


template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hTabularThermo<EquationOfState>& pt
)
{
    return hTabularThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.Hf_
    );
}


template<class EquationOfState>
inline Foam::hTabularThermo<EquationOfState> Foam::operator==
(
    const hTabularThermo<EquationOfState>& pt1,
    const hTabularThermo<EquationOfState>& pt2
)
{
    return pt2 - pt1;
}


// ************************************************************************* //
