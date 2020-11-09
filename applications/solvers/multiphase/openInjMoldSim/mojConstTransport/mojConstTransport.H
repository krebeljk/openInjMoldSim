/*---------------------------------------------------------------------------*\
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

Class
    Foam::mojConstTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    mojConstTransportI.H
    mojConstTransport.C

\*---------------------------------------------------------------------------*/

#ifndef mojConstTransport_H
#define mojConstTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class mojConstTransport;

template<class Thermo>
inline mojConstTransport<Thermo> operator+
(
    const mojConstTransport<Thermo>&,
    const mojConstTransport<Thermo>&
);

template<class Thermo>
inline mojConstTransport<Thermo> operator-
(
    const mojConstTransport<Thermo>&,
    const mojConstTransport<Thermo>&
);

template<class Thermo>
inline mojConstTransport<Thermo> operator*
(
    const scalar,
    const mojConstTransport<Thermo>&
);

template<class Thermo>
inline mojConstTransport<Thermo> operator==
(
    const mojConstTransport<Thermo>&,
    const mojConstTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const mojConstTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class mojConstTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class mojConstTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline mojConstTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline mojConstTransport(const word&, const mojConstTransport&);

        //- Construct from Istream
        mojConstTransport(Istream&);

        //- Construct from dictionary
        mojConstTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<mojConstTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<mojConstTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<mojConstTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "mojConst<" + Thermo::typeName() + '>'; //to definira vnos v dictionary!
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T, const scalar shrRate) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline mojConstTransport& operator=(const mojConstTransport&);

        inline void operator+=(const mojConstTransport&);

        inline void operator-=(const mojConstTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend mojConstTransport operator+ <Thermo>
        (
            const mojConstTransport&,
            const mojConstTransport&
        );

        friend mojConstTransport operator- <Thermo>
        (
            const mojConstTransport&,
            const mojConstTransport&
        );

        friend mojConstTransport operator* <Thermo>
        (
            const scalar,
            const mojConstTransport&
        );

        friend mojConstTransport operator== <Thermo>
        (
            const mojConstTransport&,
            const mojConstTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const mojConstTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mojConstTransportI.H"

#ifdef NoRepository
#   include "mojConstTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
