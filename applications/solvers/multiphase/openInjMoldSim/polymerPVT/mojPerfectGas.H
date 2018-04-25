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

Class
    Foam::mojPerfectGas

Description
    Perfect gas equation of state.

SourceFiles
    mojPerfectGasI.H
    mojPerfectGas.C

\*---------------------------------------------------------------------------*/

#ifndef mojPerfectGas_H
#define mojPerfectGas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class mojPerfectGas;

template<class Specie>
inline mojPerfectGas<Specie> operator+
(
    const mojPerfectGas<Specie>&,
    const mojPerfectGas<Specie>&
);

template<class Specie>
inline mojPerfectGas<Specie> operator-
(
    const mojPerfectGas<Specie>&,
    const mojPerfectGas<Specie>&
);

template<class Specie>
inline mojPerfectGas<Specie> operator*
(
    const scalar,
    const mojPerfectGas<Specie>&
);

template<class Specie>
inline mojPerfectGas<Specie> operator==
(
    const mojPerfectGas<Specie>&,
    const mojPerfectGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const mojPerfectGas<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class mojPerfectGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class mojPerfectGas
:
    public Specie
{

public:

    // Constructors

        //- Construct from components
        inline mojPerfectGas(const Specie& sp);

        //- Construct from Istream
        mojPerfectGas(Istream&);

        //- Construct from dictionary
        mojPerfectGas(const dictionary& dict);

        //- Construct as named copy
        inline mojPerfectGas(const word& name, const mojPerfectGas&);

        //- Construct and return a clone
        inline autoPtr<mojPerfectGas> clone() const;

        // Selector from Istream
        inline static autoPtr<mojPerfectGas> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<mojPerfectGas> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "mojPerfectGas<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const mojPerfectGas&);
        inline void operator-=(const mojPerfectGas&);

        inline void operator*=(const scalar);


    // Friend operators

        friend mojPerfectGas operator+ <Specie>
        (
            const mojPerfectGas&,
            const mojPerfectGas&
        );

        friend mojPerfectGas operator- <Specie>
        (
            const mojPerfectGas&,
            const mojPerfectGas&
        );

        friend mojPerfectGas operator* <Specie>
        (
            const scalar s,
            const mojPerfectGas&
        );

        friend mojPerfectGas operator== <Specie>
        (
            const mojPerfectGas&,
            const mojPerfectGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const mojPerfectGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mojPerfectGasI.H"

#ifdef NoRepository
#   include "mojPerfectGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
