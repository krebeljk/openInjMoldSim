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

Class
    Foam::hTabularThermo

Description
    Thermodynamics package templated on the equation of state, using tabular
    functions for cp, h and s

    h and s derived from cp

SourceFiles
    hTabularThermoI.H
    hTabularThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hTabularThermo_H
#define hTabularThermo_H

#include "scalar.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class hTabularThermo;

template<class EquationOfState>
inline hTabularThermo<EquationOfState> operator+
(
    const hTabularThermo<EquationOfState>&,
    const hTabularThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabularThermo<EquationOfState> operator-
(
    const hTabularThermo<EquationOfState>&,
    const hTabularThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabularThermo<EquationOfState> operator*
(
    const scalar,
    const hTabularThermo<EquationOfState>&
);

template<class EquationOfState>
inline hTabularThermo<EquationOfState> operator==
(
    const hTabularThermo<EquationOfState>&,
    const hTabularThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hTabularThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                      Class hTabularThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hTabularThermo
:
    public EquationOfState
{
    // Private data

        //- Heat of formation
        //  Note: input in [J/kg], but internally uses [J/kmol]
        scalar Hf_;
	
	interpolation2DTable<scalar> cpTable;
	//interpolation2DTable<scalar> hTable; //Kristjan: redundant, governed by TEqn


    // Private Member Functions

        //- Construct from components
        inline hTabularThermo
        (
            const EquationOfState& pt,
            const scalar Hf
        );


public:

    // Constructors

        //- Construct from Istream
        hTabularThermo(Istream& is);

        //- Construct from dictionary
        hTabularThermo(const dictionary& dict);

        //- Construct as copy
        inline hTabularThermo(const hTabularThermo&);

        //- Construct as a named copy
        inline hTabularThermo(const word&, const hTabularThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hTabular<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
            inline scalar hc() const;
	    
            inline scalar s(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline hTabularThermo& operator=(const hTabularThermo&);
        inline void operator+=(const hTabularThermo&);
        inline void operator-=(const hTabularThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend hTabularThermo operator+ <EquationOfState>
        (
            const hTabularThermo&,
            const hTabularThermo&
        );

        friend hTabularThermo operator- <EquationOfState>
        (
            const hTabularThermo&,
            const hTabularThermo&
        );

        friend hTabularThermo operator* <EquationOfState>
        (
            const scalar,
            const hTabularThermo&
        );

        friend hTabularThermo operator== <EquationOfState>
        (
            const hTabularThermo&,
            const hTabularThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hTabularThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hTabularThermoI.H"

#ifdef NoRepository
#   include "hTabularThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
