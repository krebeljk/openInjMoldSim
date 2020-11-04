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
    Foam::crossWLFTransport

Description
    Transport package using crossWLF's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    crossWLFTransportI.H
    crossWLFTransport.C

\*---------------------------------------------------------------------------*/

#ifndef crossWLFTransport_H
#define crossWLFTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class crossWLFTransport;

template<class Thermo>
inline crossWLFTransport<Thermo> operator+
(
    const crossWLFTransport<Thermo>&,
    const crossWLFTransport<Thermo>&
);

template<class Thermo>
inline crossWLFTransport<Thermo> operator-
(
    const crossWLFTransport<Thermo>&,
    const crossWLFTransport<Thermo>&
);

template<class Thermo>
inline crossWLFTransport<Thermo> operator*
(
    const scalar,
    const crossWLFTransport<Thermo>&
);

template<class Thermo>
inline crossWLFTransport<Thermo> operator==
(
    const crossWLFTransport<Thermo>&,
    const crossWLFTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const crossWLFTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class crossWLFTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class crossWLFTransport
:
    public Thermo
{
    // Private data

        // crossWLF's coefficients
        scalar n_, Tau_, D1_, D2_, D3_, A1_, A2_, kappa_;
        scalar etaMin_, etaMax_, TnoFlow_, deltaTempInterp_;


    // Private Member Functions

        //- Calculate the crossWLF coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );


public:

    // Constructors

        //- Construct from components
        inline crossWLFTransport
        (
            const Thermo& t,
            const scalar n_,
            const scalar Tau_,
            const scalar D1_,
            const scalar D2_,
            const scalar D3_,
            const scalar A1_,
            const scalar A2_,
            const scalar kappa_,
            const scalar etaMin_,
            const scalar etaMax_,
            const scalar TnoFlow_,
            const scalar deltaTempInterp_
        );

        //- Construct from two viscosities
        inline crossWLFTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Construct as named copy
        inline crossWLFTransport(const word&, const crossWLFTransport&);

        //- Construct from Istream
        crossWLFTransport(Istream&);

        //- Construct from dictionary
        crossWLFTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<crossWLFTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<crossWLFTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<crossWLFTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "crossWLF<" + Thermo::typeName() + '>';
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

        inline crossWLFTransport& operator=(const crossWLFTransport&);

        inline void operator+=(const crossWLFTransport&);

        inline void operator-=(const crossWLFTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend crossWLFTransport operator+ <Thermo>
        (
            const crossWLFTransport&,
            const crossWLFTransport&
        );

        friend crossWLFTransport operator- <Thermo>
        (
            const crossWLFTransport&,
            const crossWLFTransport&
        );

        friend crossWLFTransport operator* <Thermo>
        (
            const scalar,
            const crossWLFTransport&
        );

        friend crossWLFTransport operator== <Thermo>
        (
            const crossWLFTransport&,
            const crossWLFTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const crossWLFTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "crossWLFTransportI.H"

#ifdef NoRepository
#   include "crossWLFTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
