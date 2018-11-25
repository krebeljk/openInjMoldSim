/*---------------------------------------------------------------------------*\
modified from: OpenFOAM-3.0.0/src/thermophysicalModels/specie/transport/sutherland/sutherlandTransport.H
 =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::taitTg

Description
    Perfect gas equation of state.

SourceFiles
    taitTgI.H
    taitTg.C

\*---------------------------------------------------------------------------*/

#ifndef taitTg_H
#define taitTg_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class taitTg;

template<class Specie>
inline taitTg<Specie> operator+
(
    const taitTg<Specie>&,
    const taitTg<Specie>&
);

template<class Specie>
inline taitTg<Specie> operator-
(
    const taitTg<Specie>&,
    const taitTg<Specie>&
);

template<class Specie>
inline taitTg<Specie> operator*
(
    const scalar,
    const taitTg<Specie>&
);

template<class Specie>
inline taitTg<Specie> operator==
(
    const taitTg<Specie>&,
    const taitTg<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const taitTg<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class taitTg Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class taitTg
:
    public Specie
{
    // Private data

        // WLF coefficients
        scalar  D1_, D2_, D3_, A1_, A2_;

        //- Fluid constant
        scalar b1m_;
        scalar b2m_;
        scalar b3m_;
        scalar b4m_;
        scalar b1s_;
        scalar b2s_;
        scalar b3s_;
        scalar b4s_;
        scalar b5_;
        scalar b6_;
        scalar b10_;
        static scalar C_;// definicija v .C // univerzalna konstanta v Taitovi enacbi

public:

    // Constructors

        //- Construct from components
        inline taitTg
        (
            const Specie& sp,
            const scalar D1_,
            const scalar D2_,
            const scalar D3_,
            const scalar A1_,
            const scalar A2_,
            const scalar b1m,
            const scalar b2m,
            const scalar b3m,
            const scalar b4m,
            const scalar b1s,
            const scalar b2s,
            const scalar b3s,
            const scalar b4s,
            const scalar b5,
            const scalar b6,
            const scalar b10
        );

        //- Construct from Istream
        taitTg(Istream&);

        //- Construct from dictionary
        taitTg(const dictionary& dict);

        //- Construct as named copy
        inline taitTg(const word& name, const taitTg&);

        //- Construct and return a clone
        inline autoPtr<taitTg> clone() const;

        // Selector from Istream
        inline static autoPtr<taitTg> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<taitTg> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "taitTg<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return fluid constant [J/(kg K)]
            inline scalar R() const;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return glassy specific volume [kg/m^3]
            inline scalar vg(scalar p, scalar T) const;

            //- Return equilibrium specific volume [kg/m^3]
            inline scalar veq(scalar p, scalar T) const;

            //- Return equilibrium free specific volume [kg/m^3]
            inline scalar vfeq(scalar p, scalar T) const;

            //- Return entropy [J/(kmol K)]
            inline scalar s(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;

            //- Return coef. of vol. therm. expansion [J/(kmol K]
            inline scalar cTExp(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const taitTg&);
        inline void operator-=(const taitTg&);

        inline void operator*=(const scalar);


    // Friend operators

        friend taitTg operator+ <Specie>
        (
            const taitTg&,
            const taitTg&
        );

        friend taitTg operator- <Specie>
        (
            const taitTg&,
            const taitTg&
        );

        friend taitTg operator* <Specie>
        (
            const scalar s,
            const taitTg&
        );

        friend taitTg operator== <Specie>
        (
            const taitTg&,
            const taitTg&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const taitTg&
        );

    // Izracun Entropije
private:

            //- Solidification interpolation to continuity
            inline scalar interpolation
            (
                 scalar vS
                ,scalar vM
                ,scalar p
                ,scalar T
            ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "taitTgI.H"

#ifdef NoRepository
#   include "taitTg.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
