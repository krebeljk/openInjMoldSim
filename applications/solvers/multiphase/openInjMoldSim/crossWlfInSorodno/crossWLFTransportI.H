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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::crossWLFTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    //Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    //As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::crossWLFTransport<Thermo>::crossWLFTransport
(
    const Thermo& t,
    const scalar n,
    const scalar Tau,
    const scalar D1,
    const scalar D2,
    const scalar D3,
    const scalar A1,
    const scalar A2,
    const scalar kappa,
    const scalar etaMin,
    const scalar etaMax,
    const scalar TnoFlow,
    const scalar deltaTempInterp
)
:
    Thermo(t),
    n_(n),
    Tau_(Tau),
    D1_(D1),
    D2_(D2),
    D3_(D3),
    A1_(A1),
    A2_(A2),
    kappa_(kappa),
    etaMin_(etaMin),
    etaMax_(etaMax),
    TnoFlow_(TnoFlow),
    deltaTempInterp_(deltaTempInterp)
{}

/*
template<class Thermo>
inline Foam::crossWLFTransport<Thermo>::crossWLFTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}
*/

template<class Thermo>
inline Foam::crossWLFTransport<Thermo>::crossWLFTransport
(
    const word& name,
    const crossWLFTransport& st
)
:
    Thermo(name, st),
    n_(st.n_),
    Tau_(st.Tau_),
    D1_(st.D1_),
    D2_(st.D2_),
    D3_(st.D3_),
    A1_(st.A1_),
    A2_(st.A2_),
    kappa_(st.kappa_),
    etaMin_(st.etaMin_),
    etaMax_(st.etaMax_),
    TnoFlow_(st.TnoFlow),
    deltaTempInterp_(st.deltaTempInterp)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::crossWLFTransport<Thermo> >
Foam::crossWLFTransport<Thermo>::clone() const
{
    return autoPtr<crossWLFTransport<Thermo> >
    (
        new crossWLFTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::crossWLFTransport<Thermo> >
Foam::crossWLFTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<crossWLFTransport<Thermo> >
    (
        new crossWLFTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::crossWLFTransport<Thermo> >
Foam::crossWLFTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<crossWLFTransport<Thermo> >
    (
        new crossWLFTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::crossWLFTransport<Thermo>::mu
(
    const scalar p,
    const scalar T,
    const scalar shrRate
) const
{   
    if (T < TnoFlow_) { return etaMax_; } //leave immediately to avoid n0 = infinite
    const scalar n0 = D1_*::exp(-(A1_*(T-D2_-D3_*p))/(A2_+T-D2_));
    scalar visk = n0/(1.0+pow(n0*shrRate/Tau_,1.0-n_));
    scalar Tstart = TnoFlow_ + deltaTempInterp_; //solidification onset temperature

    //Info << "visk= " << visk  << " p=" << p << " T=" << T << " shrRate=" << shrRate << endl;   // debug Kristjan

    if (visk > etaMax_) { visk = etaMax_; }
    if (visk < etaMin_) { visk = etaMin_; }

    if (T < Tstart) //interpolated increase of viscosity
        {
            visk = etaMax_ * (Tstart-T)/deltaTempInterp_
                  +  visk *(T-TnoFlow_)/deltaTempInterp_;
        }

    return visk;
}


template<class Thermo>
inline Foam::scalar Foam::crossWLFTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    //scalar Cv_ = this->Cv(p, T);
    //return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    return kappa_;
}


template<class Thermo>
inline Foam::scalar Foam::crossWLFTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::crossWLFTransport<Thermo>&
Foam::crossWLFTransport<Thermo>::operator=
(
    const crossWLFTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    n_ = st.n_;
    Tau_ = st.Tau_;
    D1_ = st.D1_;
    D2_ = st.D2_;
    D3_ = st.D3_;
    A1_ = st.A1_;
    A2_ = st.A2_;
    kappa_ = st.kappa_;
    etaMin_ = st.etaMin_;
    etaMax_ = st.etaMax_;
    TnoFlow_ = st.TnoFlow_;
    deltaTempInterp_ = st.deltaTempInterp_;

    return *this;
}


template<class Thermo>
inline void Foam::crossWLFTransport<Thermo>::operator+=
(
    const crossWLFTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    //As_ = molr1*As_ + molr2*st.As_;
    //Ts_ = molr1*Ts_ + molr2*st.Ts_;
}


template<class Thermo>
inline void Foam::crossWLFTransport<Thermo>::operator-=
(
    const crossWLFTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    //As_ = molr1*As_ - molr2*st.As_;
    //Ts_ = molr1*Ts_ - molr2*st.Ts_;
}


template<class Thermo>
inline void Foam::crossWLFTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::crossWLFTransport<Thermo> Foam::operator+
(
    const crossWLFTransport<Thermo>& st1,
    const crossWLFTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return crossWLFTransport<Thermo>
    (
        t,
        molr1*st1.As_ + molr2*st2.As_,
        molr1*st1.Ts_ + molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::crossWLFTransport<Thermo> Foam::operator-
(
    const crossWLFTransport<Thermo>& st1,
    const crossWLFTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return crossWLFTransport<Thermo>
    (
        t,
        molr1*st1.As_ - molr2*st2.As_,
        molr1*st1.Ts_ - molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::crossWLFTransport<Thermo> Foam::operator*
(
    const scalar s,
    const crossWLFTransport<Thermo>& st
)
{
    return crossWLFTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.As_,
        st.Ts_
    );
}


template<class Thermo>
inline Foam::crossWLFTransport<Thermo> Foam::operator==
(
    const crossWLFTransport<Thermo>& st1,
    const crossWLFTransport<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
