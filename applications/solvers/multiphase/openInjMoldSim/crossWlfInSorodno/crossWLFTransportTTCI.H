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
inline void Foam::crossWLFTransportTTC<Thermo>::calcCoeffs
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
inline Foam::crossWLFTransportTTC<Thermo>::crossWLFTransportTTC
(
    const Thermo& t,
    const scalar n,
    const scalar Tau,
    const scalar D1,
    const scalar D2,
    const scalar D3,
    const scalar D3s,
    const scalar A1,
    const scalar A2,
    const interpolation2DTable<scalar> kappa,
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
    D3s_(D3s),
    A1_(A1),
    A2_(A2),
    kappa_(kappa),
    etaMin_(etaMin),
    etaMax_(etaMax),
    TnoFlow_(TnoFlow),
    deltaTempInterp_(deltaTempInterp)
{
    kappa_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}

/*
template<class Thermo>
inline Foam::crossWLFTransportTTC<Thermo>::crossWLFTransportTTC
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
inline Foam::crossWLFTransportTTC<Thermo>::crossWLFTransportTTC
(
    const word& name,
    const crossWLFTransportTTC& st
)
:
    Thermo(name, st),
    n_(st.n_),
    Tau_(st.Tau_),
    D1_(st.D1_),
    D2_(st.D2_),
    D3_(st.D3_),
    D3s_(st.D3s_),
    A1_(st.A1_),
    A2_(st.A2_),
    kappa_(st.kappa),
    etaMin_(st.etaMin_),
    etaMax_(st.etaMax_),
    TnoFlow_(st.TnoFlow),
    deltaTempInterp_(st.deltaTempInterp)
{
    kappa_.outOfBounds(interpolation2DTable<scalar>::CLAMP);
}


template<class Thermo>
inline Foam::autoPtr<Foam::crossWLFTransportTTC<Thermo> >
Foam::crossWLFTransportTTC<Thermo>::clone() const
{
    return autoPtr<crossWLFTransportTTC<Thermo> >
    (
        new crossWLFTransportTTC<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::crossWLFTransportTTC<Thermo> >
Foam::crossWLFTransportTTC<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<crossWLFTransportTTC<Thermo> >
    (
        new crossWLFTransportTTC<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::crossWLFTransportTTC<Thermo> >
Foam::crossWLFTransportTTC<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<crossWLFTransportTTC<Thermo> >
    (
        new crossWLFTransportTTC<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::crossWLFTransportTTC<Thermo>::mu
(
    const scalar p,
    const scalar T,
    const scalar strig
) const
{
    scalar TnoFlowS = TnoFlow_ + D3s_ * p;
    if (T < TnoFlowS) { return etaMax_; } //leave immediately to avoid n0 = infinite
    const scalar n0 = D1_*::exp(-(A1_*(T-D2_-D3_*p))/(A2_+T-D2_));
    scalar visk = n0/(1.0+pow(n0*strig/Tau_,1.0-n_));
    scalar Tstart = TnoFlowS + deltaTempInterp_; //solidification onset temperature .. D3s*p from TnoFlowS

    //Info << "visk= " << visk  << " p=" << p << " T=" << T << " strig=" << strig << endl;   // debug Kristjan

    if (visk > etaMax_) { visk = etaMax_; }
    if (visk < etaMin_) { visk = etaMin_; }

    if (T < Tstart) //interpolated increase of viscosity
        {
            visk = pow(etaMax_, (Tstart-T)/deltaTempInterp_)
                  *  pow(visk, (T-TnoFlowS)/deltaTempInterp_);//weighted avg in log space
        }

    return visk;
}


template<class Thermo>
inline Foam::scalar Foam::crossWLFTransportTTC<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    //scalar Cv_ = this->Cv(p, T);
    //return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    return kappa_(T, p); // this->Cpv(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::crossWLFTransportTTC<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::crossWLFTransportTTC<Thermo>&
Foam::crossWLFTransportTTC<Thermo>::operator=
(
    const crossWLFTransportTTC<Thermo>& st
)
{
    Thermo::operator=(st);

    n_ = st.n_;
    Tau_ = st.Tau_;
    D1_ = st.D1_;
    D2_ = st.D2_;
    D3_ = st.D3_;
    D3s_ = st.D3s_;
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
inline void Foam::crossWLFTransportTTC<Thermo>::operator+=
(
    const crossWLFTransportTTC<Thermo>& st
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
inline void Foam::crossWLFTransportTTC<Thermo>::operator-=
(
    const crossWLFTransportTTC<Thermo>& st
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
inline void Foam::crossWLFTransportTTC<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::crossWLFTransportTTC<Thermo> Foam::operator+
(
    const crossWLFTransportTTC<Thermo>& st1,
    const crossWLFTransportTTC<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return crossWLFTransportTTC<Thermo>
    (
        t,
        molr1*st1.As_ + molr2*st2.As_,
        molr1*st1.Ts_ + molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::crossWLFTransportTTC<Thermo> Foam::operator-
(
    const crossWLFTransportTTC<Thermo>& st1,
    const crossWLFTransportTTC<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return crossWLFTransportTTC<Thermo>
    (
        t,
        molr1*st1.As_ - molr2*st2.As_,
        molr1*st1.Ts_ - molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::crossWLFTransportTTC<Thermo> Foam::operator*
(
    const scalar s,
    const crossWLFTransportTTC<Thermo>& st
)
{
    return crossWLFTransportTTC<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.As_,
        st.Ts_
    );
}


template<class Thermo>
inline Foam::crossWLFTransportTTC<Thermo> Foam::operator==
(
    const crossWLFTransportTTC<Thermo>& st1,
    const crossWLFTransportTTC<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
