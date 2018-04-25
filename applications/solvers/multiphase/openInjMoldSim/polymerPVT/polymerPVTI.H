/*---------------------------------------------------------------------------*\
modified from: OpenFOAM-3.0.0/src/thermophysicalModels/specie/transport/sutherland/sutherlandTransportI.H
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

\*---------------------------------------------------------------------------*/

#include "polymerPVT.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::polymerPVT<Specie>::polymerPVT
(
    const Specie& sp,
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
    const scalar b7,
    const scalar b8,
    const scalar b9
)
:
    Specie(sp),
    b1m_(b1m),
    b2m_(b2m),
    b3m_(b3m),
    b4m_(b4m),
    b1s_(b1s),
    b2s_(b2s),
    b3s_(b3s),
    b4s_(b4s),
    b5_(b5),
    b6_(b6),
    b7_(b7),
    b8_(b8),
    b9_(b9)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::polymerPVT<Specie>::polymerPVT
(
    const word& name,
    const polymerPVT<Specie>& pf
)
:
    Specie(name, pf),
    b1m_(pf.b1m_),
    b2m_(pf.b2m_),
    b3m_(pf.b3m_),
    b4m_(pf.b4m_),
    b1s_(pf.b1s_),
    b2s_(pf.b2s_),
    b3s_(pf.b3s_),
    b4s_(pf.b4s_),
    b5_(pf.b5_),
    b6_(pf.b6_),
    b7_(pf.b7_),
    b8_(pf.b8_),
    b9_(pf.b9_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::polymerPVT<Specie> >
Foam::polymerPVT<Specie>::clone() const
{
    return autoPtr<polymerPVT<Specie> >(new polymerPVT<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::polymerPVT<Specie> >
Foam::polymerPVT<Specie>::New(Istream& is)
{
    return autoPtr<polymerPVT<Specie> >(new polymerPVT<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::polymerPVT<Specie> >
Foam::polymerPVT<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<polymerPVT<Specie> >(new polymerPVT<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::R() const
{
        Info << "POZOR: klic R() v polymerPVTI.H" << endl; // ce se to pjavi v izpisu simulacije bo treba preveriti, kaj so posledice te funkcije
    //return R_;
    return 1;
}


template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::rho(scalar p, scalar T) const // mathematica export fim_01_01_E_04_02
{
    scalar spec_volS;
    scalar spec_volM;
    scalar spec_vol;

    //melt
     spec_volM = (b1m_+b2m_*(-b5_+T))*(1-C_*log(1+p*pow(b3m_,-1)*exp(b4m_*(-b5_+T))));
    //solid
     spec_volS = (b1s_+b2s_*(-b5_+T))*(1-C_*log(1+p*pow(b3s_,-1)*exp(b4s_*(-b5_+T))))+b7_*exp(-(b9_*p)+b8_*(-b5_+T));

     spec_vol = interpolation(
          spec_volS
         ,spec_volM
         ,p
         ,T
     );
    return 1.0/spec_vol;
}

template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::cTExp(scalar p, scalar T) const // num. der.
{
    scalar spec_volS;
    scalar spec_volSp;
    scalar spec_volSm;
    scalar spec_volM;
    scalar spec_volMp;
    scalar spec_volMm;
    scalar coefThermExp;

    scalar dT = 2.0; //[K]

    //melt
     spec_volMp = (b1m_+b2m_*(-b5_+T+dT))*(1-C_*log(1+p*pow(b3m_,-1)*exp(b4m_*(-b5_+T+dT))));
     spec_volMm = (b1m_+b2m_*(-b5_+T-dT))*(1-C_*log(1+p*pow(b3m_,-1)*exp(b4m_*(-b5_+T-dT))));
     spec_volM  = (b1m_+b2m_*(-b5_+T   ))*(1-C_*log(1+p*pow(b3m_,-1)*exp(b4m_*(-b5_+T   ))));
    //solid
     spec_volSp = (b1s_+b2s_*(-b5_+T+dT))*(1-C_*log(1+p*pow(b3s_,-1)*exp(b4s_*(-b5_+T+dT))));
     spec_volSm = (b1s_+b2s_*(-b5_+T-dT))*(1-C_*log(1+p*pow(b3s_,-1)*exp(b4s_*(-b5_+T-dT))));
     spec_volS  = (b1s_+b2s_*(-b5_+T   ))*(1-C_*log(1+p*pow(b3s_,-1)*exp(b4s_*(-b5_+T   ))));

     coefThermExp = interpolation( // this eliminates crystallization effect on thermal volume change
          (spec_volSp-spec_volSm)/(spec_volS*2.0*dT)
         ,(spec_volMp-spec_volMm)/(spec_volM*2.0*dT)
         ,p
         ,T
     );
    return coefThermExp;
}


template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::s(scalar p, scalar T) const // POZOR - NI TESTIRANO
{
    scalar pZ = (T - b5_)/b6_; //Tlak na meji med S in M domenama
    scalar entr;

    if (T < b5_ && p > 0.0) //solid
    {
        entr = sIntS(0.0, p, T);
    }
    else if (T < b5_ && p < 0.0 && p > pZ)
    {
        entr = sIntS(0.0, p, T);
    }
    else if (T < b5_ && p < pZ)
    {
        entr = sIntS(0.0, pZ, T) + sIntM(pZ, p, T);
    }
    else if (T > b5_ && p < 0.0)  //melt
    {
        entr = sIntM(0.0, p, T);
    }
    else if (T > b5_ && p > 0.0 && p < pZ)
    {
        entr = sIntM(0.0, p, T);
    }
    else if (T > b5_ && p > pZ)
    {
        entr = sIntM(0.0, pZ, T) + sIntS(pZ, p, T);
    }

    Info << "Ni testirano: s(J/kgK)=" << entr << "s(J/kmolK)=" << entr*this->W() << " p=" << p << " T=" << T << endl;

    return entr*this->W(); // pozor, je v kmol
}

// A (b1m_+b2m_*(T-b5_))
// B (b3m_*exp(-b4m_*(T-b5_)))
 
template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::psi(scalar p, scalar T) const // mathematica export fim_01_01_E_04_02
{
    if (p < scalar(10.0)) //majhen tlak
    {
        Info << "Pozor! Majhen tlak, nevarnost deljenja z 0 v enacbi stanja, funkcija psi. p=" << p << " Pa" << endl;
    }

    scalar psiS;
    scalar psiM;
    scalar psiRet;

    psiS = -((-(b7_*b9_*exp(-(b5_*b8_)-b9_*p+b8_*T))-C_*(b1s_+b2s_*(-b5_+T))*exp(b4s_*T)*pow(b3s_*exp(b4s_*b5_)+p*exp(b4s_*T),-1))*pow((b1s_+b2s_*(-b5_+T))*(1-C_*log(1+p*pow(b3s_,-1)*exp(b4s_*(-b5_+T))))+b7_*exp(-(b5_*b8_)-b9_*p+b8_*T),-2));
    psiM = C_*exp(b4m_*(-b5_+T))*pow(b1m_+b2m_*(-b5_+T),-1)*pow(-1+C_*log(1+p*pow(b3m_,-1)*exp(b4m_*(-b5_+T))),-2)*pow(b3m_+p*exp(b4m_*(-b5_+T)),-1);
    psiRet = interpolation(
          psiS
         ,psiM
         ,p
         ,T
     );
    return psiRet;
}


template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::Z(scalar p, scalar T) const
{
    Info << "POZOR: Z(p, T) v polymerPVTI.H" << endl; // ce se to pjavi v izpisu simulacije bo treba preveriti, kaj so posledice te funkcije
    return 1;
}


template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::cpMcv(scalar p, scalar T) const //fim_01_01_E_04_01_cpMcv.txt - MMA generated
{
    scalar cpmcv;
    scalar cpmcvMol;
    scalar cpmcvM;
    scalar cpmcvS;
    //melt
    cpmcvM = T*pow(C_,-1)*pow(b1m_+b2m_*(-b5_+T),-1)*pow(p+b3m_*exp(b4m_*(b5_-T)),-1)*pow(b1m_*b4m_*C_*p+b2m_*C_*log(1+p*pow(b3m_,-1)*exp(b4m_*(-b5_+T)))*(p+b3m_*exp(b4m_*(b5_-T)))-b2m_*(p+b4m_*b5_*C_*p-b4m_*C_*p*T+b3m_*exp(b4m_*(b5_-T))),2);

    //solid
    cpmcvS = T*exp(-(b5_*b8_)-b9_*p-2*b4s_*T)*pow(p+b3s_*exp(b4s_*(b5_-T)),-1)*pow(C_*(b1s_+b2s_*(-b5_+T))*exp(b5_*b8_+b9_*p)+b7_*b9_*(p+b3s_*exp(b4s_*(b5_-T)))*exp(b8_*T),-1)*pow(b2s_*exp(b5_*b8_+b9_*p)*(b3s_*exp(b4s_*b5_)+p*(1+b4s_*C_*(b5_-T))*exp(b4s_*T))+b7_*b8_*p*exp((b4s_+b8_)*T)-b2s_*C_*log(1+p*pow(b3s_,-1)*exp(b4s_*(-b5_+T)))*exp(b9_*p)*(b3s_*exp(b5_*(b4s_+b8_))+p*exp(b5_*b8_+b4s_*T))-b1s_*b4s_*C_*p*exp(b5_*b8_+b9_*p+b4s_*T)+b3s_*b7_*b8_*exp(b4s_*b5_+b8_*T),2);

     cpmcv = interpolation(
          cpmcvS
         ,cpmcvM
         ,p
         ,T
     );
    cpmcvMol = cpmcv * this->W();
    //Info << "cpmcvM[J/(kmol K)]=" << cpmcvMol << " p=" << p << " T=" << T << endl;
    return cpmcvMol;
}


template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::sIntS(scalar p1, scalar p2, scalar T) const
{
    return //fim_01_01_E_04_02_entropija.txt - koda generirana iz mathematice
    b2s_*p1-b2s_*p2-b4s_*C_*p1*(b1s_+b2s_*(-b5_+T))+b4s_*C_*p2*(b1s_+b2s_*(-b5_+T))+b3s_*b4s_*C_*(b1s_+b2s_*(-b5_+T))*log(b3s_*exp(b4s_*b5_)+p1*exp(b4s_*T))*exp(b4s_*(b5_-T))-b3s_*b4s_*C_*(b1s_+b2s_*(-b5_+T))*log(b3s_*exp(b4s_*b5_)+p2*exp(b4s_*T))*exp(b4s_*(b5_-T))+b2s_*C_*(p1+b3s_*exp(b4s_*(b5_-T)))-b2s_*C_*log(1+p1*pow(b3s_,-1)*exp(b4s_*(-b5_+T)))*(p1+b3s_*exp(b4s_*(b5_-T)))-b2s_*C_*(p2+b3s_*exp(b4s_*(b5_-T)))+b2s_*C_*log(1+p2*pow(b3s_,-1)*exp(b4s_*(-b5_+T)))*(p2+b3s_*exp(b4s_*(b5_-T)))-b7_*b8_*pow(b9_,-1)*exp(-(b5_*b8_)-b9_*p1+b8_*T)+b7_*b8_*pow(b9_,-1)*exp(-(b5_*b8_)-b9_*p2+b8_*T);
}

template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::sIntM(scalar p1, scalar p2, scalar T) const
{
    return//fim_01_01_E_04_02_entropija.txt - koda generirana iz mathematice
    b2m_*p1-b2m_*p2-b4m_*C_*p1*(b1m_+b2m_*(-b5_+T))+b4m_*C_*p2*(b1m_+b2m_*(-b5_+T))+b3m_*b4m_*C_*(b1m_+b2m_*(-b5_+T))*log(b3m_*exp(b4m_*b5_)+p1*exp(b4m_*T))*exp(b4m_*(b5_-T))-b3m_*b4m_*C_*(b1m_+b2m_*(-b5_+T))*log(b3m_*exp(b4m_*b5_)+p2*exp(b4m_*T))*exp(b4m_*(b5_-T))+b2m_*C_*(p1+b3m_*exp(b4m_*(b5_-T)))-b2m_*C_*log(1+p1*pow(b3m_,-1)*exp(b4m_*(-b5_+T)))*(p1+b3m_*exp(b4m_*(b5_-T)))-b2m_*C_*(p2+b3m_*exp(b4m_*(b5_-T)))+b2m_*C_*log(1+p2*pow(b3m_,-1)*exp(b4m_*(-b5_+T)))*(p2+b3m_*exp(b4m_*(b5_-T)));

}

template<class Specie>
inline Foam::scalar Foam::polymerPVT<Specie>::interpolation(scalar vS ,scalar vM ,scalar p ,scalar T) const
{
    scalar Td = 1.0; //[K] half the interpolation band
    scalar T0 = b5_ + b6_*p; //Tait solidification temperature
    scalar Tp = T0 + Td; //upper band limit
    scalar Tm = T0 - Td; //lower band limit

    if ( T > Tp )//melt
    {
        return vM;
    }
    else if ( T > Tm)
    {
        return (T-Tm)/(Tp-Tm) * vM + (Tp-T)/(Tp-Tm) * vS;
    }
    else
    {
        return vS;
    }

}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
/*
template<class Specie>
inline void Foam::polymerPVT<Specie>::operator+=
(
    const polymerPVT<Specie>& pf
)
{
    scalar molr1 = this->nMoles();

    Specie::operator+=(pf);

    molr1 /= this->nMoles();
    scalar molr2 = pf.nMoles()/this->nMoles();

    R_ = 1.0/(molr1/R_ + molr2/pf.R_);
    rho0_ = molr1*rho0_ + molr2*pf.rho0_;
}


template<class Specie>
inline void Foam::polymerPVT<Specie>::operator-=
(
    const polymerPVT<Specie>& pf
)
{
    scalar molr1 = this->nMoles();

    Specie::operator-=(pf);

    molr1 /= this->nMoles();
    scalar molr2 = pf.nMoles()/this->nMoles();

    R_ = 1.0/(molr1/R_ - molr2/pf.R_);
    rho0_ = molr1*rho0_ - molr2*pf.rho0_;
}


template<class Specie>
inline void Foam::polymerPVT<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::polymerPVT<Specie> Foam::operator+
(
    const polymerPVT<Specie>& pf1,
    const polymerPVT<Specie>& pf2
)
{
    scalar nMoles = pf1.nMoles() + pf2.nMoles();
    scalar molr1 = pf1.nMoles()/nMoles;
    scalar molr2 = pf2.nMoles()/nMoles;

    return polymerPVT<Specie>
    (
        static_cast<const Specie&>(pf1)
      + static_cast<const Specie&>(pf2),
        1.0/(molr1/pf1.R_ + molr2/pf2.R_),
        molr1*pf1.rho0_ + molr2*pf2.rho0_
    );
}


template<class Specie>
inline Foam::polymerPVT<Specie> Foam::operator-
(
    const polymerPVT<Specie>& pf1,
    const polymerPVT<Specie>& pf2
)
{
    scalar nMoles = pf1.nMoles() + pf2.nMoles();
    scalar molr1 = pf1.nMoles()/nMoles;
    scalar molr2 = pf2.nMoles()/nMoles;

    return polymerPVT<Specie>
    (
        static_cast<const Specie&>(pf1)
      - static_cast<const Specie&>(pf2),
        1.0/(molr1/pf1.R_ - molr2/pf2.R_),
        molr1*pf1.rho0_ - molr2*pf2.rho0_
    );
}


template<class Specie>
inline Foam::polymerPVT<Specie> Foam::operator*
(
    const scalar s,
    const polymerPVT<Specie>& pf
)
{
    return polymerPVT<Specie>
    (
        s*static_cast<const Specie&>(pf),
        pf.R_,
        pf.rho0_
    );
}


template<class Specie>
inline Foam::polymerPVT<Specie> Foam::operator==
(
    const polymerPVT<Specie>& pf1,
    const polymerPVT<Specie>& pf2
)
{
    return pf2 - pf1;
}
*/

// ************************************************************************* //
