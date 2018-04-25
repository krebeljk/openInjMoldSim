/*---------------------------------------------------------------------------*\
kristjan modified: OpenFOAM-3.0.0/src/thermophysicalModels/basic/rhoThermo/rhoThermos.C
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "mojRhoThermo.H"
#include "mojMakeThermo.H"

#include "specie.H"
#include "mojPerfectGas.H"
#include "incompressiblePerfectGas.H"
#include "rhoConst.H"
#include "PengRobinsonGas.H"
#include "adiabaticPerfectFluid.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "mojThermo.H"

#include "constTransport.H"
#include "mojConstTransport.H"
#include "sutherlandTransport.H"
#include "polymerPVT.H"
#include "taitHCR.H"
#include "crossWLFTransport.H"
#include "crossWLFTransportTTC.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"
#include "hTabularThermo.H"


#include "heRhoThermo.H"
#include "mojHeRhoThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */



mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    mojConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    mojPerfectGas,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    mojConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    polymerPVT,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    mojConstTransport,
    sensibleInternalEnergy,
    hConstThermo,
    taitHCR,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransport,
    sensibleInternalEnergy,
    hConstThermo,
    polymerPVT,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransport,
    sensibleInternalEnergy,
    hConstThermo,
    taitHCR,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransport,
    sensibleInternalEnergy,
    hPolynomialThermo,//new
    polymerPVT,
    specie
);


mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransport,
    sensibleInternalEnergy,
    hPolynomialThermo,//new
    taitHCR,
    specie
);
mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransport,
    sensibleInternalEnergy,
    hTabularThermo,//new
    polymerPVT,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransport,
    sensibleInternalEnergy,
    hTabularThermo,//new
    taitHCR,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransportTTC,//new
    sensibleInternalEnergy,
    hPolynomialThermo,
    polymerPVT,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransportTTC,//new
    sensibleInternalEnergy,
    hPolynomialThermo,
    taitHCR,
    specie
);
mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransportTTC,//new
    sensibleInternalEnergy,
    hTabularThermo,//new
    polymerPVT,
    specie
);

mojMakeThermo
(
    mojRhoThermo,
    mojHeRhoThermo,
    pureMixture,
    crossWLFTransportTTC,//new
    sensibleInternalEnergy,
    hTabularThermo,//new
    taitHCR,
    specie
);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
