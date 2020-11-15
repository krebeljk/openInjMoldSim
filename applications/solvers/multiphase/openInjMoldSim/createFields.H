
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Constructing mojTwoPhaseMixtureThermo\n" << endl;
mojTwoPhaseMixtureThermo mixture(mesh);
mixture.correct();

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

volScalarField& p = mixture.p();
volScalarField& T = mixture.T();
volScalarField& shrRate = mixture.shrRate();
volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& psi1 = mixture.thermo1().psi();
volScalarField& rho2 = mixture.thermo2().rho();
const volScalarField& psi2 = mixture.thermo2().psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);


dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());

//kristjan - za HTC

volScalarField mojKappaOut
(
    IOobject
    (
        "mojKappaOut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.kappa()
);

//kristjan - pregled viskoznosti

volScalarField visc
(
    IOobject
    (
        "visc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.mu()
);

//kristjan - Calc Celsius
dimensionedScalar twoSevenThree("twoSevenThree", dimensionSet(0,0,0,1,0,0,0), 273.15);
volScalarField Tc
(
    IOobject
    (
        "Tc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    T - twoSevenThree
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt
(
    IOobject
    (
        "dgdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
);

// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture);

// Construct compressible turbulence model
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New(rho, U, rhoPhi, mixture)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));



Info<< "Reading solidification properties\n" << endl;
    const word sldDictName("solidificationProperties");

    // Create and input-output object - this holds the path to the dict and its name
    IOobject sldDictIO
    (
        sldDictName, // name of the file
        mesh.time().constant(), // path to where the file is
        mesh, // reference to the mesh needed by the constructor
        IOobject::READ_IF_PRESENT // indicate that reading this dictionary is optional
    );

    // Check the if the dictionary is present and follows the OF format
    if (!sldDictIO.typeHeaderOk<IOdictionary>()) //VER MELHOR
    {
        Info << "No constant/solidificationProperties dictionary found." << endl;
        Info << "Flow cessation modeled as visc(TnoFlow)=etaMax." << nl << nl << endl;
    }

    // Get access to the dictionary
    dictionary sldDict;
    // Initialise the dictionary object
    sldDict = IOdictionary(sldDictIO);

    // This template method needs to know the type of the variable and can provide
    // a default value if the entry is not found in the dictionary
    dimensionedScalar shrMod("shearModulus", dimensionSet(1,-1,-2,0,0,0,0),
        sldDict.lookupOrDefault<scalar>("shearModulus", 900e6)
    );
    dimensionedScalar shrRateLimEl("shrRateLimEl", dimensionSet(0,0,-1,0,0,0,0),
        sldDict.lookupOrDefault<scalar>("shrRateLimEl", 1.0)
    );
    dimensionedScalar viscLimEl("viscLimEl", dimensionSet(1,-1,-1,0,0,0,0),
        sldDict.lookupOrDefault<scalar>("viscLimEl", 1e6)
    );

    if (sldDictIO.typeHeaderOk<IOdictionary>()) //VER MELHOR
    {
        // Summarise what's been read and print in the console
        Info << nl << "Solidification properties: " << nl
             << "Shear modulus: " << shrMod << nl
             << "Highest shear rate with elasticity: " << shrRateLimEl << nl
             << "Lowest viscosity with elasticity: " << viscLimEl << nl
             << nl << endl;
    }


    volSymmTensorField elSigDev
    (
        IOobject
        (
            "elSigDev",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE
        )
        , mesh
        , dimensionedSymmTensor("zero", dimensionSet(1,-1,-2,0,0,0,0), symmTensor::zero)
        , "zeroGradient"
    );

// pressure preventing material collapse bellow pMin
volScalarField pAux
(
    IOobject
    (
        "pAux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeropAux", dimPressure, 0.0)
);
