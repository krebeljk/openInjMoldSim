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

volScalarField shrRate
(
    IOobject
    (
        "shrRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(symm(fvc::grad(U))&&symm(fvc::grad(U)))
);


#include "createPhi.H"

Info<< "Constructing mojTwoPhaseMixtureThermo\n" << endl;
mojTwoPhaseMixtureThermo mixture(mesh);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

volScalarField& p = mixture.p();
volScalarField& T = mixture.T();
volScalarField& strig = mixture.thermo1().strig();
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
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*mixture.thermo1().mu() + alpha2*mixture.thermo2().mu()
);

//kristjan - residuals p_rgh

volScalarField p_rgh_resid
(
    IOobject
    (
        "p_rgh_resid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("initOne", dimensionSet(0,0,-1,0,0,0,0), 1.0)
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
