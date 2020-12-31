#include "createTimeControls.H"

scalar maxSolidCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxSolidCo", 1.0); // Kristjan
scalar pAuxRlx =
    runTime.controlDict().lookupOrDefault<scalar>("pAuxRlx", 0.1); // Kristjan
