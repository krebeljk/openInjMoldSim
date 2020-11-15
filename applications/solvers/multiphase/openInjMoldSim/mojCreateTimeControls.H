#include "createTimeControls.H"

scalar maxSolidCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxSolidCo", 0.005); // Kristjan
scalar pAuxRlx =
    runTime.controlDict().lookupOrDefault<scalar>("pAuxRlx", 0.1); // Kristjan
