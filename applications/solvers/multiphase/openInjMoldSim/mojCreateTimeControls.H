#include "createTimeControls.H"

scalar maxSolidCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxSolidCo", 1.0); // Kristjan
