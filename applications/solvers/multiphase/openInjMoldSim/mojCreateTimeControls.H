#include "createTimeControls.H"

scalar maxSolidCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxSolidCo", 0.005); // Kristjan
