CorrectPhi
(
    UFluid[i],
    phiFluid[i],
    p_rghFluid[i],
    surfaceScalarField("rAUf", fvc::interpolate(rAUFluid[i])),
    geometricZeroField(),
    pimpleFluid[i]
);

#include "incompressibleContinuityErrors.H"
