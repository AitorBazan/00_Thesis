// tmp<volScalarField> rAU;

PtrList<volScalarField> rAUFluid(fluidRegions.size());
forAll(rAUFluid, i)
{
    const fvMesh& mesh = fluidRegions[i];

    if (correctPhiFluid[i])
    {
        rAUFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rAU",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("rAU", dimTime/dimDensity, 1)
            )
        );

        #include "incompCorrectPhi.H"
    }
    else
    {
        CorrectPhi
        (
            UFluid[i],
            phiFluid[i],
            p_rghFluid[i],
            dimensionedScalar("rAUf", dimTime/rhoFluid[i].dimensions(), 1),
            geometricZeroField(),
            pimpleFluid[i]
        );

        #include "incompressibleContinuityErrors.H"
    }

}

UPtrList<volScalarField> rAU(rAUFluid.size());
forAll(rAUFluid, i)
{
    rAU.set(i, rAUFluid[i]);
}