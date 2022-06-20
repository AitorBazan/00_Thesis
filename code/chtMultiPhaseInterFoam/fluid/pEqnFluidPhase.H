{
    bool closedVolume = p_rgh.needReference();
    rho = thermol.rho();

    rhok = rho - rhoR;

    if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }

    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));

    volVectorField HbyA("HbyA", U);

    HbyA = rAU()*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
        + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi)
    );
    

    if (p_rgh.needReference())
    {
        fvc::makeRelative(phiHbyA, U);
        adjustPhi(phiHbyA, U, p_rgh);
        fvc::makeAbsolute(phiHbyA, U);
    }

    surfaceScalarField phig
    (
        (
            thermol.surfaceTensionForce()
            - ghf*fvc::snGrad(rhok)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {

        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA)
            - fvm::laplacian(rAUf, p_rgh)
        );


        if (thermol.includeVolChange())
        {
            p_rghEqn += thermol.volTransfer(p_rgh);
        }

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve
        (
            mesh.solver
            (
                p_rgh.select
                (
                    (
                        oCorr == nOuterCorr-1
                        && corr == nCorr-1
                        && nonOrth == nNonOrthCorr
                    )
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {

            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU()*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);

            U.correctBoundaryConditions();

            fvOptions.correct(U);

            K = 0.5*magSqr(U);

        }

    }

    #include "incompressibleContinuityErrors.H"

    fvc::makeRelative(phi, U);
    p == p_rgh + rhok*gh;


    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }

    if (!correctPhi)
    {
        rAU.clear();
    }
    rho = thermol.rho();

// Update pressure time derivative if needed
    if (thermol.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}
