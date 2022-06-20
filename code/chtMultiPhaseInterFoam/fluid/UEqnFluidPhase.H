
rhok = rho - rhoR;

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + turbulence.divDevRhoReff(U)
  ==
    fvOptions(rho, U)
);

UEqn.relax();

thermol.addInterfacePorosity(UEqn);

if (momentumPredictor)
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
                thermol.surfaceTensionForce()
                - ghf*fvc::snGrad(rhok)
                - fvc::snGrad(p_rgh) 
            ) * mesh.magSf()
        ),
        mesh.solver(U.select(finalIter))
    );

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
