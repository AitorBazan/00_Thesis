if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rhosol, h)
          - (
               therms.isotropic()
             ? fvm::laplacian(betav*therms.alpha(), h, "laplacian(alpha,h)")
             : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rhosol, h)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve(mesh.solver(h.select(finalIter)));

        fvOptions.correct(h);
    }

    therms.correct();

    Info<< "Min/max T:" << min(therms.T()).value() << ' '
        << max(therms.T()).value() << endl;
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
