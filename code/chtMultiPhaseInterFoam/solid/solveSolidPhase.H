if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

// rhoCps =  rho*cp;


{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        // fvScalarMatrix TEqn
        // (
        //     fvm::ddt(rho*cp, T) - fvm::laplacian(kappa, T)
        //     ==
        //     fvOptions(rho*cp, T)
        // );

        fvScalarMatrix TEqn
        (
            fvm::ddt(betav*rhosol*cpsol, Tsol)
          - (
            //    thermo.isotropic()
               fvm::laplacian(betav*therms.alpha()*cpsol, Tsol, "laplacian(alpha,T)")
            //  : fvm::laplacian(betav*taniAlpha()*cp, T, "laplacian(alpha,T)")
            )
          ==
            fvOptions(rhosol*cpsol, Tsol)
        );

        TEqn.relax();

        fvOptions.constrain(TEqn);

        // TEqn.solve(mesh.solver(Tsol.select(finalIter)));
        TEqn.solve();

        fvOptions.correct(Tsol);

        therms.correct();

    }


    Info<< "Min/max T:" << min(Tsol).value() << ' '
        << max(Tsol).value() << endl;
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
