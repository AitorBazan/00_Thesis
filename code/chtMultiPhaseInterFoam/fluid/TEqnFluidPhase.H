{
    rhoCp =  rho*thermol.Cp();

    rhok = rho - rhoR;
    const surfaceScalarField rhoCpPhi(fvc::interpolate(thermol.Cp())*rhoPhi);
    volTensorField gradU = fvc::grad(U);
    volTensorField tau = turbulence.muEff() * (gradU + gradU.T());
    const volScalarField kappaEff
    (
        "kappaEff",
        thermol.kappa() + thermol.Cp()*turbulence.mut()/thermol.Prt()
    );

    
    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T, "div(phi,T)")
      + fvc::div(rhoPhi/fvc::interpolate(rho),p,"div(phiv,p)")
      + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)      
      - fvm::laplacian(kappaEff, T,  "laplacian(kappa,T)")
      ==
        thermol.heatTransfer(T)
        + fvc::div(tau & U, "div(tau,U)")
        + rho*(U&g)
        + fvOptions(rhoCp, T)
    );
    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve(mesh.solver(T.select(finalIter)));

    fvOptions.correct(T);
    thermol.correct();

    Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;

}

