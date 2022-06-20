scalar DiNum = -GREAT;

forAll(solidRegions, i)
{
    //- Note: do not use setRegionSolidFields.H to avoid double registering Cp
    //#include "setRegionSolidFields.H"
    const solidThermo& therms = thermos[i];

    tmp<volScalarField> magKappa;
    if (therms.isotropic())
    {
        magKappa = therms.kappa();
    }
    else
    {
        magKappa = mag(therms.Kappa());
    }

    tmp<volScalarField> tcp = therms.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> trho = therms.rho();
    const volScalarField& rho = trho();

    DiNum = max
    (
        solidRegionDiffNo
        (
            solidRegions[i],
            runTime,
            rho*cp,
            magKappa()
        ),
        DiNum
    );

}
