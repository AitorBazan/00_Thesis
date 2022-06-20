fvMesh& mesh = solidRegions[i];
solidThermo& therms = thermos[i];

// volScalarField& rho = rhos[i];
// volScalarField& cp = cps[i];
// volScalarField& rhoCps = rhoCpSolid[i];
// volScalarField& K = Ks[i];
// volScalarField& T = Ts[i];

tmp<volScalarField> trho = therms.rho();
const volScalarField& rhosol = trho();

tmp<volScalarField> tcp = therms.Cp();
const volScalarField& cpsol = tcp();

// tmp<volScalarField> tkappa = thermo.kappa();
// const volScalarField& kappa = tkappa();

// volScalarField& Tsol = Ts[i];

tmp<volSymmTensorField> taniAlpha;
if (!therms.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCp = therms.Kappa()/cpsol;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.transformPrincipal
        (
            mesh.cellCentres(),
            tkappaByCp()
        );

    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];
volScalarField& h = therms.he();
