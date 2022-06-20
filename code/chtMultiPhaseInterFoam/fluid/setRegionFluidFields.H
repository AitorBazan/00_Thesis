    fvMesh& mesh = fluidRegions[i];

    multiphaseSystem& thermol = thermoFluid[i];
    volScalarField& kappa = kappaLK[i];
    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& T = TFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& Cp = CpFluid[i];
    volScalarField& rhoCp = rhoCpFluid[i]; 
    surfaceScalarField& rhoPhi = rhoPhiFluid[i];
    volScalarField& rhok = rhokFluid[i];
    dimensionedScalar& rhoR = rhoRFluid[i];
    CompressibleTurbulenceModel<multiphaseSystem>& turbulence = turbulenceFluid[i];
    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];
    volScalarField& p = pFluid[i];
    volScalarField& p_rgh = p_rghFluid[i];

    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    fv::options& fvOptions = fluidFvOptions[i];

    bool frozenFlow = frozenFlowFluid[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    const dimensionedScalar rhoMax = rhoMaxFluid[i];
    const dimensionedScalar rhoMin = rhoMinFluid[i];

    const pressureControl& pressureControl = pressureControls[i];

    pimpleControl& pimple = pimpleFluid[i];

    bool correctPhi = correctPhiFluid[i];
    bool ddtCorr = ddtCorrFluid[i];
    
