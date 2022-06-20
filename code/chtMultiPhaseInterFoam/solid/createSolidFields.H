    // Initialise solid field pointer lists
    PtrList<coordinateSystem> coordinates(solidRegions.size());
    PtrList<solidThermo> thermos(solidRegions.size());
    PtrList<radiation::radiationModel> radiations(solidRegions.size());
    PtrList<fv::options> solidHeatSources(solidRegions.size());
    // PtrList<volScalarField> rhos(solidRegions.size());
    // PtrList<volScalarField> cps(solidRegions.size());
    // PtrList<volScalarField> rhoCpSolid(solidRegions.size());
    // PtrList<volScalarField> Ks(solidRegions.size());
    // PtrList<volScalarField> Ts(solidRegions.size());
    PtrList<volScalarField> betavSolid(solidRegions.size());
    PtrList<volSymmTensorField> aniAlphas(solidRegions.size());

    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;
            
        // Info<< "    Adding to Ts\n" << endl;
        // Ts.set
        // (
        //     i,
        //     new volScalarField
        //     (
        //         IOobject
        //         (
        //             "T",
        //             runTime.timeName(),
        //             solidRegions[i],
        //             IOobject::MUST_READ,
        //             IOobject::AUTO_WRITE
        //         ),
        //         solidRegions[i]
        //     )
        // );      
        Info<< "    Adding to thermos\n" << endl;
        thermos.set(i, solidThermo::New(solidRegions[i]));

        Info<< "    Adding to radiations\n" << endl;
        radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

        Info<< "    Adding fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::options(solidRegions[i])
        );

        if (!thermos[i].isotropic())
        {
            Info<< "    Adding coordinateSystems\n" << endl;
            coordinates.set
            (
                i,
                coordinateSystem::New
                (
                    solidRegions[i],
                    thermos[i],
                    coordinateSystem::typeName_()
                )
            );

            tmp<volVectorField> tkappaByCp =
                thermos[i].Kappa()/thermos[i].Cp();

            aniAlphas.set
            (
                i,
                new volSymmTensorField
                (
                    IOobject
                    (
                        "Anialpha",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedSymmTensor(tkappaByCp().dimensions(), Zero),
                    zeroGradientFvPatchSymmTensorField::typeName
                )
            );

            aniAlphas[i].primitiveFieldRef() =
                coordinates[i].transformPrincipal
                (
                    solidRegions[i].cellCentres(),
                    tkappaByCp()
                );
            aniAlphas[i].correctBoundaryConditions();

        }

        IOobject betavSolidIO
        (
            "betavSolid",
            runTime.timeName(),
            solidRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (betavSolidIO.typeHeaderOk<volScalarField>(true))
        {
            betavSolid.set
            (
                i,
                new volScalarField(betavSolidIO, solidRegions[i])
            );
        }
        else
        {
            betavSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betavSolid",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("1", dimless, scalar(1))
                )
            );
        }        
        // Info<< "    Adding to rhos\n" << endl;
        // rhos.set
        // (
        //     i,
        //     new volScalarField
        //     (
        //         IOobject
        //         (
        //             "rho",
        //             runTime.timeName(),
        //             solidRegions[i],
        //             IOobject::NO_READ,
        //             IOobject::AUTO_WRITE
        //         ),
        //         thermos[i].rho()
        //     )
        // );

        // Info<< "    Adding to Cps\n" << endl;
        // cps.set
        // (
        //     i,
        //     new volScalarField
        //     (
        //         IOobject
        //         (
        //             "cp",
        //             runTime.timeName(),
        //             solidRegions[i],
        //             IOobject::NO_READ,
        //             IOobject::AUTO_WRITE
        //         ),
        //         thermos[i].Cp()
        //     )
        // );

        // Info<< "Calculating field rhoCps\n" << endl;
        // rhoCpSolid.set
        // (
        //     i,
        //     new volScalarField
        //     (
        //         IOobject
        //         (
        //             "rhoCp",
        //             runTime.timeName(),
        //             solidRegions[i],
        //             IOobject::NO_READ,
        //             IOobject::NO_WRITE
        //         ),
        //         thermos[i].rho()*thermos[i].Cp()
        //     )
        // );

        // Info<< "    Adding to Ks\n" << endl;
        // Ks.set
        // (
        //     i,
        //     new volScalarField
        //     (
        //         IOobject
        //         (
        //             "K",
        //             runTime.timeName(),
        //             solidRegions[i],
        //             IOobject::NO_READ,
        //             IOobject::AUTO_WRITE
        //         ),
        //         thermos[i].kappa()
        //     )
        // );

  


    }
