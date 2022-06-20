PtrList<surfaceScalarField> phiFluid(fluidRegions.size());

forAll(fluidRegions, i)
{    
    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].phi()
        )
    );
}