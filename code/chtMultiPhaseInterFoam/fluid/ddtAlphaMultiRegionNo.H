    scalar ddtAlphaNum = -GREAT;
    forAll(fluidRegions, regioni)
    {
        ddtAlphaNum = max
        (
            ddtAlphaNo
            (
                fluidRegions[regioni],
                runTime,
                thermoFluid[regioni],
                phiFluid[regioni]
            ),
            ddtAlphaNum
        );
    }   
